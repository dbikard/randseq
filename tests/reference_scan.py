"""The readable reference implementation of the pattern scan.

This is the original, string-based scan: for every sequence it builds a Python set of the motif
strings matching a pattern, then aggregates fold changes per motif. `randseq.core` computes the
same table with integer codes and `np.bincount`, roughly a hundred times faster.

It lives here, outside the package, because it is not part of the API -- nothing should import
it to do analysis. It is kept because it is the only *independent* check that vectorising the
scan changed no number: it is written to be obviously correct rather than fast, so when the two
disagree, the fast one is what is wrong. `test_reference_scan.py` is what uses it.

Extracted from randseq/core.py, unchanged, when the package stopped exporting it.
"""
import re
from collections import defaultdict
from typing import Callable, Dict, List, Set, Tuple

import numpy as np
import pandas as pd

from randseq.core import score
from randseq.utils import revcomp

def _make_regex_group_str(length: int, content_char_class: str = "[ATGC]") -> str:
    """
    Creates a regex capturing group string.
    If length is 0, creates an empty capturing group "()".
    Otherwise, creates a group like "([ATGC]{length})".
    """
    if length > 0:
        return f"({content_char_class}{{{length}}})"
    else: # length == 0
        return "()"


def get_sites_in_seq(
    sequence_library: List[str],
    pattern: Tuple[int, int, int] = (3, 6, 4),
    no_ori: bool = True,
    custom_revcomp: Callable[[str], str] = revcomp
) -> List[Set[str]]:
    """
    Computes the list of unique sites matching a defined pattern within each sequence
    of a given library.

    The pattern specifies a motif structure: (defined_bases_part1, spacer_length, defined_bases_part2).
    For example, a pattern (3, 6, 4) looks for sites like XXXNNNNNNXXXX, where X
    represents a defined base (A, T, G, C) and N represents an undefined base in the motif.

    Args:
        sequence_library: A list of DNA sequence strings. They do NOT have to be the same
                          length -- real libraries are ragged (in the published screen only
                          ~91% of inserts were exactly the nominal 150 nt). The pattern is
                          only required to fit inside the shortest one.
        pattern: A tuple of three non-negative integers:
            1. The number of defined bases in the first part of the site.
            2. The length of the spacer (undefined bases) in the middle.
            3. The number of defined bases in the second part of the site.
        no_ori: If True (default), the search for sites is performed on both the
                provided sequence and its reverse complement. If False, only the
                provided sequence strand is searched.
        custom_revcomp: A function to compute the reverse complement. Defaults to
                        a basic internal `revcomp` function.

    Returns:
        A list of sets. Each set contains unique site strings found in the
        corresponding input sequence. Sites are represented with 'N' for the
        spacer region (e.g., "ATCNNNNGGC").

    Raises:
        ValueError: If the sequence library is empty.
        ValueError: If the pattern is not a tuple of three non-negative integers.
        ValueError: If the total length of the site defined by the pattern is zero.
        ValueError: If the total length of the site pattern exceeds the length of the
                    sequences in the library.
    """
    # --- Input Validation ---
    if not sequence_library:
        raise ValueError("Sequence library cannot be empty.")

    if not (
        isinstance(pattern, tuple) and
        len(pattern) == 3 and
        all(isinstance(p_num, int) and p_num >= 0 for p_num in pattern)
    ):
        raise ValueError(
            "Pattern must be a tuple of three non-negative integers: "
            "(defined_bases_part1, spacer_length, defined_bases_part2)."
        )

    defined_len1, spacer_len, defined_len2 = pattern
    pattern_total_length = defined_len1 + spacer_len + defined_len2

    if pattern_total_length == 0:
        # Finding motifs of total length 0 (e.g., pattern (0,0,0)) is problematic
        # as it would match everywhere. It's better to disallow it.
        raise ValueError("The total length of the site pattern (sum of its parts) cannot be zero.")

    shorter_seq = min([len(s) for s in sequence_library])
    if pattern_total_length > shorter_seq:
        raise ValueError(
            f"Site pattern's total length ({pattern_total_length}) is longer than "
            f"the shorter sequence length ({shorter_seq})."
        )

    # --- Regex Preparation ---
    # Uses the helper _make_regex_group_str to create "()" for zero-length parts.
    # This ensures three capturing groups are always present for consistent unpacking later.
    group1_str = _make_regex_group_str(defined_len1)
    # The spacer part in the original regex also used [ATGC] for the characters
    # it matched in the sequence, even though these are replaced by 'N's in the motif.
    group2_str = _make_regex_group_str(spacer_len) 
    group3_str = _make_regex_group_str(defined_len2)
    
    # The lookahead `(?=...)` is crucial for finding all overlapping sites.
    regex_pattern_str = f"(?={group1_str}{group2_str}{group3_str})"
    
    # Using re.IGNORECASE for robustness, assuming 'atgc' is equivalent to 'ATGC'.
    compiled_regex = re.compile(regex_pattern_str, re.IGNORECASE)

    # --- Site Finding ---
    all_sequence_results: List[Set[str]] = []
    
    for seq_original_case in sequence_library:
        current_sequence_sites: Set[str] = set()
        
        # Standardize sequence to uppercase for consistent processing and regex matching
        processed_seq = seq_original_case.upper()

        sequences_to_scan = [processed_seq]
        if no_ori:
            # Ensure the revcomp function also standardizes or expects uppercase
            reverse_complement_seq = custom_revcomp(processed_seq)
            sequences_to_scan.append(reverse_complement_seq)

        for target_sequence in sequences_to_scan:
            # Find all matches in the current target sequence (forward or reverse complement)
            # matches will be a list of 3-tuples due to the three capturing groups.
            matches = compiled_regex.findall(target_sequence)
            for part1, actual_spacer_bases, part2 in matches: # Unpack the 3 captured groups
                # Construct the site motif string using 'N's for the spacer region.
                # The `actual_spacer_bases` captured by the regex are intentionally
                # ignored here for the motif string, which uses a generic 'N' spacer.
                # If spacer_len is 0, 'N' * 0 is an empty string.
                site_motif = f"{part1}{'N' * spacer_len}{part2}"
                current_sequence_sites.add(site_motif)
        
        all_sequence_results.append(current_sequence_sites)

    return all_sequence_results


def get_fold_change_values_per_site(site_sets_list, fold_changes_list):
    """
    Computes a dictionary mapping each unique site to a list of its associated
    log2 Fold Change (FC) values using collections.defaultdict for optimization.

    Args:
        site_sets_list (list): A list of sets, where each set contains the 
                               sequence motifs found in a member of the library.
        fold_changes_list (list): A list of log2FC values, corresponding to each
                                  sequence in the library.

    Returns:
        dict: A dictionary where keys are sequence motifs and values are numpy arrays
              containing the log2FC values associated with that site.
    """
    # Step 1: Initialize a defaultdict to store lists of FC values for each site.
    # defaultdict(list) will automatically create a new list for a site
    # the first time it's encountered.
    site_to_fcs_map_intermediate = defaultdict(list)

    # Step 2: Iterate through the site sets and their corresponding fold changes.
    # zip() pairs each site_set with its fc_value.
    for site_set, fc_value in zip(site_sets_list, fold_changes_list):
        # For each site in the current set, append the fc_value
        # to the list associated with that site in the dictionary.
        for site in site_set:
            site_to_fcs_map_intermediate[site].append(fc_value)

    # Step 3: Convert the lists of FC values to NumPy arrays.
    # This is done after all values are collected to efficiently create the arrays.
    site_to_fcs_map_final = {
        site: np.array(fcs_list)
        for site, fcs_list in site_to_fcs_map_intermediate.items()
    }

    return site_to_fcs_map_final


def get_sites_scores(site_FCs, pattern, log2FC_thr=-1):
    """ 
    Scores sites based on depletion and occurrence.
    Returns a pandas DataFrame with 'site', 'fraction_depleted', 'num_sequences', 'avg_log2fc'.
    """
    site_scores_list = []
    
    for site_key, fc_values in site_FCs.items():
        if not isinstance(fc_values, (list, np.ndarray)) or len(fc_values) == 0:
            continue # Skip if fc_values is not a list/array or is empty
            
        current_score = score(fc_values, log2FC_thr)
        num_occurrences = len(fc_values)
        avg_fc = np.mean(fc_values) if num_occurrences > 0 else np.nan
        
        site_scores_list.append({
            'motif': site_key,
            'fraction_depleted': current_score,
            'num_sequences': num_occurrences,
            'avg_log2fc': avg_fc
        })
            
    if not site_scores_list:
        return pd.DataFrame(columns=['motif', 'fraction_depleted', 'num_sequences', 'avg_log2fc'])
    return pd.DataFrame(site_scores_list)
