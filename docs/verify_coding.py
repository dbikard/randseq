"""Empirical check of the guarantees the motif coding relies on."""
import numpy as np, itertools, random
from randseq.core import (encode_library, _pattern_window_codes, decode_motif_codes,
                          get_pattern_scores, get_sites_in_seq,
                          get_fold_change_values_per_site, get_sites_scores)

rng = random.Random(0)

# --- 1. encode/decode is a bijection on the whole code space -----------------
print("1. round-trip over the ENTIRE code space")
for pat in [(6,0,0),(3,4,4),(2,1,3),(1,0,1),(4,8,4)]:
    n = 4**(pat[0]+pat[2])
    strs = decode_motif_codes(np.arange(n), pat)
    assert len(set(strs)) == n, f"{pat}: not injective"
    # re-encode by hand from the string and check we get the code back
    B = {b:i for i,b in enumerate("ACGT")}
    d1,sp,d2 = pat
    for c in rng.sample(range(n), min(n, 500)):
        s = strs[c]
        assert len(s) == d1+sp+d2 and s[d1:d1+sp] == "N"*sp, s
        defined = s[:d1] + s[d1+sp:]
        back = 0
        for ch in defined: back = back*4 + B[ch]
        assert back == c, (pat, c, s, back)
    print(f"   {pat}: {n} codes, all distinct, re-encode exact")

# --- 2. the rejection sentinel cannot collide with a real motif --------------
print("2. sentinel is exactly one past the top code")
for pat in [(6,0,0),(3,4,4),(4,8,4)]:
    n = 4**(pat[0]+pat[2])
    fwd,_ = encode_library(["N"*(pat[0]+pat[1]+pat[2])])
    assert _pattern_window_codes(fwd, pat)[0,0] == n
    assert decode_motif_codes(np.arange(n), pat)  # every code below it is a real motif
    print(f"   {pat}: sentinel {n}, real codes 0..{n-1}")

# --- 3. reverse complement is an involution and matches string revcomp -------
print("3. revcomp")
COMP = str.maketrans("ACGT", "TGCA")
B = {b: i for i, b in enumerate("ACGT")}
for _ in range(300):
    L = rng.randint(5, 60)
    s_ = "".join(rng.choice("ACGTN" if rng.random() < .3 else "ACGT") for _ in range(L))
    fwd, rc = encode_library([s_])
    want = s_.translate(COMP)[::-1]                      # N is not in the table, stays N
    assert [B.get(ch, 4) for ch in s_] == fwd[0, :L].tolist()
    assert [B.get(ch, 4) for ch in want] == rc[0, :L].tolist(), (s_, want)
    assert encode_library([want])[1][0, :L].tolist() == fwd[0, :L].tolist()   # involution
print("   300 sequences (incl. N): rc == string revcomp, and rc(rc) == fwd")

# --- 4. equivalence with the regex reference, on nasty input ----------------
print("4. vectorized == regex reference, on ragged/N-containing libraries")
tot = 0
for trial in range(40):
    seqs, fcs = [], []
    for _ in range(rng.randint(3, 25)):
        L = rng.randint(16, 45)   # reference impl requires len >= pattern width (max 13)
        alpha = "ACGT" if rng.random() < .7 else "ACGTN"
        seqs.append("".join(rng.choice(alpha) for _ in range(L)))
        fcs.append(rng.uniform(-6, 2))
    pat = rng.choice([(6,0,0),(3,4,4),(2,1,3),(3,6,4),(1,2,2),(4,4,3)])
    ref = get_sites_scores(get_fold_change_values_per_site(
              get_sites_in_seq(seqs, pattern=pat, no_ori=True), fcs), pat).set_index('motif').sort_index()
    fast = get_pattern_scores(encode_library(seqs), fcs, pat).set_index('motif').sort_index()
    assert fast.index.tolist() == ref.index.tolist(), (pat, trial)
    assert fast.num_sequences.tolist() == ref.num_sequences.tolist()
    assert np.allclose(fast.fraction_depleted, ref.fraction_depleted, atol=1e-12)
    assert np.allclose(fast.avg_log2fc, ref.avg_log2fc, atol=1e-9)
    tot += len(ref)
print(f"   40 randomized libraries, {tot} motifs, identical on every column")

# --- 5. the guarantee has a stated edge, and it fails loudly ---------------
print("5. the int32 boundary is checked, not assumed")
fwd, _ = encode_library(["ACGT" * 20])
for k in (15, 16):
    try:
        _pattern_window_codes(fwd, (k, 0, 0))
        print(f"   d1+d2={k}: accepted ({4**k:,} codes)")
    except ValueError as e:
        print(f"   d1+d2={k}: rejected -> {str(e)[:70]}...")
print("\nALL CHECKS PASSED")
