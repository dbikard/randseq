"""Compatibility shim for code written against `randseq.core_v2`.

`core_v2` was merged into `core` (commit b758f8f). This module re-exports the merged
implementations so that existing scripts — notably `code/fig3/generate_fig3_motifs.py` in the
randseq_paper repo, which does `from randseq.core_v2 import find_restricted_motifs_mp` — keep
importing successfully.

The functions you get from here are the ones in `randseq.core`, not the old `core_v2` ones, so
results will differ from a pre-merge run. See CHANGELOG.md for what moved and why. New code
should import from `randseq.core` directly; this module will be removed in a future release.
"""

import warnings as _warnings

from .core import *  # noqa: F401,F403
from .core import (  # noqa: F401
    find_restricted_motifs,
    find_restricted_motifs_mp,
    TooManyCandidatesError,
)

_warnings.warn(
    "randseq.core_v2 is deprecated: core_v2 was merged into randseq.core. This shim re-exports "
    "the merged implementations, whose results differ from the old core_v2 (see CHANGELOG.md). "
    "Import from randseq.core instead.",
    DeprecationWarning,
    stacklevel=2,
)
