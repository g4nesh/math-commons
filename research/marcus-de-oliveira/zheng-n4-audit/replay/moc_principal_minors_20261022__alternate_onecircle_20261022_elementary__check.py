"""Replay for the supplement copy moc_principal_minors_20261022/alternate_onecircle_20261022_elementary.pdf.

This copy differs from the __workspace__ copy only in one sentence of the introduction (a
scope disclaimer); Lemma 1, Theorems 1 and 2, all displayed formulas and proofs are the same
(compared line by line on the text extracts and the page images).  The checks are therefore
the ones in the sibling script, which is own code in this directory and is run here unchanged.
"""
import os
import runpy
import sys

here = os.path.dirname(os.path.abspath(__file__))
mod = runpy.run_path(os.path.join(here, "__workspace____alternate_onecircle_20261022_elementary__check.py"))
print("document copy: moc_principal_minors_20261022__alternate_onecircle_20261022_elementary")
sys.exit(mod["run_all"]())
