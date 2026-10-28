"""Assign VALENCIA CSTs to every swab.

Prevotella, Atopobium and Sneathia species are merged to genus in the centroids (V4 cannot
separate them), and Valencia.py is run with a one-line fix for current pandas
(positional indexing via .iloc).
In:  data/processed/valencia/valencia_input.csv, tools/VALENCIA
Out: data/processed/valencia/{CST_centroids_merged, valencia_output}.csv
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[1]
vdir = root / "tools" / "VALENCIA"
d11 = root / "data" / "processed" / "valencia"

# Centroids with species merged to genus
c = pd.read_csv(vdir / "CST_centroids_012920.csv")
for genus in ["Prevotella", "Atopobium", "Sneathia"]:
    cols = [x for x in c.columns if x.startswith(genus + "_")]
    gcol = "g_" + genus
    c[gcol] = c[cols].sum(axis=1) + (c[gcol] if gcol in c.columns else 0)
    c = c.drop(columns=cols)
cent_merged = d11 / "CST_centroids_merged.csv"
c.to_csv(cent_merged, index=False)

# Patched copy of Valencia.py
src = (vdir / "Valencia.py").read_text()
assert src.count("median[taxon_count]") == 2
with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
    f.write(src.replace("median[taxon_count]", "median.iloc[taxon_count]"))
    patched = f.name

subprocess.run([sys.executable, patched, "-ref", str(cent_merged),
                "-i", str(d11 / "valencia_input.csv"), "-o", str(d11 / "valencia_output")], check=True)
Path(patched).unlink()
