"""Parse OMBA project files."""
from __future__ import annotations
import io
import pandas as pd

def parse_table(uploaded_file) -> pd.DataFrame:
    if uploaded_file is None:
        return pd.DataFrame()
    text = uploaded_file.getvalue().decode("utf-8", errors="replace")
    return pd.read_csv(io.StringIO(text), sep="\t")
