"""Weighted length calculations."""
from __future__ import annotations
import pandas as pd

def weighted_mean(df: pd.DataFrame, sample_col: str) -> float:
    total = df[sample_col].sum()
    if total == 0:
        return float("nan")
    return float((df["Length"] * df[sample_col]).sum() / total)
