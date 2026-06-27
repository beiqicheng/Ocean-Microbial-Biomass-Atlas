"""Taxonomy utilities."""
from __future__ import annotations

def split_taxonomy(taxonomy: str):
    return [x.strip() for x in str(taxonomy).split(";") if x.strip()]

def get_rank_label(taxonomy: str, rank: str) -> str:
    parts = split_taxonomy(taxonomy)
    idx_map = {"phylum": 1, "class": 2, "order": 3, "family": 4, "genus": 5, "species": 6}
    idx = idx_map.get(rank)
    if not parts:
        return "Unclassified"
    if idx is not None and idx < len(parts):
        return parts[idx]
    return parts[-1]
