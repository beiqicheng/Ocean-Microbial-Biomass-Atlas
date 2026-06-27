"""Shared helpers."""
from __future__ import annotations

def fmt(n):
    try:
        return f"{n:,.1f}"
    except Exception:
        return "—"
