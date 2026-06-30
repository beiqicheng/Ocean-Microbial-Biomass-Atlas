import io
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st


st.set_page_config(page_title="OMBA - LSU Explorer", layout="wide")

RANKS = ["phylum", "class", "order", "family", "genus", "species"]


@st.cache_data(show_spinner=False)
def parse_tsv(uploaded_file) -> pd.DataFrame:
    """Parse a SILVA-like LSU matrix.

    Expected columns:
    Taxonomy | Length | Sample1 | Sample2 | ...
    """
    if uploaded_file is None:
        return pd.DataFrame()

    text = uploaded_file.getvalue().decode("utf-8", errors="replace")
    df = pd.read_csv(io.StringIO(text), sep="\t")

    if df.shape[1] < 3:
        raise ValueError("The LSU matrix must have at least three columns: Taxonomy, Length, and one sample column.")

    if "Taxonomy" not in df.columns or "Length" not in df.columns:
        cols = list(df.columns)
        cols[0] = "Taxonomy"
        cols[1] = "Length"
        df.columns = cols

    df["Taxonomy"] = df["Taxonomy"].astype(str)
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce").fillna(0).astype(int)
    for c in df.columns[2:]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    return df


@st.cache_data(show_spinner=False)
def split_taxonomy(taxonomy: str):
    return [x.strip() for x in str(taxonomy).split(";") if x.strip()]


@st.cache_data(show_spinner=False)
def get_rank_label(taxonomy: str, rank: str) -> str:
    parts = split_taxonomy(taxonomy)
    if not parts:
        return "Unclassified"
    idx_map = {"phylum": 1, "class": 2, "order": 3, "family": 4, "genus": 5, "species": 6}
    idx = idx_map.get(rank)
    if idx is not None and idx < len(parts):
        return parts[idx]
    return parts[-1]


@st.cache_data(show_spinner=False)
def weighted_mean(df: pd.DataFrame, sample_col: str) -> float:
    total = df[sample_col].sum()
    if total == 0:
        return float("nan")
    return float((df["Length"] * df[sample_col]).sum() / total)


@st.cache_data(show_spinner=False)
def weighted_mean_all(df: pd.DataFrame, sample_cols: list[str]):
    counts = df[sample_cols].sum(axis=1)
    total = int(counts.sum())
    if total == 0:
        return float("nan"), 0
    return float((df["Length"] * counts).sum() / total), total


@st.cache_data(show_spinner=False)
def filter_rows(df: pd.DataFrame, query: str, mode: str, rank: str) -> pd.DataFrame:
    q = query.strip().lower()
    if not q or df.empty:
        return df.iloc[0:0].copy()
    if mode == "taxon":
        mask = df["Taxonomy"].str.lower().str.contains(q, na=False)
    else:
        mask = df["Taxonomy"].apply(lambda x: q in get_rank_label(x, rank).lower())
    return df.loc[mask].copy()


@st.cache_data(show_spinner=False)
def build_plot_frame(df: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    rows = []
    for sample in sample_cols:
        total = int(df[sample].sum())
        weighted = weighted_mean(df, sample) if total else float("nan")
        rows.append({"Sample": sample, "Weighted length (bp)": weighted, "Total counts": total})
    return pd.DataFrame(rows)


st.title("LSU Explorer")
st.caption("Upload an LSU matrix, choose a taxon or rank, and inspect the weighted LSU length across samples.")

with st.sidebar:
    st.header("Input")
    uploaded = st.file_uploader("Upload LSU_matrix.txt", type=["txt", "tsv", "csv"])
    st.divider()
    st.subheader("Controls")
    mode = st.radio("Search mode", ["taxon", "rank"], horizontal=True)
    rank = st.selectbox("Taxonomic rank", RANKS, index=0, disabled=(mode != "rank"))
    query = st.text_input("Taxon or rank term", value=st.session_state.get("lsu_query", "Dinoflagellata"))

if uploaded is None:
    st.info("Upload LSU_matrix.txt to begin.")
    st.stop()

try:
    df = parse_tsv(uploaded)
except Exception as e:
    st.error(str(e))
    st.stop()

sample_cols = list(df.columns[2:])
if not sample_cols:
    st.error("No sample columns were found in the LSU matrix.")
    st.stop()

st.session_state["lsu_query"] = query
filtered = filter_rows(df, query, mode, rank)
summary_weighted, summary_total = weighted_mean_all(filtered, sample_cols)
plot_df = build_plot_frame(filtered, sample_cols)

left, right = st.columns([1.4, 1])

with left:
    st.subheader("Weighted LSU length across samples")
    if filtered.empty:
        st.warning("No matching taxa. Try a broader term such as Metazoa, Dinoflagellata, Alveolata, or Bacteria.")
    else:
        fig = px.line(
            plot_df,
            x="Sample",
            y="Weighted length (bp)",
            markers=True,
            title=f"Weighted LSU length for {query}",
        )
        fig.update_layout(height=520, margin=dict(l=20, r=20, t=60, b=120))
        fig.update_xaxes(tickangle=-45)
        st.plotly_chart(fig, use_container_width=True)

    st.download_button(
        "Download plot data as CSV",
        data=plot_df.to_csv(index=False).encode("utf-8"),
        file_name=f"LSU_{query.replace(' ', '_')}_weighted_length.csv",
        mime="text/csv",
        disabled=filtered.empty,
    )

with right:
    st.subheader("Summary")
    c1, c2 = st.columns(2)
    c1.metric("Matched rows", f"{len(filtered):,}")
    c2.metric("Total counts", f"{summary_total:,}")
    st.metric("Weighted length (all samples)", f"{summary_weighted:,.1f} bp" if pd.notna(summary_weighted) else "—")

    st.markdown("### Top matching rows")
    display_cols = ["Taxonomy", "Length"] + sample_cols
    if mode == "rank":
        filtered = filtered.copy()
        filtered[rank] = filtered["Taxonomy"].apply(lambda x: get_rank_label(x, rank))
        display_cols = ["Taxonomy", rank, "Length"] + sample_cols
    st.dataframe(filtered[display_cols].head(200), use_container_width=True, height=340)

st.divider()
st.markdown("### Notes")
st.write(
    "This page calculates weighted LSU length using: sum(length × counts) / sum(counts). "
    "Later we can connect this page to saved projects, taxonomy-tree clicks, and RNA/biomass pages."
)
