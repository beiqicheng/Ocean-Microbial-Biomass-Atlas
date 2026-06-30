import io
from pathlib import Path
from typing import List

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


st.set_page_config(page_title="OMBA - SSU Explorer", layout="wide")

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
CONFIG_PATH = Path("config/marine_groups.yaml")


@st.cache_data(show_spinner=False)
def load_marine_groups() -> dict:
    if not CONFIG_PATH.exists():
        return {
            "Prokaryotes": ["Bacteria", "Archaea"],
            "Cyanobacteria": ["Cyanobacteria"],
            "Chlorophyta": ["Chlorophyta"],
            "Dinoflagellata": ["Dinoflagellata"],
            "Diatoms": ["Bacillariophyta", "Diatomea"],
            "Pelagophyceae": ["Pelagophyceae"],
            "Haptophyta": ["Haptophyta"],
            "Ciliophora": ["Ciliophora"],
            "Ochrophyta": ["Ochrophyta"],
            "Metazoa": ["Metazoa"],
        }
    import yaml

    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data


@st.cache_data(show_spinner=False)
def parse_tsv(uploaded_file) -> pd.DataFrame:
    if uploaded_file is None:
        return pd.DataFrame()

    text = uploaded_file.getvalue().decode("utf-8", errors="replace")
    df = pd.read_csv(io.StringIO(text), sep="	")

    if df.shape[1] < 3:
        raise ValueError("The SSU matrix must have at least three columns: Taxonomy, Length, and one sample column.")

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
def split_taxonomy(taxonomy: str) -> List[str]:
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
def weighted_mean_all(df: pd.DataFrame, sample_cols: List[str]):
    counts = df[sample_cols].sum(axis=1)
    total = int(counts.sum())
    if total == 0:
        return float("nan"), 0
    return float((df["Length"] * counts).sum() / total), total


@st.cache_data(show_spinner=False)
def summarize_rows(df: pd.DataFrame, sample_cols: List[str]):
    counts = df[sample_cols].sum(axis=1) if not df.empty else pd.Series(dtype=float)
    if df.empty:
        return {
            "refs": 0,
            "weighted": float("nan"),
            "mean": float("nan"),
            "minimum": float("nan"),
            "maximum": float("nan"),
            "total_counts": 0,
        }
    total = int(counts.sum())
    return {
        "refs": int(len(df)),
        "weighted": float((df["Length"] * counts).sum() / total) if total else float("nan"),
        "mean": float(df["Length"].mean()),
        "minimum": float(df["Length"].min()),
        "maximum": float(df["Length"].max()),
        "total_counts": total,
    }


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
def build_plot_frame(df: pd.DataFrame, sample_cols: List[str]) -> pd.DataFrame:
    rows = []
    for sample in sample_cols:
        total = int(df[sample].sum())
        weighted = weighted_mean(df, sample) if total else float("nan")
        rows.append({"Sample": sample, "Weighted length (bp)": weighted, "Total counts": total})
    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False)
def build_rank_tree(df: pd.DataFrame, rank: str):
    if df.empty:
        return pd.DataFrame(columns=["Group", "References", "Total counts", "Weighted length (all samples)"])
    grouped = []
    for gname, gdf in df.groupby(df["Taxonomy"].apply(lambda x: get_rank_label(x, rank))):
        w, total = weighted_mean_all(gdf, list(df.columns[2:]))
        grouped.append({
            "Group": gname,
            "References": int(len(gdf)),
            "Total counts": total,
            "Weighted length (all samples)": w,
        })
    return pd.DataFrame(grouped).sort_values("Total counts", ascending=False)


@st.cache_data(show_spinner=False)
def make_figure(plot_df: pd.DataFrame, query: str, color: str = "#1f4e78"):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=plot_df["Sample"],
            y=plot_df["Weighted length (bp)"],
            mode="lines+markers",
            line=dict(color=color, width=3),
            marker=dict(size=8, color=color),
            name="Weighted length",
            hovertemplate="Sample=%{x}<br>Weighted length=%{y:.1f} bp<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Weighted SSU length for {query}",
        template="plotly_white",
        height=520,
        margin=dict(l=20, r=20, t=60, b=120),
        hovermode="x unified",
        xaxis_title="Sample",
        yaxis_title="Weighted length (bp)",
    )
    fig.update_xaxes(tickangle=-45)
    return fig


st.title("SSU Explorer")
st.caption("Explore weighted SSU length by taxon or rank. Use the marine-groups shortcuts for faster navigation.")

uploaded = st.file_uploader("Upload SSU_matrix.txt", type=["txt", "tsv", "csv"])

if uploaded is None:
    st.info("Upload SSU_matrix.txt to begin.")
    st.stop()

try:
    df = parse_tsv(uploaded)
except Exception as e:
    st.error(str(e))
    st.stop()

sample_cols = list(df.columns[2:])
if not sample_cols:
    st.error("No sample columns were found in the SSU matrix.")
    st.stop()

marine_groups = load_marine_groups()

# Session defaults
if "ssu_mode" not in st.session_state:
    st.session_state.ssu_mode = "taxon"
if "ssu_rank" not in st.session_state:
    st.session_state.ssu_rank = "phylum"
if "ssu_query" not in st.session_state:
    st.session_state.ssu_query = "Dinoflagellata"

left, right = st.columns([1.25, 1.05], gap="large")

with left:
    st.subheader("Taxonomy browser")
    st.markdown("#### Marine groups")

    cols = st.columns(2)
    group_names = list(marine_groups.keys())
    for i, group_name in enumerate(group_names):
        col = cols[i % 2]
        if col.button(group_name, use_container_width=True):
            st.session_state.ssu_query = marine_groups[group_name][0]
            st.session_state.ssu_mode = "taxon"

    st.markdown("#### Full taxonomy")
    st.session_state.ssu_mode = st.radio(
        "Search mode",
        ["taxon", "rank"],
        horizontal=True,
        index=0 if st.session_state.ssu_mode == "taxon" else 1,
        key="ssu_mode_radio",
        label_visibility="collapsed",
    )
    if st.session_state.ssu_mode_radio == "taxon":
        query = st.selectbox(
            "Choose a taxon",
            options=sorted(pd.unique(df["Taxonomy"].apply(lambda x: get_rank_label(x, "phylum")))),
            index=0,
            help="Pick a higher-level taxon. You can later replace this with a clickable tree.",
        )
        st.session_state.ssu_query = query
    else:
        st.session_state.ssu_rank = st.selectbox(
            "Taxonomic rank",
            RANKS,
            index=RANKS.index(st.session_state.ssu_rank) if st.session_state.ssu_rank in RANKS else 0,
            key="ssu_rank_select",
        )
        query = st.text_input(
            "Rank term",
            value=st.session_state.ssu_query,
            help="Examples: Alveolata, Bacteria, Dinoflagellata, Metazoa.",
        )
        st.session_state.ssu_query = query

    filtered = filter_rows(df, st.session_state.ssu_query, st.session_state.ssu_mode, st.session_state.ssu_rank)
    plot_df = build_plot_frame(filtered, sample_cols)
    fig = make_figure(plot_df, st.session_state.ssu_query, color="#1f4e78")
    st.plotly_chart(fig, use_container_width=True)

    fig_bytes = None
    try:
        fig_bytes = fig.to_image(format="png", scale=2)
    except Exception:
        fig_bytes = fig.to_html(include_plotlyjs="cdn").encode("utf-8")

    st.download_button(
        "Download figure",
        data=fig_bytes,
        file_name=f"SSU_{st.session_state.ssu_query.replace(' ', '_')}.png" if isinstance(fig_bytes, (bytes, bytearray)) else f"SSU_{st.session_state.ssu_query.replace(' ', '_')}.html",
        mime="image/png" if isinstance(fig_bytes, (bytes, bytearray)) and fig_bytes[:8] == b"PNG

" else "text/html",
    )

with right:
    st.subheader("Summary statistics")
    summary = summarize_rows(filtered, sample_cols)
    m1, m2 = st.columns(2)
    m1.metric("Matched references", f"{summary['refs']:,}")
    m2.metric("Total counts", f"{summary['total_counts']:,}")
    m3, m4 = st.columns(2)
    m3.metric("Weighted length", f"{summary['weighted']:,.1f} bp" if pd.notna(summary["weighted"]) else "—")
    m4.metric("Mean length", f"{summary['mean']:,.1f} bp" if pd.notna(summary["mean"]) else "—")
    m5, m6 = st.columns(2)
    m5.metric("Min length", f"{summary['minimum']:,.0f} bp" if pd.notna(summary["minimum"]) else "—")
    m6.metric("Max length", f"{summary['maximum']:,.0f} bp" if pd.notna(summary["maximum"]) else "—")

    st.markdown("### Top matching rows")
    display_cols = ["Taxonomy", "Length"] + sample_cols
    if st.session_state.ssu_mode == "rank":
        filtered = filtered.copy()
        filtered[st.session_state.ssu_rank] = filtered["Taxonomy"].apply(lambda x: get_rank_label(x, st.session_state.ssu_rank))
        display_cols = ["Taxonomy", st.session_state.ssu_rank, "Length"] + sample_cols
    st.dataframe(filtered[display_cols].head(200), use_container_width=True, height=360)

st.divider()

rank_df = build_rank_tree(filtered, st.session_state.ssu_rank)
if st.session_state.ssu_mode == "rank" and not rank_df.empty:
    st.subheader(f"{st.session_state.ssu_rank.title()} groups")
    st.dataframe(rank_df.head(20), use_container_width=True, height=260)
    rank_fig = go.Figure()
    top = rank_df.head(15)
    rank_fig.add_trace(
        go.Bar(
            x=top["Group"],
            y=top["Weighted length (all samples)"],
            marker_color="#4c78a8",
            hovertemplate="Group=%{x}<br>Weighted length=%{y:.1f} bp<extra></extra>",
        )
    )
    rank_fig.update_layout(
        title=f"Top {st.session_state.ssu_rank} groups by weighted SSU length",
        template="plotly_white",
        height=420,
        margin=dict(l=20, r=20, t=60, b=120),
        xaxis_title=st.session_state.ssu_rank.title(),
        yaxis_title="Weighted length (bp)",
    )
    rank_fig.update_xaxes(tickangle=-45)
    st.plotly_chart(rank_fig, use_container_width=True)

st.markdown("### Notes")
st.write(
    "Weighted SSU length is calculated as sum(length × counts) / sum(counts). "
    "The Marine Groups buttons are shortcuts for the dominant eukaryotic and prokaryotic groups you care about most."
)
