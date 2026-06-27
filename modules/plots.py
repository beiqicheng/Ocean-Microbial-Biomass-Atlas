"""Plot helpers."""
from __future__ import annotations
import plotly.express as px

def line_plot(df, x_col="Sample", y_col="Weighted length (bp)", title=""):
    fig = px.line(df, x=x_col, y=y_col, markers=True, title=title)
    fig.update_layout(height=480, margin=dict(l=20, r=20, t=60, b=80))
    return fig
