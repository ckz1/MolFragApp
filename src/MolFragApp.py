r"""MolFragApp (Scientific UI)
===============

Molecular dynamics data processing - v2.2 (integrated logo and metrics header)
Version: v2.2
Date: 2026.05.22

UI updates relative to v1.2/v1.3:
- Adds a short trajectory ID for navigation while keeping the full HDF5 path in metadata.
- Adds filter reset, active-filter summary, and filtered-data export.
- Adds a global selected-trajectory status bar above the analysis tabs.
- Removes the second dataframe from the statistics tab to avoid conflicting selection logic.
- Adds event-level summary cards, improved empty/error states, and clearer metadata layout.
- Adds energy-view presets, viewer-size control, current-frame export, and CSV export buttons.
- Reorganizes Dalitz/Newton controls as momentum-space event-map controls.
- Adopts the CPC figure style: white scientific dashboard, blue/teal/purple/green accent palette, rounded cards, pill tabs, selected-record banner, workflow focus, and export/reuse panels.
- v1.7 adds schema validation, persistent selected-record identity, current-subset momentum maps, richer Notebook/CLI manifests, and table-first three-body export.
- v1.8 unifies selection/slider accent colors and adds trajectory GIF export from HDF5 geometry frames.
- v1.9 replaces the large hero panel with a compact logo/status bar for manuscript-ready screenshots.
- v2.0 removes duplicated dataset/schema/descriptor chips from the top bar; dataset checks stay in the sidebar and Metadata tab.
- v2.2 uses a larger logo block and a 2x2 metrics grid, removing the stray HTML artifact and improving header proportions.
"""

import os
import base64
import json
import io
import html
import ast
from glob import glob

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

VERSION = "v2.2 (compact logo + 2x2 metrics UI)"
st.set_page_config(page_icon="⚗️", page_title="MolFragApp", layout="wide")

LANG_OPTIONS = ["English", "中文"]
FILTER_KEYS = ["gf_state", "gf_frag", "gf_ker", "gf_minstep", "gf_nfrag"]


def get_lang() -> str:
    return st.session_state.get("ui_lang", "English")


def tr(zh: str, en: str) -> str:
    return zh if get_lang() == "中文" else en


def help_icon(help_text: str) -> str:
    safe = html.escape(help_text, quote=True)
    return f'<span class="mf-help-icon" data-tooltip="{safe}" tabindex="0">i</span>'


def normalize_matplotlib_figure(fig_like, size=(4.6, 4.6)):
    fig = getattr(fig_like, "figure", fig_like)
    try:
        fig.set_size_inches(*size, forward=True)
        fig.tight_layout()
    except Exception:
        pass
    return fig


# Try importing utils
try:
    from utils import *
except ImportError:
    st.error(
        tr(
            "无法导入 utils.py，请确保 utils.py 在当前目录下且包含更新后的代码。",
            "Cannot import utils.py. Please make sure utils.py is in the current directory and contains the updated code.",
        )
    )
    st.stop()


# ================= CPC 图件风格样式注入 =================
def inject_custom_css(publication_mode: bool = False):
    """Inject a publication-style scientific dashboard theme.

    The color palette and card layout follow the CPC UI figure set:
    clean white canvas, pale gray panels, blue active controls, and
    cyan/purple/green secondary accents.
    """
    publication_css = """
        footer {visibility: hidden;}
        .stDeployButton {display: none;}
        header[data-testid="stHeader"] {background: rgba(255,255,255,0.92);}
    """ if publication_mode else """
        header[data-testid="stHeader"] {background: rgba(255,255,255,0.86);}
    """

    st.markdown(
        f"""
        <style>
        :root {{
            --mf-blue: #6EA8D5;
            --mf-blue-2: #9BC3E2;
            --mf-cyan: #89D1DA;
            --mf-control: #5D8FBA;
            --mf-control-dark: #315F86;
            --mf-control-bg: #E7F1FA;
            --mf-purple: #B7AFDF;
            --mf-green: #93CFA7;
            --mf-orange: #EDCB89;
            --mf-red: #DA8A8A;
            --mf-ink: #2C3E50;
            --mf-muted: #6C757D;
            --mf-light: #8B95A1;
            --mf-bg: #FFFFFF;
            --mf-canvas: #F7F9FC;
            --mf-panel: #FFFFFF;
            --mf-sidebar: #F8F9FA;
            --mf-line: #E9ECEF;
            --mf-line-2: #DDE3EA;
            --mf-soft-blue: #F1F6FF;
            --mf-soft-cyan: #F1FBFD;
            --mf-soft-green: #F1FAF4;
            --mf-soft-purple: #F6F3FF;
            --mf-soft-orange: #FFF8EA;
            --mf-gray-tab: #F1F3F5;
        }}

        .stApp {{
            background: linear-gradient(180deg, #ffffff 0%, #fbfcfe 60%, #f7f9fc 100%);
            color: var(--mf-ink);
        }}

        .block-container {{
            max-width: 96vw;
            padding-top: 2.0rem;
            padding-bottom: 2.0rem;
        }}
        @media (max-width: 1400px) {{
            .block-container {{
                max-width: 98vw;
            }}
        }}

        section[data-testid="stSidebar"] {{
            background-color: var(--mf-sidebar);
            border-right: 1px solid var(--mf-line);
        }}
        section[data-testid="stSidebar"] .block-container {{
            padding-top: 1.6rem;
        }}

        h1, h2, h3, h4, h5, h6 {{
            font-family: 'Segoe UI', 'Source Sans Pro', 'Microsoft YaHei', sans-serif;
            color: var(--mf-ink) !important;
            font-weight: 700;
            letter-spacing: -0.01em;
        }}
        p, span, label, div {{
            font-family: 'Segoe UI', 'Source Sans Pro', 'Microsoft YaHei', sans-serif;
        }}

        /* CPC-style custom cards */
        .mf-hero {{
            border: 1px solid var(--mf-line-2);
            background: rgba(255,255,255,0.98);
            border-radius: 18px;
            padding: 22px 26px;
            margin: 0 0 18px 0;
            box-shadow: 0 8px 26px rgba(79, 131, 204, 0.045);
        }}
        .mf-fig-label {{
            color: #3D6F9F;
            font-weight: 800;
            font-size: 0.95rem;
            margin-bottom: 2px;
        }}
        .mf-title {{
            color: var(--mf-ink);
            font-size: 2.0rem;
            line-height: 1.15;
            font-weight: 800;
            margin: 0;
        }}
        .mf-subtitle {{
            color: var(--mf-muted);
            font-size: 1.02rem;
            margin-top: 6px;
        }}
        .mf-topbar {{
            border: 1px solid var(--mf-line);
            background: rgba(255,255,255,0.88);
            border-radius: 14px;
            padding: 12px 14px;
            margin: 0 0 12px 0;
            box-shadow: 0 2px 12px rgba(44, 62, 80, 0.025);
            display: flex;
            align-items: stretch;
            justify-content: space-between;
            gap: 18px;
            overflow: visible;
        }}
        .mf-topbar-left {{
            display: flex;
            align-items: center;
            gap: 14px;
            min-width: 260px;
            padding-right: 18px;
            border-right: 1px solid var(--mf-line);
        }}
        .mf-logo-slot {{
            width: 64px;
            height: 64px;
            border-radius: 14px;
            border: 1px solid rgba(93,143,186,0.24);
            background: #101014;
            display: flex;
            align-items: center;
            justify-content: center;
            overflow: hidden;
            flex: 0 0 auto;
            box-shadow: 0 3px 14px rgba(44,62,80,0.10);
        }}
        .mf-logo-slot img {{
            width: 100%;
            height: 100%;
            object-fit: cover;
            display: block;
        }}
        .mf-logo-fallback {{
            color: #ffffff;
            font-weight: 900;
            letter-spacing: -0.04em;
            font-size: 1.18rem;
        }}
        .mf-brand-title {{
            font-size: 1.18rem;
            line-height: 1.10;
            font-weight: 850;
            color: var(--mf-ink);
            margin: 0;
        }}
        .mf-brand-subtitle {{
            color: var(--mf-muted);
            font-size: 0.84rem;
            margin-top: 3px;
        }}
        .mf-topbar-metrics {{
            flex: 1 1 auto;
            display: grid;
            grid-template-columns: repeat(4, minmax(150px, 1fr));
            gap: 10px;
            align-items: stretch;
            min-width: 0;
        }}
        .mf-top-metric {{
            border: 1px solid var(--mf-line);
            background: rgba(248, 251, 255, 0.78);
            border-radius: 12px;
            padding: 10px 13px;
            min-height: 58px;
            display: flex;
            flex-direction: column;
            justify-content: center;
        }}
        .mf-top-metric-label {{
            color: var(--mf-muted);
            font-size: 0.80rem;
            font-weight: 700;
            margin-bottom: 3px;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        .mf-top-metric-value {{
            color: var(--mf-control-dark);
            font-family: Consolas, 'Courier New', monospace;
            font-weight: 850;
            font-size: 1.10rem;
            line-height: 1.12;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        @media (max-width: 1200px) {{
            .mf-topbar {{ align-items: flex-start; flex-direction: column; }}
            .mf-topbar-left {{ border-right: none; padding-right: 0; }}
            .mf-topbar-metrics {{ width: 100%; grid-template-columns: repeat(2, minmax(150px, 1fr)); }}
        }}
        @media (max-width: 700px) {{
            .mf-topbar-metrics {{ grid-template-columns: 1fr; }}
        }}

        /* v2.2 compact header: large logo + 2x2 metrics */
        .mf-header-v22 {{
            border: 1px solid var(--mf-line);
            background: rgba(255,255,255,0.92);
            border-radius: 14px;
            padding: 14px 16px;
            margin: 0 0 12px 0;
            box-shadow: 0 2px 12px rgba(44, 62, 80, 0.025);
            display: grid;
            grid-template-columns: minmax(260px, 330px) minmax(420px, 1fr);
            gap: 16px;
            align-items: stretch;
        }}
        .mf-header-brand-v22 {{
            display: flex;
            align-items: center;
            gap: 14px;
            min-width: 0;
            padding-right: 16px;
            border-right: 1px solid var(--mf-line);
        }}
        .mf-logo-slot-v22 {{
            width: 84px;
            height: 84px;
            border-radius: 13px;
            border: 1px solid rgba(93,143,186,0.22);
            background: #08080a;
            display: flex;
            align-items: center;
            justify-content: center;
            overflow: hidden;
            flex: 0 0 84px;
            box-shadow: 0 3px 14px rgba(44,62,80,0.12);
        }}
        .mf-logo-slot-v22 img {{
            width: 100%;
            height: 100%;
            object-fit: cover;
            display: block;
        }}
        .mf-brand-text-v22 {{ min-width: 0; }}
        .mf-brand-title-v22 {{
            font-size: 1.30rem;
            line-height: 1.12;
            font-weight: 850;
            color: var(--mf-ink);
            margin: 0;
            white-space: nowrap;
        }}
        .mf-brand-subtitle-v22 {{
            color: var(--mf-muted);
            font-size: 0.88rem;
            margin-top: 4px;
            white-space: nowrap;
        }}
        .mf-metric-grid-v22 {{
            display: grid;
            grid-template-columns: repeat(2, minmax(180px, 1fr));
            gap: 10px;
            align-items: stretch;
        }}
        .mf-metric-v22 {{
            border: 1px solid var(--mf-line);
            background: rgba(248, 251, 255, 0.78);
            border-radius: 12px;
            padding: 9px 12px;
            min-height: 38px;
            display: flex;
            flex-direction: column;
            justify-content: center;
        }}
        .mf-metric-label-v22 {{
            color: var(--mf-muted);
            font-size: 0.78rem;
            font-weight: 700;
            margin-bottom: 2px;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        .mf-metric-value-v22 {{
            color: var(--mf-control-dark);
            font-family: Consolas, 'Courier New', monospace;
            font-weight: 850;
            font-size: 1.05rem;
            line-height: 1.10;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        @media (max-width: 1050px) {{
            .mf-header-v22 {{ grid-template-columns: 1fr; }}
            .mf-header-brand-v22 {{ border-right: none; padding-right: 0; }}
        }}
        @media (max-width: 620px) {{
            .mf-metric-grid-v22 {{ grid-template-columns: 1fr; }}
            .mf-logo-slot-v22 {{ width: 72px; height: 72px; flex-basis: 72px; }}
        }}
        .mf-card {{
            border: 1px solid var(--mf-line);
            background: var(--mf-panel);
            border-radius: 14px;
            padding: 14px 16px;
            margin: 8px 0 12px 0;
            box-shadow: 0 2px 10px rgba(44, 62, 80, 0.025);
        }}
        .mf-card-title {{
            color: var(--mf-ink);
            font-weight: 800;
            font-size: 1.05rem;
            margin-bottom: 4px;
        }}
        .mf-card-subtitle {{
            color: var(--mf-muted);
            font-size: 0.92rem;
        }}
        .mf-soft-blue {{ background: var(--mf-soft-blue); border-color: #DDEAFF; }}
        .mf-soft-green {{ background: var(--mf-soft-green); border-color: #CDEED8; }}
        .mf-soft-orange {{ background: var(--mf-soft-orange); border-color: #F5DFB7; }}
        .mf-soft-purple {{ background: var(--mf-soft-purple); border-color: #DED8FF; }}
        .mf-pill {{
            display: inline-block;
            padding: 5px 12px;
            border: 1px solid var(--mf-line-2);
            border-radius: 999px;
            color: var(--mf-muted);
            background: #fff;
            font-size: 0.82rem;
            font-weight: 700;
            margin-right: 6px;
            margin-top: 6px;
        }}
        .mf-help-icon {{
            position: relative;
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 18px;
            height: 18px;
            margin-left: 8px;
            border-radius: 50%;
            border: 1px solid rgba(110, 168, 213, 0.35);
            color: #3D6F9F;
            font-size: 0.72rem;
            font-weight: 800;
            background: rgba(255,255,255,0.88);
            cursor: help;
            vertical-align: middle;
            overflow: visible;
            z-index: 20;
        }}
        .mf-help-icon::after {{
            content: attr(data-tooltip);
            position: absolute;
            top: 125%;
            left: 50%;
            transform: translateX(-50%);
            width: 300px;
            max-width: 45vw;
            padding: 9px 11px;
            border-radius: 10px;
            border: 1px solid var(--mf-line-2);
            background: rgba(255, 255, 255, 0.98);
            color: var(--mf-ink);
            box-shadow: 0 8px 24px rgba(44, 62, 80, 0.12);
            font-size: 0.80rem;
            line-height: 1.45;
            font-weight: 500;
            text-align: left;
            white-space: normal;
            opacity: 0;
            visibility: hidden;
            pointer-events: none;
            transition: opacity 120ms ease-in-out;
            z-index: 9999;
        }}
        .mf-help-icon:hover::after {{
            opacity: 1;
            visibility: visible;
        }}
        .mf-card, .mf-hero {{ overflow: visible; }}
        /* Unified selected-filter tags: visible but less aggressive than Streamlit red */
        div[data-baseweb="tag"],
        span[data-baseweb="tag"] {{
            background-color: var(--mf-control-bg) !important;
            color: var(--mf-control-dark) !important;
            border: 1px solid rgba(93, 143, 186, 0.46) !important;
            box-shadow: none !important;
        }}
        div[data-baseweb="tag"] span,
        span[data-baseweb="tag"] span {{
            color: var(--mf-control-dark) !important;
            font-weight: 650 !important;
        }}
        div[data-baseweb="tag"] svg,
        span[data-baseweb="tag"] svg {{
            color: var(--mf-control-dark) !important;
            fill: var(--mf-control-dark) !important;
        }}
        /* Multiselect dropdown and expander accent cleanup */
        .stMultiSelect [data-baseweb="select"] svg,
        .stSelectbox [data-baseweb="select"] svg {{
            color: var(--mf-control-dark) !important;
            fill: var(--mf-control-dark) !important;
        }}
        div[data-testid="stExpander"] summary,
        div[data-testid="stExpander"] summary span,
        div[data-testid="stExpander"] summary svg {{
            color: var(--mf-ink) !important;
            fill: var(--mf-control-dark) !important;
        }}
        /* Slider accent: medium blue, not red and not too pale */
        div[data-testid="stSlider"] [role="slider"] {{
            background-color: var(--mf-control) !important;
            border: 2px solid #FFFFFF !important;
            box-shadow: 0 0 0 1px rgba(93, 143, 186, 0.42) !important;
        }}
        div[data-testid="stSlider"] div[data-baseweb="slider"] div {{
            color: var(--mf-control-dark) !important;
        }}
        div[data-testid="stSlider"] div[data-baseweb="slider"] div[style*="background"] {{
            background-color: rgba(93, 143, 186, 0.72) !important;
        }}
        .mf-pill-active {{
            background: rgba(110, 168, 213, 0.16);
            color: #3D6F9F;
            border-color: rgba(110, 168, 213, 0.42);
        }}
        .mf-badge {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            height: 28px;
            min-width: 28px;
            padding: 0 9px;
            border-radius: 999px;
            background: var(--mf-blue);
            color: #fff;
            font-weight: 800;
            font-size: 0.78rem;
        }}
        .mf-workflow-step {{
            display: grid;
            grid-template-columns: 34px auto;
            gap: 10px;
            align-items: start;
            margin: 10px 0;
        }}
        .mf-step-dot {{
            width: 26px;
            height: 26px;
            border-radius: 50%;
            background: var(--mf-blue);
            color: #fff;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: 800;
            font-size: 0.78rem;
            box-shadow: 0 0 0 4px #fff;
        }}
        .mf-step-title {{
            color: var(--mf-ink);
            font-weight: 800;
            line-height: 1.1;
        }}
        .mf-step-note {{
            color: var(--mf-muted);
            font-size: 0.86rem;
            margin-top: 2px;
        }}
        .selected-card {{
            border: 1px solid #DDEAFF;
            border-left: 4px solid rgba(79, 131, 204, 0.48);
            background: rgba(241, 246, 255, 0.82);
            border-radius: 12px;
            padding: 12px 14px;
            margin: 0 0 14px 0;
        }}
        .selected-title {{
            font-weight: 800;
            color: #3D6F9F;
            margin-bottom: 4px;
            font-size: 1.02rem;
        }}
        .selected-subtitle {{
            color: var(--mf-ink);
            font-size: 0.92rem;
        }}
        .small-muted {{
            color: var(--mf-muted);
            font-size: 0.84rem;
        }}
        .mf-codebox {{
            border: 1px solid var(--mf-line-2);
            background: #F8FAFC;
            border-radius: 12px;
            padding: 12px 14px;
            font-family: Consolas, 'Courier New', monospace;
            color: var(--mf-ink);
            white-space: pre-wrap;
            font-size: 0.88rem;
        }}

        /* Metrics mimic Fig. 3 cards */
        div[data-testid="stMetric"] {{
            background-color: #ffffff;
            padding: 12px 14px;
            border-radius: 12px;
            border: 1px solid var(--mf-line);
            box-shadow: 0 2px 10px rgba(44,62,80,0.025);
        }}
        div[data-testid="stMetricLabel"] label {{
            color: var(--mf-muted) !important;
            font-size: 0.82rem;
            font-weight: 700;
        }}
        div[data-testid="stMetricValue"] div {{
            color: #3D6F9F !important;
            font-family: 'Consolas', monospace;
            font-weight: 800;
            font-size: 1.28rem;
        }}

        /* Tabs as rounded pills */
        .stTabs [data-baseweb="tab-list"] {{
            gap: 8px;
            border-bottom: none;
            margin-bottom: 6px;
        }}
        .stTabs [data-baseweb="tab"] {{
            background-color: var(--mf-gray-tab);
            color: var(--mf-ink);
            border-radius: 999px;
            border: 1px solid var(--mf-line);
            padding: 8px 18px;
            font-weight: 700;
        }}
        .stTabs [data-baseweb="tab"][aria-selected="true"] {{
            background-color: rgba(79, 131, 204, 0.13);
            border-color: rgba(110, 168, 213, 0.35);
            color: #3D6F9F !important;
            box-shadow: inset 0 -2px 0 rgba(110, 168, 213, 0.45);
        }}

        /* Inputs */
        .stTextInput input,
        .stNumberInput input,
        .stSelectbox div[data-baseweb="select"],
        .stMultiSelect div[data-baseweb="select"],
        textarea {{
            color: var(--mf-ink) !important;
            border-radius: 9px !important;
        }}
        button[kind="primary"] {{
            background-color: rgba(110, 168, 213, 0.18);
            border: 1px solid rgba(110, 168, 213, 0.36);
            color: #3D6F9F !important;
            border-radius: 999px;
            font-weight: 800;
        }}
        button[kind="primary"]:hover {{
            background-color: rgba(110, 168, 213, 0.26);
            border-color: rgba(110, 168, 213, 0.52);
            color: #3D6F9F !important;
        }}
        button[kind="secondary"] {{
            border-radius: 999px;
            border-color: var(--mf-line-2);
        }}
        div[data-testid="stDataFrame"] {{
            border: 1px solid var(--mf-line);
            border-radius: 12px;
            overflow: hidden;
        }}
        .stDownloadButton button {{
            border-radius: 999px;
            color: #3D6F9F !important;
            border-color: rgba(79, 131, 204, 0.28) !important;
            background-color: rgba(255,255,255,0.76) !important;
            font-weight: 700;
        }}
        {publication_css}
        </style>
        """,
        unsafe_allow_html=True,
    )


CPC_COLORS = {
    "blue": "#6EA8D5",
    "blue2": "#9BC3E2",
    "cyan": "#89D1DA",
    "purple": "#B7AFDF",
    "green": "#93CFA7",
    "orange": "#EDCB89",
    "red": "#DA8A8A",
    "ink": "#2C3E50",
    "muted": "#6C757D",
}

CPC_SEQUENCE = [
    "rgba(110, 168, 213, 0.58)",
    "rgba(137, 209, 218, 0.56)",
    "rgba(183, 175, 223, 0.56)",
    "rgba(147, 207, 167, 0.56)",
    "rgba(237, 203, 137, 0.58)",
    "rgba(155, 195, 226, 0.56)",
]

# ================= 通用辅助函数 =================
def decode_h5_value(value):
    """Decode common HDF5 attribute values into Python-native objects."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "S":
            return [x.decode("utf-8") for x in value]
        if value.ndim == 0:
            return value.item()
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def to_jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float) and np.isnan(obj):
        return None
    return obj


def make_traj_id(internal_path: str) -> str:
    path = str(internal_path).replace("\\", "/").strip("/")
    if not path:
        return "Unknown"
    tail = path.split("/")[-1]
    return tail if tail else path.replace("/", "_")


def format_ker(value) -> str:
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.2f}"
    except Exception:
        return "NA"


def get_selection_rows(event) -> list:
    try:
        return list(event.selection.get("rows", []))
    except Exception:
        return []


def get_selected_row(event, df: pd.DataFrame):
    rows = get_selection_rows(event)
    if len(rows) != 1 or df.empty:
        return None
    pos = rows[0]
    if 0 <= pos < len(df):
        return df.iloc[pos]
    return None


def persist_selected_record(row) -> None:
    """Store a stable HDF5 record path so reruns/tabs do not drop the selected trajectory."""
    if row is not None and "file" in row:
        st.session_state["selected_record_path"] = str(row["file"])


def resolve_selected_row(event_row, filtered: pd.DataFrame):
    """Resolve the active selected row from a new UI event, previous session state, or fallback row."""
    if filtered is None or filtered.empty:
        st.session_state.pop("selected_record_path", None)
        return None

    if event_row is not None:
        persist_selected_record(event_row)
        return event_row

    selected_path = st.session_state.get("selected_record_path")
    if selected_path and "file" in filtered.columns:
        matched = filtered[filtered["file"].astype(str) == str(selected_path)]
        if not matched.empty:
            return matched.iloc[0]

    fallback = filtered.iloc[0]
    persist_selected_record(fallback)
    return fallback


def reset_filter_widgets():
    for key in FILTER_KEYS:
        st.session_state.pop(key, None)


def clear_invalid_filter_state(all_states, all_frag_string, delta_energy_min, delta_energy_max, nfrag_max):
    """Remove stale widget states when a new HDF5 file has different options."""
    if "gf_state" in st.session_state:
        current = st.session_state["gf_state"]
        if any(x not in all_states for x in current):
            st.session_state.pop("gf_state", None)

    if "gf_frag" in st.session_state:
        current = st.session_state["gf_frag"]
        if any(x not in all_frag_string for x in current):
            st.session_state.pop("gf_frag", None)

    if "gf_ker" in st.session_state:
        try:
            lo, hi = st.session_state["gf_ker"]
            if lo < delta_energy_min or hi > delta_energy_max or lo > hi:
                st.session_state.pop("gf_ker", None)
        except Exception:
            st.session_state.pop("gf_ker", None)

    if "gf_nfrag" in st.session_state:
        try:
            lo, hi = st.session_state["gf_nfrag"]
            if lo < 0 or hi > int(nfrag_max) or lo > hi:
                st.session_state.pop("gf_nfrag", None)
        except Exception:
            st.session_state.pop("gf_nfrag", None)


def build_filter_summary(state_sel, frag_sel, de_range, min_steps, nfrag_range, all_states, all_frag_string, nfrag_max, has_ker=True):
    items = []
    if len(state_sel) != len(all_states):
        items.append(f"{tr('电子态', 'State')}={len(state_sel)}")
    if len(frag_sel) != len(all_frag_string):
        items.append(f"{tr('产物', 'Products')}={len(frag_sel)}")
    if has_ker:
        items.append(f"KER={de_range[0]:.2f}-{de_range[1]:.2f} eV")
    if min_steps > 1:
        items.append(f"{tr('最小步数', 'Min steps')}={min_steps}")
    if nfrag_range[0] != 0 or nfrag_range[1] != int(nfrag_max):
        items.append(f"{tr('片段数', 'Fragments')}={nfrag_range[0]}-{nfrag_range[1]}")
    return "; ".join(items) if items else tr("默认筛选", "Default filters")



def make_filter_context(state_sel, frag_sel, de_range, min_steps, nfrag_range, source="current_filtered_subset"):
    return {
        "record_source": source,
        "state": list(map(str, state_sel)),
        "products": list(map(str, frag_sel)),
        "energy_descriptor_range": [float(de_range[0]), float(de_range[1])] if de_range is not None else None,
        "min_steps": int(min_steps),
        "fragment_count_range": [int(nfrag_range[0]), int(nfrag_range[1])] if nfrag_range is not None else None,
    }


@st.cache_data(show_spinner=False)
def get_energy_descriptor_metadata(h5_path: str):
    """Read file-level descriptor/unit metadata. Falls back to delta_energy/eV."""
    meta = {
        "name": "delta_energy",
        "display_name": "delta_energy",
        "unit": "eV",
        "definition": "",
        "has_definition": False,
    }
    if not os.path.exists(h5_path):
        return meta

    name_keys = ["energy_descriptor_name", "energy_descriptor", "delta_energy_name", "ker_descriptor_name"]
    unit_keys = ["energy_unit", "delta_energy_unit", "ker_unit", "KER_unit"]
    def_keys = ["energy_descriptor_definition", "delta_energy_definition", "ker_definition", "KER_definition"]

    try:
        with h5py.File(h5_path, "r") as f:
            attrs = f.attrs
            for key in name_keys:
                if key in attrs:
                    meta["name"] = str(decode_h5_value(attrs[key]))
                    meta["display_name"] = meta["name"]
                    break
            for key in unit_keys:
                if key in attrs:
                    meta["unit"] = str(decode_h5_value(attrs[key]))
                    break
            for key in def_keys:
                if key in attrs:
                    definition = str(decode_h5_value(attrs[key])).strip()
                    meta["definition"] = definition
                    meta["has_definition"] = bool(definition)
                    break
    except Exception as exc:
        meta["definition"] = f"Unable to read descriptor metadata: {exc}"
        meta["has_definition"] = False
    return meta


def energy_axis_label(energy_meta: dict) -> str:
    name = energy_meta.get("display_name") or energy_meta.get("name") or "delta_energy"
    unit = energy_meta.get("unit") or ""
    return f"{name} ({unit})" if unit else str(name)


def energy_metric_label(energy_meta: dict) -> str:
    name = energy_meta.get("display_name") or energy_meta.get("name") or "delta_energy"
    if str(name).lower() in {"ker", "kinetic energy release"}:
        return tr("平均 KER", "Average KER")
    return tr(f"平均 {name}", f"Average {name}")


@st.cache_data(show_spinner=False)
def validate_hdf5_schema(h5_path: str):
    """Validate the minimal fragment-event HDF5 contract used by MolFragApp."""
    required_group_attrs = ["state", "nstep", "nfrag", "frag_list", "delta_energy"]
    required_datasets = ["geometry"]
    recommended_file_attrs = [
        "schema_name",
        "schema_version",
        "length_unit",
        "velocity_unit",
        "time_unit",
        "energy_unit",
        "fragment_ordering_rule",
        "energy_descriptor_definition",
    ]
    recommended_group_datasets = ["velocity", "velocities", "masses", "time", "energy", "energies", "expec_out_raw"]

    report = {
        "status": "missing",
        "hdf5_path": h5_path,
        "trajectory_group_count": 0,
        "file_attributes": {},
        "missing_recommended_file_attributes": recommended_file_attrs,
        "missing_required_by_group": [],
        "dataset_coverage": {},
        "notes": [],
    }
    if not os.path.exists(h5_path):
        report["notes"].append("HDF5 file does not exist.")
        return report

    try:
        with h5py.File(h5_path, "r") as f:
            report["file_attributes"] = {str(k): to_jsonable(decode_h5_value(v)) for k, v in f.attrs.items()}
            report["missing_recommended_file_attributes"] = [k for k in recommended_file_attrs if k not in f.attrs]
            coverage = {key: 0 for key in sorted(set(required_datasets + recommended_group_datasets))}

            def visitor(name, node):
                if isinstance(node, h5py.Group) and "geometry" in node:
                    report["trajectory_group_count"] += 1
                    missing_attrs = [k for k in required_group_attrs if k not in node.attrs]
                    missing_dsets = [k for k in required_datasets if k not in node]
                    if missing_attrs or missing_dsets:
                        report["missing_required_by_group"].append(
                            {"group": name, "missing_attrs": missing_attrs, "missing_datasets": missing_dsets}
                        )
                    for key in coverage:
                        if key in node:
                            coverage[key] += 1
                        elif key in node.attrs:
                            coverage[key] += 1

            f.visititems(visitor)
            report["dataset_coverage"] = coverage

        if report["trajectory_group_count"] == 0:
            report["status"] = "missing"
            report["notes"].append("No trajectory group containing a geometry dataset was found.")
        elif report["missing_required_by_group"]:
            report["status"] = "partial"
        elif report["missing_recommended_file_attributes"]:
            report["status"] = "usable"
        else:
            report["status"] = "complete"
    except Exception as exc:
        report["status"] = "error"
        report["notes"].append(str(exc))
    return report


def render_schema_status(report: dict):
    status = report.get("status", "missing")
    labels = {
        "complete": tr("Schema complete", "Schema complete"),
        "usable": tr("Schema usable", "Schema usable"),
        "partial": tr("Schema partial", "Schema partial"),
        "missing": tr("Schema missing", "Schema missing"),
        "error": tr("Schema error", "Schema error"),
    }
    card_class = "mf-soft-green" if status in {"complete", "usable"} else "mf-soft-orange"
    n_groups = report.get("trajectory_group_count", 0)
    missing_req = len(report.get("missing_required_by_group", []))
    missing_rec = len(report.get("missing_recommended_file_attributes", []))
    st.markdown(
        f"""
        <div class="mf-card {card_class}">
            <div class="mf-card-title">{labels.get(status, status)}</div>
            <div class="mf-card-subtitle">{n_groups} trajectory groups · {missing_req} required-issue groups · {missing_rec} recommended file attributes missing</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    with st.expander(tr("Schema details", "Schema details"), expanded=False):
        st.json(to_jsonable(report))


def summarize_energy_descriptor(energy_meta: dict):
    definition = energy_meta.get("definition") or tr(
        "未在 HDF5 file attributes 中找到 energy_descriptor_definition；界面将其作为 stored delta_energy descriptor 显示。",
        "No energy_descriptor_definition was found in HDF5 file attributes; the UI displays it as a stored delta_energy descriptor.",
    )
    return {
        tr("描述符", "Descriptor"): energy_meta.get("name", "delta_energy"),
        tr("单位", "Unit"): energy_meta.get("unit", "eV"),
        tr("定义", "Definition"): definition,
    }


def parse_expec_table(input_data, is_content=False) -> pd.DataFrame:
    """Parse SHARC-like expec.out text into a DataFrame."""
    if is_content:
        f_obj = io.StringIO(input_data)
    else:
        if not os.path.exists(input_data):
            return pd.DataFrame()
        f_obj = open(input_data, "r")

    try:
        _ = f_obj.readline()
        title_row = f_obj.readline()
        names = [col_header.strip().replace(" ", "") for col_header in title_row.lstrip("#! ").split("|")]
        names = [n for n in names if n]
        f_obj.seek(0)
        data = pd.read_csv(f_obj, sep=r"\s+", skiprows=3, names=names, engine="python")
        if data.empty:
            return data
        if "Time" not in data.columns:
            data.insert(0, "Time", np.arange(len(data), dtype=float))
        return data
    finally:
        if not is_content:
            f_obj.close()


def get_traj_energy_column(data: pd.DataFrame):
    if "Epot" in data.columns:
        return "Epot"
    if "TotEnergy" in data.columns:
        return "TotEnergy"
    return None


# ================= 能量曲线绘制 =================
def show_energy_curve_interactive(
    input_data,
    incoord_data=None,
    is_content=False,
    height=350,
    highlight_idx=-1,
    key_suffix="",
):
    """Interactive energy curve with preset modes and a single current-frame highlight."""
    try:
        if isinstance(input_data, pd.DataFrame):
            data = input_data.copy()
        else:
            data = parse_expec_table(input_data, is_content=is_content)
        if data.empty:
            st.warning(tr("能量数据为空。", "Energy data is empty."))
            return data

        all_data_cols = [c for c in data.columns if c not in ["Time", "Iter", "Step"]]
        energy_cols = [c for c in all_data_cols if "Energy" in c]
        traj_col = get_traj_energy_column(data)

        mode_labels = {
            "adiabatic": tr("绝热能级 + 当前态轨迹", "Adiabatic energies + current-state path"),
            "trajectory": tr("仅当前态轨迹", "Current-state path only"),
            "custom": tr("自定义列", "Custom columns"),
        }
        view_mode = st.selectbox(
            tr("能量显示模式", "Energy view mode"),
            options=list(mode_labels.keys()),
            format_func=lambda k: mode_labels[k],
            key=f"energy_mode_{key_suffix}",
        )

        if view_mode == "adiabatic":
            options = energy_cols if energy_cols else all_data_cols[:2]
            st.caption(tr("默认显示绝热能级，并叠加真实动力学路径。", "Adiabatic energy curves are shown with the propagated trajectory path overlaid."))
        elif view_mode == "trajectory":
            options = []
            st.caption(tr("只显示当前态能量路径，适合检查单条轨迹的传播。", "Only the propagated current-state energy path is shown."))
        else:
            default_cols = energy_cols if energy_cols else all_data_cols[:2]
            options = st.multiselect(
                tr("选择展示数据", "Select data to display"),
                options=all_data_cols,
                default=default_cols,
                key=f"energy_sel_{key_suffix}",
            )

        if incoord_data is None:
            fig = go.Figure()
        else:
            from plotly.subplots import make_subplots

            fig = make_subplots(specs=[[{"secondary_y": True}]])
            fig.add_trace(
                go.Scatter(
                    x=incoord_data[0],
                    y=incoord_data[1],
                    name=incoord_data[2],
                    line=dict(dash="dot", color="gray", width=1),
                ),
                secondary_y=True,
            )

        for key in options:
            if key in data.columns:
                fig.add_trace(
                    go.Scatter(
                        x=data["Time"],
                        y=data[key],
                        mode="lines",
                        name=key,
                        line=dict(width=1.5),
                        opacity=0.66,
                    )
                )

        if traj_col:
            fig.add_trace(
                go.Scatter(
                    x=data["Time"],
                    y=data[traj_col],
                    mode="markers",
                    marker=dict(color="black", size=4),
                    name=tr("当前态轨迹", "Current-state trajectory"),
                    legendrank=1,
                )
            )

        if highlight_idx >= 0 and highlight_idx < len(data):
            current_time = data["Time"].iloc[highlight_idx]
            fig.add_vline(x=current_time, line_width=1, line_dash="dash", line_color="red", opacity=0.50)
            if traj_col:
                cur_val = data[traj_col].iloc[highlight_idx]
                fig.add_trace(
                    go.Scatter(
                        x=[current_time],
                        y=[cur_val],
                        mode="markers",
                        marker=dict(color="red", size=14, symbol="circle-open", line=dict(width=3)),
                        name=tr("当前帧", "Current frame"),
                        showlegend=False,
                        hoverinfo="skip",
                    )
                )

        fig.update_layout(
            xaxis_title=tr("时间 (fs)", "Time (fs)"),
            yaxis_title=tr("能量 (eV)", "Energy (eV)"),
            height=height,
            margin=dict(l=0, r=0, t=30, b=0),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            hovermode="x unified",
        )
        st.plotly_chart(fig, use_container_width=True)
        return data

    except Exception as e:
        st.warning(tr(f"无法解析能量数据: {e}", f"Unable to parse energy data: {e}"))
        return pd.DataFrame()


# ================= 3D Mol 播放器 =================
def get_total_frames(mol_block):
    if not mol_block:
        return 0
    try:
        lines = mol_block.splitlines()
        if not lines:
            return 0
        natoms = int(lines[0].strip())
        lines_per_frame = natoms + 2
        return len(lines) // lines_per_frame
    except Exception:
        return 1


def extract_xyz_frame(mol_block: str, frame_idx: int) -> str:
    if not mol_block:
        return ""
    lines = mol_block.splitlines()
    if not lines:
        return ""
    try:
        natoms = int(lines[0].strip())
        lines_per_frame = natoms + 2
        total = len(lines) // lines_per_frame
        if total <= 0:
            return ""
        frame_idx = min(max(frame_idx, 0), total - 1)
        start = frame_idx * lines_per_frame
        return "\n".join(lines[start : start + lines_per_frame]) + "\n"
    except Exception:
        return mol_block


def show_mol_advanced(input_data, mode="auto", frame_idx=0, is_content=False, width=680, height=420):
    if is_content:
        mol_block = input_data
    else:
        if os.path.exists(input_data):
            with open(input_data, "r") as f:
                mol_block = f.read()
        else:
            st.warning(tr(f"未找到文件: {input_data}", f"File not found: {input_data}"))
            return

    try:
        view = py3Dmol.view(width=int(width), height=int(height))
        view.addModelsAsFrames(mol_block, "xyz")
        view.setStyle(
            {
                "stick": {"radius": 0.14, "colorscheme": "Jmol"},
                "sphere": {"radius": 0.36, "colorscheme": "Jmol"},
            }
        )

        if mode == "auto":
            view.animate({"loop": "forward", "step": 1, "retime": True})
        else:
            view.setFrame(frame_idx)

        view.setBackgroundColor("white")
        view.zoomTo()
        try:
            view.zoom(1.55)
        except Exception:
            pass
        html_block = f"""
        <div style="width:100%; display:flex; justify-content:center; align-items:center;">
            {view._make_html()}
        </div>
        """
        components.html(html_block, height=int(height) + 10)

    except Exception as e:
        st.error(tr(f"3D 视图错误: {e}", f"3D view error: {e}"))


# ================= HDF5 读取 =================
@st.cache_data
def load_hdf5_summary(h5_path: str):
    data_list = []
    if not os.path.exists(h5_path):
        return pd.DataFrame()

    def visitor_func(name, node):
        if isinstance(node, h5py.Group) and "geometry" in node:
            attrs = node.attrs

            def get_attr(k, default=None):
                if k not in attrs:
                    return default
                return decode_h5_value(attrs[k])

            f_list_raw = get_attr("frag_list", "")
            if isinstance(f_list_raw, list):
                f_list = [str(x).strip() for x in f_list_raw if str(x).strip()]
            else:
                f_list = [x.strip() for x in str(f_list_raw).split(",")] if f_list_raw else []

            f_json = str(get_attr("frag", "[]"))

            try:
                delta_energy = float(get_attr("delta_energy", np.nan))
            except Exception:
                delta_energy = np.nan

            try:
                nstep = int(get_attr("nstep", 0))
            except Exception:
                nstep = 0

            try:
                nfrag = int(get_attr("nfrag", 0))
            except Exception:
                nfrag = 0

            item = {
                "traj_id": make_traj_id(name),
                "file": name,
                "state": str(get_attr("state", "Unknown")),
                "nstep": nstep,
                "nfrag": nfrag,
                "frag_list": f_list,
                "frag_string": " + ".join(f_list) if f_list else "Bound",
                "frag": f_json,
                "delta_energy": delta_energy,
                "expec_out": "internal",
            }
            data_list.append(item)

    with h5py.File(h5_path, "r") as f:
        f.visititems(visitor_func)

    df = pd.DataFrame(data_list)
    if not df.empty:
        df = df.sort_values(["state", "traj_id"]).reset_index(drop=True)
    return df


@st.cache_data(show_spinner=False)
def get_h5_frame_count(h5_path, internal_path):
    if not os.path.exists(h5_path):
        return 0
    with h5py.File(h5_path, "r") as f:
        if internal_path not in f or "geometry" not in f[internal_path]:
            return 0
        return int(f[internal_path]["geometry"].shape[0])


@st.cache_data(show_spinner=False)
def get_xyz_content_from_h5(h5_path, internal_path, stride: int = 1, frame_idx=None):
    """Return an XYZ string. Use frame_idx for single-frame access or stride for preview playback."""
    stride = max(1, int(stride or 1))
    with h5py.File(h5_path, "r") as f:
        if internal_path not in f:
            return ""
        grp = f[internal_path]
        geo_ds = grp["geometry"]
        n_frames = int(geo_ds.shape[0])

        if "elements" in grp.attrs:
            atoms_raw = decode_h5_value(grp.attrs["elements"])
            atoms = atoms_raw if isinstance(atoms_raw, list) else list(atoms_raw)
            atoms = [a.decode("utf-8") if isinstance(a, bytes) else str(a) for a in atoms]
        else:
            atoms = ["X"] * int(geo_ds.shape[1])

        if frame_idx is not None:
            idx = int(frame_idx)
            if idx < 0:
                idx = n_frames + idx
            idx = min(max(idx, 0), max(n_frames - 1, 0))
            frame_indices = [idx]
        else:
            frame_indices = list(range(0, n_frames, stride))
            if frame_indices and frame_indices[-1] != n_frames - 1:
                frame_indices.append(n_frames - 1)

        out = io.StringIO()
        n_atoms = len(atoms)
        for i in frame_indices:
            out.write(f"{n_atoms}\nFrame {i}\n")
            frame = geo_ds[i]
            for j, atom in enumerate(atoms):
                x, y, z = frame[j]
                out.write(f"{atom:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n")
        return out.getvalue()


COVALENT_RADII_FOR_GIF = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}

ATOM_COLORS_FOR_GIF = {
    "H": "#F2F2F2", "C": "#404040", "N": "#4F6BD8", "O": "#D94A4A", "F": "#7BC77A",
    "P": "#E68632", "S": "#E6C84F", "Cl": "#7BC77A", "Br": "#A65A3A", "I": "#8A63A8",
}


def _clean_element_symbol(atom_label: str) -> str:
    text = str(atom_label).strip()
    if not text:
        return "X"
    if len(text) >= 2 and text[:2].capitalize() in COVALENT_RADII_FOR_GIF:
        return text[:2].capitalize()
    return text[0].upper()


@st.cache_data(show_spinner=False)
def get_geometry_frames_for_gif(h5_path, internal_path, start_frame: int, end_frame: int, stride: int = 1):
    """Return atom labels, HDF5 frame indices and geometry frames for GIF export."""
    stride = max(1, int(stride or 1))
    with h5py.File(h5_path, "r") as f:
        if internal_path not in f or "geometry" not in f[internal_path]:
            return [], [], np.zeros((0, 0, 3), dtype=float)
        grp = f[internal_path]
        geo_ds = grp["geometry"]
        n_frames = int(geo_ds.shape[0])
        if n_frames <= 0:
            return [], [], np.zeros((0, 0, 3), dtype=float)

        start = int(start_frame)
        end = int(end_frame)
        if start < 0:
            start = n_frames + start
        if end < 0:
            end = n_frames + end
        start = min(max(start, 0), n_frames - 1)
        end = min(max(end, 0), n_frames - 1)
        if start > end:
            start, end = end, start
        frame_indices = list(range(start, end + 1, stride))
        if frame_indices and frame_indices[-1] != end:
            frame_indices.append(end)

        if "elements" in grp.attrs:
            atoms_raw = decode_h5_value(grp.attrs["elements"])
            atoms = atoms_raw if isinstance(atoms_raw, list) else list(atoms_raw)
            atoms = [a.decode("utf-8") if isinstance(a, bytes) else str(a) for a in atoms]
        else:
            atoms = ["X"] * int(geo_ds.shape[1])

        frames = np.asarray(geo_ds[frame_indices], dtype=float)
        return atoms, frame_indices, frames


def infer_bonds_for_gif(atoms, coords, scale: float = 1.22, max_distance: float = 2.4):
    """Infer simple covalent-radius bonds for a frame; intended only for visualization."""
    bonds = []
    natom = len(atoms)
    for i in range(natom):
        ei = _clean_element_symbol(atoms[i])
        ri = COVALENT_RADII_FOR_GIF.get(ei, 0.75)
        for j in range(i + 1, natom):
            ej = _clean_element_symbol(atoms[j])
            rj = COVALENT_RADII_FOR_GIF.get(ej, 0.75)
            threshold = min(max_distance, scale * (ri + rj))
            dist = float(np.linalg.norm(coords[i] - coords[j]))
            if 0.05 < dist <= threshold:
                bonds.append((i, j))
    return bonds


def build_trajectory_gif(
    h5_path,
    internal_path,
    start_frame: int,
    end_frame: int,
    stride: int = 1,
    fps: int = 12,
    image_size: int = 640,
    show_frame_label: bool = True,
):
    """Render a lightweight 3D molecule GIF from HDF5 geometry frames using Matplotlib/Pillow."""
    from PIL import Image

    atoms, frame_indices, frames = get_geometry_frames_for_gif(
        h5_path, internal_path, start_frame=start_frame, end_frame=end_frame, stride=stride
    )
    if len(frame_indices) == 0 or frames.size == 0:
        return b"", 0

    all_xyz = frames.reshape(-1, 3)
    mins = np.nanmin(all_xyz, axis=0)
    maxs = np.nanmax(all_xyz, axis=0)
    center = (mins + maxs) / 2.0
    span = float(np.nanmax(maxs - mins))
    if not np.isfinite(span) or span < 1.0:
        span = 1.0
    span *= 1.18

    image_size = int(image_size)
    fps = max(1, int(fps or 12))
    duration_ms = int(round(1000.0 / fps))
    marker_size = 170 if len(atoms) <= 5 else 85
    colors = [ATOM_COLORS_FOR_GIF.get(_clean_element_symbol(a), "#8C8C8C") for a in atoms]

    pil_frames = []
    for local_idx, frame_no in enumerate(frame_indices):
        coords = frames[local_idx]
        fig = plt.figure(figsize=(image_size / 100.0, image_size / 100.0), dpi=100)
        ax = fig.add_subplot(111, projection="3d")
        ax.set_facecolor("white")
        fig.patch.set_facecolor("white")

        bonds = infer_bonds_for_gif(atoms, coords)
        for i, j in bonds:
            ax.plot(
                [coords[i, 0], coords[j, 0]],
                [coords[i, 1], coords[j, 1]],
                [coords[i, 2], coords[j, 2]],
                color="#555555",
                linewidth=2.0,
                alpha=0.82,
                zorder=1,
            )

        ax.scatter(
            coords[:, 0], coords[:, 1], coords[:, 2],
            s=marker_size,
            c=colors,
            edgecolors="#202020",
            linewidths=0.6,
            depthshade=True,
            zorder=3,
        )

        if len(atoms) <= 8:
            for atom, xyz in zip(atoms, coords):
                ax.text(xyz[0], xyz[1], xyz[2], f" {atom}", fontsize=9, color="#303030")

        ax.set_xlim(center[0] - span / 2, center[0] + span / 2)
        ax.set_ylim(center[1] - span / 2, center[1] + span / 2)
        ax.set_zlim(center[2] - span / 2, center[2] + span / 2)
        ax.view_init(elev=22, azim=42)
        try:
            ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass
        ax.set_axis_off()
        if show_frame_label:
            ax.set_title(f"Frame {frame_no}", fontsize=12, color="#315F86", pad=8)
        plt.tight_layout(pad=0.04)

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=100, bbox_inches="tight", pad_inches=0.04)
        plt.close(fig)
        buf.seek(0)
        pil_frames.append(Image.open(buf).convert("RGB"))

    out = io.BytesIO()
    pil_frames[0].save(
        out,
        format="GIF",
        save_all=True,
        append_images=pil_frames[1:],
        duration=duration_ms,
        loop=0,
        optimize=True,
    )
    out.seek(0)
    return out.getvalue(), len(pil_frames)


@st.cache_data(show_spinner=False)
def get_energy_dataframe_from_h5(h5_path, internal_path):
    """Read structured HDF5 energy arrays first; fall back to expec_out_raw text."""
    if not os.path.exists(h5_path):
        return pd.DataFrame()
    with h5py.File(h5_path, "r") as f:
        if internal_path not in f:
            return pd.DataFrame()
        grp = f[internal_path]

        time_values = None
        for tname in ["time", "Time", "times"]:
            if tname in grp:
                time_values = np.asarray(grp[tname][()]).reshape(-1)
                break

        for key in ["expec_data", "expectation", "expectations"]:
            if key in grp:
                arr = grp[key][()]
                if getattr(arr.dtype, "names", None):
                    df = pd.DataFrame({name: arr[name] for name in arr.dtype.names})
                else:
                    arr = np.asarray(arr)
                    if arr.ndim == 1:
                        df = pd.DataFrame({key: arr})
                    else:
                        df = pd.DataFrame(arr, columns=[f"{key}_{i+1}" for i in range(arr.shape[1])])
                if "Time" not in df.columns:
                    df.insert(0, "Time", time_values[: len(df)] if time_values is not None and len(time_values) >= len(df) else np.arange(len(df), dtype=float))
                return df

        for key in ["energy", "energies", "energy_data", "adiabatic_energies"]:
            if key in grp:
                arr = np.asarray(grp[key][()])
                if arr.ndim == 1:
                    df = pd.DataFrame({"Energy": arr})
                elif arr.ndim == 2:
                    df = pd.DataFrame(arr, columns=[f"Energy{i+1}" for i in range(arr.shape[1])])
                else:
                    continue
                df.insert(0, "Time", time_values[: len(df)] if time_values is not None and len(time_values) >= len(df) else np.arange(len(df), dtype=float))
                return df

        if "expec_out_raw" in grp:
            raw = grp["expec_out_raw"][()]
            raw_text = raw.decode("utf-8") if isinstance(raw, bytes) else str(raw)
            return parse_expec_table(raw_text, is_content=True)
    return pd.DataFrame()


@st.cache_data(show_spinner=False)
def get_expec_text_from_h5(h5_path, internal_path):
    with h5py.File(h5_path, "r") as f:
        if internal_path not in f:
            return ""
        grp = f[internal_path]
        if "expec_out_raw" in grp:
            raw = grp["expec_out_raw"][()]
            if isinstance(raw, bytes):
                return raw.decode("utf-8")
            return str(raw)
    return ""


@st.cache_data(show_spinner=False)
def get_h5_group_overview(h5_path, internal_path):
    overview = {"datasets": [], "attributes": []}
    if not os.path.exists(h5_path):
        return overview
    with h5py.File(h5_path, "r") as f:
        if internal_path not in f:
            return overview
        grp = f[internal_path]
        overview["attributes"] = sorted(list(grp.attrs.keys()))
        for key, obj in grp.items():
            if isinstance(obj, h5py.Dataset):
                overview["datasets"].append(
                    {
                        "name": key,
                        "type": "dataset",
                        "shape": str(tuple(obj.shape)),
                        "dtype": str(obj.dtype),
                    }
                )
            elif isinstance(obj, h5py.Group):
                overview["datasets"].append(
                    {"name": key, "type": "group", "shape": "", "dtype": ""}
                )
    return overview



# ================= 三体 observable table helpers =================
def parse_fragment_indices(raw_frag):
    """Best-effort parser for fragment index sets stored as JSON/Python-like strings."""
    if raw_frag is None or (isinstance(raw_frag, float) and np.isnan(raw_frag)):
        return []
    value = raw_frag
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return []
        for parser in (json.loads, ast.literal_eval):
            try:
                value = parser(s)
                break
            except Exception:
                value = None
        if value is None:
            return []
    if isinstance(value, dict):
        items = list(value.values())
    else:
        items = list(value) if isinstance(value, (list, tuple, np.ndarray)) else []
    out = []
    for item in items:
        if isinstance(item, dict) and "atoms" in item:
            item = item["atoms"]
        if isinstance(item, (list, tuple, np.ndarray)):
            try:
                out.append([int(x) for x in item])
            except Exception:
                out.append([])
        else:
            try:
                out.append([int(item)])
            except Exception:
                out.append([])
    return out


def get_dataset_or_attr_array(grp, names):
    for name in names:
        if name in grp:
            return np.asarray(grp[name][()])
        if name in grp.attrs:
            return np.asarray(decode_h5_value(grp.attrs[name]))
    return None


@st.cache_data(show_spinner=False)
def build_three_body_observable_table(h5_path: str, records_json: str, orders_json: str, frame_idx: int):
    """Build a table-first export for three-body maps.

    The function attempts to compute Dalitz/Newton coordinates when velocity, masses,
    and fragment index sets are available. If not, it still returns a metadata table
    preserving record identity and plotting parameters.
    """
    records = json.loads(records_json)
    orders = json.loads(orders_json)
    rows = []
    if not os.path.exists(h5_path):
        return pd.DataFrame(rows)

    with h5py.File(h5_path, "r") as f:
        for rec in records:
            path = rec.get("file", "")
            base = {
                "record_key": path,
                "traj_id": rec.get("traj_id", make_traj_id(path)),
                "state": rec.get("state", "Unknown"),
                "product_channel": rec.get("frag_string", ""),
                "stored_delta_energy": rec.get("delta_energy", np.nan),
                "analysis_frame": int(frame_idx),
                "fragment_order_zero_based": str(orders),
                "observable_status": "metadata-only",
            }
            if path not in f:
                base["observable_status"] = "missing HDF5 group"
                rows.append(base)
                continue
            grp = f[path]
            try:
                frag_raw = rec.get("frag", None)
                if not frag_raw or str(frag_raw) in {"[]", "None", "nan"}:
                    frag_raw = decode_h5_value(grp.attrs.get("frag", "[]")) if "frag" in grp.attrs else "[]"
                frag_indices = parse_fragment_indices(frag_raw)
                if len(frag_indices) < 3:
                    base["observable_status"] = "fragment index sets unavailable"
                    rows.append(base)
                    continue

                velocities = get_dataset_or_attr_array(grp, ["velocity", "velocities", "veloc"])
                masses = get_dataset_or_attr_array(grp, ["masses", "mass"])
                if velocities is None or masses is None:
                    base["observable_status"] = "velocity or masses unavailable"
                    rows.append(base)
                    continue

                velocities = np.asarray(velocities, dtype=float)
                masses = np.asarray(masses, dtype=float).reshape(-1)
                n_frames = velocities.shape[0] if velocities.ndim == 3 else 1
                idx = int(frame_idx)
                if idx < 0:
                    idx = n_frames + idx
                idx = min(max(idx, 0), max(n_frames - 1, 0))
                vframe = velocities[idx] if velocities.ndim == 3 else velocities

                ordered_frags = [frag_indices[i] for i in orders]
                frag_masses, momenta, kinetic = [], [], []
                for atom_ids in ordered_frags:
                    atom_ids = [i for i in atom_ids if 0 <= i < len(masses)]
                    if not atom_ids:
                        raise ValueError("empty fragment after applying atom indices")
                    m = masses[atom_ids]
                    p = (m[:, None] * vframe[atom_ids]).sum(axis=0)
                    M = float(m.sum())
                    Ek = float(np.dot(p, p) / (2.0 * M)) if M > 0 else np.nan
                    frag_masses.append(M)
                    momenta.append(p)
                    kinetic.append(Ek)

                P = np.asarray(momenta, dtype=float)
                E = np.asarray(kinetic, dtype=float)
                EK = float(np.nansum(E))
                if EK > 0:
                    dalitz_x = float((E[0] - E[1]) / (np.sqrt(3.0) * EK))
                    dalitz_y = float(E[2] / EK - 1.0 / 3.0)
                else:
                    dalitz_x = np.nan
                    dalitz_y = np.nan

                p0_norm = float(np.linalg.norm(P[0]))
                if p0_norm > 0:
                    ux = P[0] / p0_norm
                    # Choose a stable in-plane y direction from P2; fallback to a Cartesian perpendicular.
                    p2_perp = P[1] - np.dot(P[1], ux) * ux
                    p2_perp_norm = float(np.linalg.norm(p2_perp))
                    if p2_perp_norm == 0:
                        trial = np.array([1.0, 0.0, 0.0]) if abs(ux[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
                        p2_perp = trial - np.dot(trial, ux) * ux
                        p2_perp_norm = float(np.linalg.norm(p2_perp))
                    uy = p2_perp / p2_perp_norm
                    coords = []
                    for p in P:
                        coords.append([float(np.dot(p, ux) / p0_norm), float(np.dot(p, uy) / p0_norm)])
                else:
                    coords = [[np.nan, np.nan], [np.nan, np.nan], [np.nan, np.nan]]

                base.update(
                    {
                        "observable_status": "computed",
                        "computed_frame": idx,
                        "E1": kinetic[0],
                        "E2": kinetic[1],
                        "E3": kinetic[2],
                        "computed_KER_like_sum": EK,
                        "dalitz_x": dalitz_x,
                        "dalitz_y": dalitz_y,
                        "newton_x1": coords[0][0],
                        "newton_y1": coords[0][1],
                        "newton_x2": coords[1][0],
                        "newton_y2": coords[1][1],
                        "newton_x3": coords[2][0],
                        "newton_y3": coords[2][1],
                    }
                )
                rows.append(base)
            except Exception as exc:
                base["observable_status"] = f"failed: {exc}"
                rows.append(base)
    return pd.DataFrame(rows)


# ================= UI 渲染函数 =================
def render_selected_summary(row):
    if row is None:
        st.info(tr("请在左侧轨迹列表中选择一条轨迹。统计页仍会使用当前筛选集。", "Select one trajectory from the left table. The statistics tab still uses the current filtered set."))
        return

    st.markdown(
        f"""
        <div class="selected-card">
            <div class="selected-title">{tr('当前选中轨迹', 'Selected trajectory')}: {row.get('traj_id', 'Unknown')}</div>
            <div class="selected-subtitle">
                {tr('电子态', 'State')}: <b>{row.get('state', 'Unknown')}</b> &nbsp; | &nbsp;
                {tr('产物', 'Products')}: <b>{row.get('frag_string', 'Bound')}</b> &nbsp; | &nbsp;
                KER: <b>{format_ker(row.get('delta_energy'))} eV</b> &nbsp; | &nbsp;
                {tr('步数', 'Steps')}: <b>{row.get('nstep', 0)}</b>
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_empty_data_message(h5_path: str):
    st.info(tr("数据为空或加载失败。", "Dataset is empty or failed to load."))
    st.markdown(
        tr(
            f"""
请检查：

1. HDF5 路径是否正确：`{h5_path}`；
2. 文件中是否存在包含 `geometry` 数据集的轨迹 group；
3. 每条轨迹 group 是否具有 `state`, `nstep`, `nfrag`, `frag_list`, `delta_energy` 等属性。
            """,
            f"""
Please check:

1. Whether the HDF5 path is correct: `{h5_path}`;
2. Whether the file contains trajectory groups with a `geometry` dataset;
3. Whether each trajectory group has attributes such as `state`, `nstep`, `nfrag`, `frag_list`, and `delta_energy`.
            """,
        )
    )


def render_filtered_download(filtered: pd.DataFrame):
    if filtered.empty:
        return
    st.download_button(
        tr("下载当前筛选轨迹 CSV", "Download filtered trajectories as CSV"),
        data=filtered.to_csv(index=False).encode("utf-8"),
        file_name="molfragapp_filtered_trajectories.csv",
        mime="text/csv",
        use_container_width=True,
        key="download_filtered_csv",
    )





def resolve_logo_path(preferred: str = "MolFragApp_logo.jpg") -> str | None:
    """Find a local logo image for the compact header.

    The preferred layout is to keep MolFragApp_logo.jpg in the same directory as
    this app. A small text fallback is used if the file is absent.
    """
    candidates = []
    if preferred:
        candidates.append(preferred)
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        candidates.extend([
            os.path.join(script_dir, preferred),
            os.path.join(script_dir, "assets", preferred),
        ])
    except Exception:
        pass
    candidates.extend([
        "assets/MolFragApp_logo.jpg",
        "/mnt/data/MolFragApp_logo.jpg",
    ])
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return None


def image_to_data_uri(path: str | None) -> str:
    if not path or not os.path.exists(path):
        return ""
    ext = os.path.splitext(path)[1].lower()
    mime = "image/png" if ext == ".png" else "image/jpeg"
    with open(path, "rb") as f:
        payload = base64.b64encode(f.read()).decode("ascii")
    return f"data:{mime};base64,{payload}"


def _short_path(path: str, max_chars: int = 34) -> str:
    name = os.path.basename(str(path)) or str(path)
    if len(name) <= max_chars:
        return name
    return name[: max_chars - 1] + "…"

def render_cpc_page_header(
    h5_path: str,
    subtitle_extra: str = "",
    total_records: int | None = None,
    n_states: int | None = None,
    avg_descriptor: float | None = None,
    max_fragments: int | None = None,
    schema_report: dict | None = None,
    energy_meta: dict | None = None,
    logo_path: str = "MolFragApp_logo.jpg",
):
    """Compact top header with a large logo and a 2x2 metric grid.

    The header keeps only software identity and high-level dataset summary.
    Schema details and energy-descriptor definitions stay in the sidebar and
    Metadata tab, avoiding repeated developer-oriented information on top.
    """
    logo_uri = image_to_data_uri(resolve_logo_path(logo_path))
    if logo_uri:
        logo_html = f'<img src="{logo_uri}" alt="MolFragApp logo" />'
    else:
        logo_html = '<span class="mf-logo-fallback">MF</span>'

    unit = (energy_meta or {}).get("unit", "eV")
    if avg_descriptor is None or not np.isfinite(avg_descriptor):
        avg_text = "NA"
    else:
        avg_text = f"{avg_descriptor:.2f} {unit}"

    metric_items = [
        (tr("总轨迹数", "Total trajectories"), "0" if total_records is None else f"{int(total_records)}"),
        (tr("电子态数量", "Number of states"), "NA" if n_states is None else f"{int(n_states)}"),
        (tr("平均能量描述符", "Mean energy descriptor"), avg_text),
        (tr("最大片段数", "Maximum fragment count"), "NA" if max_fragments is None else f"{int(max_fragments)}"),
    ]
    metrics_html = "".join(
        f'<div class="mf-metric-v22"><div class="mf-metric-label-v22">{html.escape(label)}</div><div class="mf-metric-value-v22">{html.escape(value)}</div></div>'
        for label, value in metric_items
    )

    header_html = f"""
<div class="mf-header-v22">
  <div class="mf-header-brand-v22">
    <div class="mf-logo-slot-v22">{logo_html}</div>
    <div class="mf-brand-text-v22">
      <div class="mf-brand-title-v22">MolFragApp</div>
      <div class="mf-brand-subtitle-v22">{tr('HDF5 碎片事件分析', 'HDF5 fragment-event analysis')}</div>
    </div>
  </div>
  <div class="mf-metric-grid-v22">{metrics_html}</div>
</div>
"""
    st.markdown(header_html, unsafe_allow_html=True)

def render_sidebar_caption(version: str, connected: bool):
    status_class = "mf-soft-green" if connected else "mf-soft-orange"
    status_text = tr("数据集已连接", "Dataset connected") if connected else tr("等待数据集", "Waiting for dataset")
    status_color = "#2FA866" if connected else "#B07A20"
    st.markdown(
        f"""
        <div class="mf-card {status_class}">
            <div class="mf-card-title">{tr('数据连接', 'Data connection')}</div>
            <div class="mf-card-subtitle">{version}</div>
            <div style="margin-top:8px;color:{status_color};font-weight:800;">{status_text}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )



def render_analysis_note(title: str, subtitle: str, accent: str = "blue"):
    cls = {"blue": "mf-soft-blue", "green": "mf-soft-green", "orange": "mf-soft-orange", "purple": "mf-soft-purple"}.get(accent, "mf-soft-blue")
    st.markdown(
        f"""
        <div class="mf-card {cls}">
            <div class="mf-card-title">{title} {help_icon(subtitle)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def render_export_reuse_panel(
    filtered: pd.DataFrame,
    selected_row,
    h5_path: str = "",
    filter_context: dict | None = None,
    schema_report: dict | None = None,
    energy_meta: dict | None = None,
    plot_context: dict | None = None,
):
    export_help = tr(
        "导出当前筛选记录、单轨迹元数据，或下载 notebook/CLI 接口 manifest；后续脚本可读取该 manifest 生成论文图。",
        "Download selected records, trajectory metadata, or a notebook/CLI interface manifest; later scripts can read the manifest to generate publication figures.",
    )
    st.markdown(
        f"""
        <div class="mf-card">
            <div class="mf-card-title">{tr('导出与复用', 'Export and reuse')} {help_icon(export_help)}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    d1, d2, d3 = st.columns(3)
    with d1:
        if not filtered.empty:
            st.download_button(
                tr("CSV table", "CSV table"),
                data=filtered.to_csv(index=False).encode("utf-8"),
                file_name="molfragapp_selected_records.csv",
                mime="text/csv",
                use_container_width=True,
                key="stats_download_selected_records",
            )
    with d2:
        if selected_row is not None:
            raw_dict = to_jsonable(selected_row.to_dict())
            st.download_button(
                tr("Metadata JSON", "Metadata JSON"),
                data=json.dumps(raw_dict, ensure_ascii=False, indent=2).encode("utf-8"),
                file_name=f"{selected_row.get('traj_id', 'trajectory')}_metadata.json",
                mime="application/json",
                use_container_width=True,
                key="stats_download_selected_metadata",
            )
        else:
            st.button(tr("Metadata JSON", "Metadata JSON"), disabled=True, use_container_width=True)
    with d3:
        manifest = {
            "source_hdf5": h5_path,
            "schema_status": schema_report.get("status") if schema_report else None,
            "schema_version": (schema_report or {}).get("file_attributes", {}).get("schema_version"),
            "energy_descriptor": energy_meta or {},
            "active_filters": filter_context or {},
            "plot_settings": plot_context or {},
            "selected_record_count": int(len(filtered)),
            "selected_record_ids": filtered["traj_id"].astype(str).tolist() if "traj_id" in filtered else [],
            "selected_hdf5_paths": filtered["file"].astype(str).tolist() if "file" in filtered else [],
            "active_selected_trajectory": to_jsonable(selected_row.to_dict()) if selected_row is not None else None,
            "figure_export_interface": {
                "notebook_entry": "Read this manifest, reopen source_hdf5, and rebuild statistics, Dalitz plots, Newton diagrams, and publication figures from selected_hdf5_paths.",
                "cli_entry": "python molfrag_cli.py --h5 <source_hdf5> --manifest molfragapp_analysis_manifest.json --export figures tables",
                "notes": "The manifest preserves record identity, active filters, descriptor metadata, and plotting context for notebook/CLI reproduction.",
            },
        }
        st.download_button(
            tr("Notebook/CLI 接口", "Notebook/CLI interface"),
            data=json.dumps(manifest, ensure_ascii=False, indent=2).encode("utf-8"),
            file_name="molfragapp_analysis_manifest.json",
            mime="application/json",
            use_container_width=True,
            key="stats_download_analysis_manifest",
            help=tr("下载当前筛选上下文，供 notebook 或 CLI 后续生成论文图。", "Download the current selection context for later notebook/CLI figure generation."),
        )


# ================= 主界面 =================
def main():
    # --- 1. 侧边栏 ---
    with st.sidebar:
        st.selectbox("Language / 语言", options=LANG_OPTIONS, index=0, key="ui_lang")
        publication_mode = st.checkbox(
            tr("论文截图模式", "Publication / screenshot mode"),
            value=False,
            key="publication_mode",
            help=tr("隐藏部分浏览器装饰并使用更适合论文截图的紧凑布局。", "Use a cleaner layout for manuscript screenshots."),
        )
        st.markdown("---")

        default_h5 = "Data.2.hdf5"
        h5_path = st.text_input(tr("HDF5 文件路径", "HDF5 file path"), value=default_h5, key="sidebar_h5path")
        connected = os.path.exists(h5_path)
        render_sidebar_caption(VERSION, connected)
        schema_report_sidebar = validate_hdf5_schema(h5_path) if connected else {"status": "missing", "hdf5_path": h5_path}
        render_schema_status(schema_report_sidebar)

        template_path = st.text_input(tr("Template 路径", "Template path"), value="*.template", key="sidebar_template")
        try:
            files = glob(template_path)
            cal_method = read_template(files[0]) if files else tr("未知", "Unknown")
        except Exception:
            cal_method = tr("读取错误", "Read error")
        st.markdown(
            f"""
            <div class="mf-card">
                <div class="mf-card-title">{tr('计算级别', 'Level of theory')}</div>
                <div class="mf-card-subtitle">{cal_method}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )

    inject_custom_css(publication_mode=publication_mode)

    # --- 数据加载 ---
    if not os.path.exists(h5_path):
        st.warning(tr(f"等待连接数据源: {h5_path}", f"Waiting for data source: {h5_path}"))
        trajs_data = pd.DataFrame(
            columns=["traj_id", "file", "state", "nstep", "nfrag", "frag_list", "frag_string", "delta_energy"]
        )
    else:
        with st.spinner(tr("正在加载 HDF5 数据...", "Loading HDF5 data...")):
            trajs_data = load_hdf5_summary(h5_path)

    if trajs_data.empty:
        empty_schema_report = validate_hdf5_schema(h5_path) if os.path.exists(h5_path) else {"status": "missing"}
        empty_energy_meta = get_energy_descriptor_metadata(h5_path) if os.path.exists(h5_path) else {"name": "delta_energy", "unit": "eV"}
        render_cpc_page_header(
            h5_path,
            total_records=0,
            n_states=0,
            avg_descriptor=None,
            max_fragments=0,
            schema_report=empty_schema_report,
            energy_meta=empty_energy_meta,
        )
        render_empty_data_message(h5_path)
        return

    schema_report = validate_hdf5_schema(h5_path)
    energy_meta = get_energy_descriptor_metadata(h5_path)
    descriptor_axis = energy_axis_label(energy_meta)
    descriptor_metric = energy_metric_label(energy_meta)

    render_cpc_page_header(
        h5_path,
        total_records=len(trajs_data),
        n_states=trajs_data["state"].nunique(),
        avg_descriptor=float(trajs_data["delta_energy"].mean()) if trajs_data["delta_energy"].notna().any() else None,
        max_fragments=int(trajs_data.nfrag.max()),
        schema_report=schema_report,
        energy_meta=energy_meta,
    )

    nfrag_max = int(trajs_data.nfrag.max())
    all_states = sorted(set(trajs_data.state.astype(str)))
    all_frag_string = sorted(set(trajs_data.frag_string.astype(str)))
    has_ker = trajs_data["delta_energy"].notna().any()
    if has_ker:
        delta_energy_min = float(trajs_data.delta_energy.min())
        delta_energy_max = float(trajs_data.delta_energy.max())
    else:
        delta_energy_min = 0.0
        delta_energy_max = 0.0

    clear_invalid_filter_state(all_states, all_frag_string, delta_energy_min, delta_energy_max, nfrag_max)

    st.markdown("---")

    # --- 核心工作区 ---
    col_nav, col_detail = st.columns([2.8, 5.2], gap="medium")

    # === 左栏：数据导航 ===
    with col_nav:
        st.subheader(tr("轨迹列表", "Trajectory list"))
        st.markdown(
            f"""
            <div class="mf-card mf-soft-blue">
                <div class="mf-card-title">{tr('Trajectory navigation', 'Trajectory navigation')} {help_icon(tr('先在轻量级 metadata 表上定义事件集合，再按需读取 geometry、velocity 和 energy datasets。', 'Define an event subset from lightweight metadata, then load geometry, velocity, and energy datasets on demand.'))}</div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        with st.expander(tr("全局筛选", "Global filters"), expanded=not publication_mode):
            if st.button(tr("重置筛选", "Reset filters"), use_container_width=True, key="reset_filters_btn"):
                reset_filter_widgets()
                st.rerun()

            state_sel = st.multiselect(
                tr("电子态", "Electronic state"),
                options=all_states,
                default=all_states,
                key="gf_state",
            )
            frag_sel = st.multiselect(
                tr("产物类型", "Product type"),
                options=all_frag_string,
                default=all_frag_string,
                key="gf_frag",
            )

            if has_ker and delta_energy_min < delta_energy_max:
                de_range = st.slider(
                    tr("能量描述符范围", "Energy-descriptor range"),
                    min_value=delta_energy_min,
                    max_value=delta_energy_max,
                    value=(delta_energy_min, delta_energy_max),
                    key="gf_ker",
                )
            else:
                de_range = (delta_energy_min, delta_energy_max)
                st.caption(tr("KER 范围固定或不可用。", "KER range is fixed or unavailable."))

            min_steps = st.number_input(tr("最小步数", "Minimum steps"), value=1, min_value=1, key="gf_minstep")

            if nfrag_max > 0:
                nfrag_range = st.slider(
                    tr("片段数量", "Fragment count"),
                    min_value=0,
                    max_value=nfrag_max,
                    value=(0, nfrag_max),
                    key="gf_nfrag",
                )
            else:
                nfrag_range = (0, 0)
                st.caption(tr("片段数范围固定。", "Fragment-count range is fixed."))

        ker_mask = pd.Series(True, index=trajs_data.index)
        if has_ker:
            ker_mask = (trajs_data.delta_energy >= de_range[0]) & (trajs_data.delta_energy <= de_range[1])

        filtered = trajs_data[
            (trajs_data.nstep >= min_steps)
            & (trajs_data.nfrag >= nfrag_range[0])
            & (trajs_data.nfrag <= nfrag_range[1])
            & (trajs_data.state.isin(state_sel))
            & ker_mask
            & (trajs_data.frag_string.isin(frag_sel))
        ].reset_index(drop=True)

        filter_summary = build_filter_summary(
            state_sel, frag_sel, de_range, min_steps, nfrag_range, all_states, all_frag_string, nfrag_max, has_ker=has_ker
        )
        filter_context = make_filter_context(state_sel, frag_sel, de_range, min_steps, nfrag_range)
        st.caption(tr(f"当前筛选：{filter_summary}", f"Active filters: {filter_summary}"))
        st.caption(tr(f"显示 {len(filtered)} / {len(trajs_data)} 条", f"Showing {len(filtered)} / {len(trajs_data)} trajectories"))

        render_filtered_download(filtered)

        event = st.dataframe(
            filtered,
            key="data_table_main",
            use_container_width=True,
            height=520 if not publication_mode else 420,
            hide_index=True,
            column_order=("traj_id", "state", "nstep", "nfrag", "frag_string", "delta_energy"),
            column_config={
                "traj_id": st.column_config.TextColumn(label=tr("轨迹", "Trajectory"), width="medium"),
                "state": st.column_config.TextColumn(label=tr("电子态", "State"), width="small"),
                "nstep": st.column_config.NumberColumn(label=tr("步数", "Steps"), width="small"),
                "nfrag": st.column_config.NumberColumn(label=tr("片段", "Fragments"), width="small"),
                "frag_string": st.column_config.TextColumn(label=tr("产物", "Products"), width="medium"),
                "delta_energy": st.column_config.NumberColumn(label=descriptor_axis, format="%.2f"),
            },
            on_select="rerun",
            selection_mode=["single-row"],
        )

    event_selected_row = get_selected_row(event, filtered)
    selected_row = resolve_selected_row(event_selected_row, filtered)

    # === 右栏：分析工作台 ===
    with col_detail:
        st.subheader(tr("分析工作台", "Analysis workspace"))
        render_selected_summary(selected_row)

        tab1, tab2, tab3, tab4 = st.tabs(
            [
                tr("事件统计", "Event statistics"),
                tr("轨迹检查", "Trajectory inspection"),
                tr("动量空间图", "Momentum-space maps"),
                tr("元数据", "Metadata"),
            ]
        )

        # === Tab 1: 统计 ===
        with tab1:
            if filtered.empty:
                st.warning(tr("当前筛选集为空，请调整左侧筛选条件。", "The current filtered set is empty. Adjust filters on the left."))
            else:
                s1, s2, s3, s4 = st.columns(4)
                s1.metric(tr("筛选轨迹", "Filtered trajectories"), f"{len(filtered)}")
                s2.metric(tr("产物通道", "Product channels"), f"{filtered['frag_string'].nunique()}")
                s3.metric(tr("三体事件", "Three-body events"), f"{int((filtered['nfrag'] == 3).sum())}")
                s4.metric(descriptor_metric, f"{filtered['delta_energy'].mean():.2f} {energy_meta.get('unit', 'eV')}")

                c_chart1, c_chart2 = st.columns(2)
                with c_chart1:
                    st.markdown(f"**{tr('产物分布', 'Product distribution')}**")
                    chart_data = filtered.groupby(["frag_string", "state"]).size().reset_index(name="count")
                    try:
                        colors_seq = generate_colors(len(sorted(chart_data["state"].unique())))
                    except Exception:
                        colors_seq = px.colors.qualitative.Plotly

                    fig = px.bar(
                        chart_data,
                        x="frag_string",
                        y="count",
                        color="state",
                        color_discrete_sequence=CPC_SEQUENCE,
                        labels={"frag_string": tr("产物", "Products"), "count": tr("数量", "Count"), "state": tr("电子态", "State")},
                    )
                    fig.update_traces(opacity=0.66, marker_line_width=0)
                    fig.update_layout(
                        paper_bgcolor="rgba(0,0,0,0)",
                        plot_bgcolor="rgba(0,0,0,0)",
                        margin=dict(l=20, r=20, t=58, b=20),
                        legend=dict(orientation="h", yanchor="bottom", y=1.07, xanchor="left", x=0, bgcolor="rgba(255,255,255,0.86)"),
                        font=dict(color=CPC_COLORS["ink"]),
                        showlegend=chart_data["state"].nunique() > 1,
                    )
                    fig.update_xaxes(gridcolor="rgba(221, 227, 234, 0.65)")
                    fig.update_yaxes(gridcolor="rgba(221, 227, 234, 0.65)")
                    st.plotly_chart(fig, use_container_width=True)
                    st.download_button(
                        tr("下载产物统计 CSV", "Download product statistics as CSV"),
                        data=chart_data.to_csv(index=False).encode("utf-8"),
                        file_name="molfragapp_product_statistics.csv",
                        mime="text/csv",
                        use_container_width=True,
                        key="download_product_stats",
                    )

                with c_chart2:
                    st.markdown(f"**{tr('能量描述符分布', 'Energy-descriptor distribution')}**")
                    with st.expander(tr("绘图参数", "Plot parameters"), expanded=False):
                        pk1, pk2 = st.columns(2)
                        ker_bin = pk1.number_input(tr("分箱数", "Bin count"), value=30, min_value=5, step=5, key="ker_bin")
                        ker_view = pk2.selectbox(
                            tr("显示模式", "View mode"),
                            options=["interactive", "publication"],
                            format_func=lambda x: tr("交互直方图", "Interactive histogram") if x == "interactive" else tr("论文风格图", "Publication-style figure"),
                            key="ker_view_mode",
                        )

                    ker_df = filtered.dropna(subset=["delta_energy"])
                    if len(ker_df) <= 1:
                        st.info(tr("当前筛选集不足以绘制 KER 分布。", "The current subset is too small for a KER distribution."))
                    elif ker_view == "interactive":
                        fig_ker = px.histogram(
                            ker_df,
                            x="delta_energy",
                            color="state",
                            nbins=int(ker_bin),
                            color_discrete_sequence=CPC_SEQUENCE,
                            labels={"delta_energy": descriptor_axis, "state": tr("电子态", "State")},
                        )
                        fig_ker.update_traces(opacity=0.64, marker_line_width=0)
                        fig_ker.update_layout(
                            paper_bgcolor="rgba(0,0,0,0)",
                            plot_bgcolor="rgba(0,0,0,0)",
                            margin=dict(l=20, r=20, t=58, b=20),
                            legend=dict(orientation="h", yanchor="bottom", y=1.06, xanchor="left", x=0, bgcolor="rgba(255,255,255,0.78)"),
                            showlegend=ker_df["state"].nunique() > 1,
                        )
                        st.plotly_chart(fig_ker, use_container_width=True)
                    else:
                        try:
                            fig_pub = plot_KER_state(ker_df, bin_count=int(ker_bin), kde_bandwidth=0.2)
                            fig_pub.patch.set_facecolor("#ffffff")
                            st.pyplot(fig_pub, use_container_width=True)
                        except Exception as e:
                            st.warning(tr(f"KER 论文风格图不可用：{e}", f"Publication-style KER plot is unavailable: {e}"))


                st.divider()
                plot_context = {
                    "statistics_view": {"ker_bin": int(ker_bin), "ker_view": ker_view},
                    "selected_record_path": selected_row.get("file") if selected_row is not None else None,
                }
                render_export_reuse_panel(
                    filtered,
                    selected_row,
                    h5_path=h5_path,
                    filter_context=filter_context,
                    schema_report=schema_report,
                    energy_meta=energy_meta,
                    plot_context=plot_context,
                )

        # === Tab 2: 轨迹检查 ===
        with tab2:
            if selected_row is None:
                st.info(tr("请先在左侧选择一条轨迹。", "Please select one trajectory on the left."))
            else:
                internal_path = selected_row["file"]

                with st.spinner(tr("读取轨迹索引与能量数据...", "Reading trajectory index and energy data...")):
                    total_frames_h5 = get_h5_frame_count(h5_path, internal_path)
                    expec_text = get_expec_text_from_h5(h5_path, internal_path)
                    energy_df = get_energy_dataframe_from_h5(h5_path, internal_path)

                if total_frames_h5 <= 0:
                    st.error(tr("无法读取该轨迹的 geometry 数据。", "Unable to read geometry data for this trajectory."))
                else:
                    c_mol, c_graph = st.columns([1, 1], gap="large")

                    with c_mol:
                        st.markdown(f"**{tr('结构演化', 'Structure evolution')}**")
                        ctrl1, ctrl2, ctrl3 = st.columns([1.3, 1.0, 1.0])
                        play_mode = ctrl1.radio(
                            tr("播放模式", "Playback mode"),
                            ["auto", "manual"],
                            format_func=lambda x: tr("自动循环", "Auto loop") if x == "auto" else tr("手动精查", "Manual inspection"),
                            horizontal=True,
                            label_visibility="collapsed",
                            key="play_mode",
                        )
                        time_step_fs = ctrl2.number_input(tr("步长 (fs)", "Time step (fs)"), value=0.5, step=0.1, key="viz_timestep")
                        viewer_size = ctrl3.selectbox(
                            tr("窗口尺寸", "Viewer size"),
                            options=["compact", "standard", "large"],
                            index=1,
                            format_func=lambda x: {
                                "compact": tr("紧凑", "Compact"),
                                "standard": tr("标准", "Standard"),
                                "large": tr("大图", "Large"),
                            }[x],
                            key="viewer_size",
                        )
                        size_map = {"compact": (560, 340), "standard": (700, 430), "large": (860, 540)}
                        viewer_width, viewer_height = size_map[viewer_size]

                        current_frame = 0
                        if play_mode == "manual":
                            current_frame = st.slider(
                                tr("进度", "Progress"),
                                0,
                                max(0, total_frames_h5 - 1),
                                0,
                                key="frame_slider",
                                label_visibility="collapsed",
                            )
                            cur_time = current_frame * time_step_fs
                            traj_col = get_traj_energy_column(energy_df) if not energy_df.empty else None
                            current_energy = "NA"
                            if traj_col and current_frame < len(energy_df):
                                try:
                                    current_energy = f"{float(energy_df[traj_col].iloc[current_frame]):.3f} eV"
                                except Exception:
                                    current_energy = "NA"

                            fm1, fm2, fm3 = st.columns(3)
                            fm1.metric(tr("当前帧", "Frame"), f"{current_frame} / {max(total_frames_h5 - 1, 0)}")
                            fm2.metric(tr("时间", "Time"), f"{cur_time:.1f} fs")
                            fm3.metric(tr("当前能量", "Current energy"), current_energy)
                            xyz_content = get_xyz_content_from_h5(h5_path, internal_path, frame_idx=current_frame)
                            show_mol_advanced(xyz_content, mode="manual", frame_idx=0, is_content=True, width=viewer_width, height=viewer_height)
                        else:
                            auto_stride = max(1, total_frames_h5 // 450)
                            xyz_content = get_xyz_content_from_h5(h5_path, internal_path, stride=auto_stride)
                            if auto_stride > 1:
                                st.caption(tr(f"预览模式：每 {auto_stride} 帧加载一帧。", f"Preview mode: one frame is loaded every {auto_stride} frames."))
                            show_mol_advanced(xyz_content, mode="auto", is_content=True, width=viewer_width, height=viewer_height)

                        cur_frame_xyz = get_xyz_content_from_h5(h5_path, internal_path, frame_idx=current_frame)
                        d1, d2 = st.columns(2)
                        with d1:
                            st.download_button(
                                tr("导出当前帧 XYZ", "Export current frame XYZ"),
                                data=cur_frame_xyz.encode("utf-8"),
                                file_name=f"{selected_row['traj_id']}_frame_{current_frame}.xyz",
                                mime="chemical/x-xyz",
                                use_container_width=True,
                                key="download_current_xyz",
                            )
                        with d2:
                            prepare_full_xyz = st.checkbox(
                                tr("准备完整轨迹导出", "Prepare full-trajectory export"),
                                value=False,
                                key="prepare_full_xyz_export",
                                help=tr("完整 XYZ 可能较大；勾选后才会读取全部帧。", "The full XYZ can be large; all frames are loaded only after enabling this option."),
                            )
                            if prepare_full_xyz:
                                full_xyz_content = get_xyz_content_from_h5(h5_path, internal_path, stride=1)
                                st.download_button(
                                    tr("导出完整轨迹 XYZ", "Export full trajectory XYZ"),
                                    data=full_xyz_content.encode("utf-8"),
                                    file_name=f"{selected_row['traj_id']}_trajectory.xyz",
                                    mime="chemical/x-xyz",
                                    use_container_width=True,
                                    key="download_full_xyz",
                                )

                        with st.expander(tr("GIF 动画导出", "GIF animation export"), expanded=False):
                            st.caption(
                                tr(
                                    "GIF 由 HDF5 geometry 帧渲染生成，用于快速保存结构演化示意；论文级三维渲染可继续使用外部可视化程序。",
                                    "The GIF is rendered from HDF5 geometry frames for quick trajectory snapshots; publication-grade 3D rendering can still be handled by external visualization tools.",
                                )
                            )
                            default_gif_end = max(0, total_frames_h5 - 1)
                            default_gif_stride = max(1, total_frames_h5 // 80) if total_frames_h5 > 0 else 1
                            g1, g2, g3, g4 = st.columns(4)
                            gif_start = g1.number_input(
                                tr("起始帧", "Start frame"),
                                min_value=0,
                                max_value=max(0, total_frames_h5 - 1),
                                value=0,
                                step=1,
                                key="gif_start_frame",
                            )
                            gif_end = g2.number_input(
                                tr("结束帧", "End frame"),
                                min_value=0,
                                max_value=max(0, total_frames_h5 - 1),
                                value=default_gif_end,
                                step=1,
                                key="gif_end_frame",
                            )
                            gif_stride = g3.number_input(
                                tr("帧间隔", "Frame stride"),
                                min_value=1,
                                max_value=max(1, total_frames_h5),
                                value=default_gif_stride,
                                step=1,
                                key="gif_stride",
                            )
                            gif_fps = g4.number_input(
                                tr("FPS", "FPS"),
                                min_value=1,
                                max_value=30,
                                value=12,
                                step=1,
                                key="gif_fps",
                            )
                            g5, g6 = st.columns([1, 1])
                            gif_size = g5.selectbox(
                                tr("图像尺寸", "Image size"),
                                options=[480, 640, 800],
                                index=1,
                                format_func=lambda x: f"{x} px",
                                key="gif_image_size",
                            )
                            show_gif_frame_label = g6.checkbox(
                                tr("显示帧号", "Show frame label"),
                                value=True,
                                key="gif_show_frame_label",
                            )
                            estimated_frames = len(range(min(int(gif_start), int(gif_end)), max(int(gif_start), int(gif_end)) + 1, max(1, int(gif_stride))))
                            if estimated_frames > 220:
                                st.warning(tr(f"当前设置预计生成约 {estimated_frames} 帧，可能较慢。建议增大帧间隔。", f"The current settings will render about {estimated_frames} frames and may be slow. Consider increasing the stride."))

                            if st.button(tr("生成 GIF", "Generate GIF"), type="primary", use_container_width=True, key="generate_traj_gif"):
                                with st.spinner(tr("正在生成 GIF 动画...", "Generating GIF animation...")):
                                    gif_bytes, rendered_frames = build_trajectory_gif(
                                        h5_path=h5_path,
                                        internal_path=internal_path,
                                        start_frame=int(gif_start),
                                        end_frame=int(gif_end),
                                        stride=int(gif_stride),
                                        fps=int(gif_fps),
                                        image_size=int(gif_size),
                                        show_frame_label=bool(show_gif_frame_label),
                                    )
                                if gif_bytes:
                                    st.session_state["trajectory_gif_bytes"] = gif_bytes
                                    st.session_state["trajectory_gif_name"] = f"{selected_row['traj_id']}_frames_{int(gif_start)}_{int(gif_end)}_stride_{int(gif_stride)}.gif"
                                    st.session_state["trajectory_gif_rendered_frames"] = rendered_frames
                                else:
                                    st.error(tr("GIF 生成失败：未读取到有效 geometry 帧。", "GIF generation failed: no valid geometry frames were read."))

                            gif_bytes_saved = st.session_state.get("trajectory_gif_bytes")
                            if gif_bytes_saved:
                                st.success(tr(f"GIF 已生成：{st.session_state.get('trajectory_gif_rendered_frames', '?')} 帧", f"GIF generated: {st.session_state.get('trajectory_gif_rendered_frames', '?')} frames"))
                                st.download_button(
                                    tr("下载 GIF 动画", "Download GIF animation"),
                                    data=gif_bytes_saved,
                                    file_name=st.session_state.get("trajectory_gif_name", f"{selected_row['traj_id']}_trajectory.gif"),
                                    mime="image/gif",
                                    use_container_width=True,
                                    key="download_trajectory_gif",
                                )

                    with c_graph:
                        st.markdown(f"**{tr('能量曲线', 'Energy curves')}**")
                        if not energy_df.empty:
                            highlight = current_frame if play_mode == "manual" else -1
                            show_energy_curve_interactive(
                                energy_df,
                                None,
                                is_content=True,
                                height=420,
                                highlight_idx=highlight,
                                key_suffix="detail_tab",
                            )
                        else:
                            st.warning(tr("无能量数据。", "No energy data."))

                    st.divider()
                    with st.expander(tr("内坐标计算工具", "Internal-coordinate tool"), expanded=False):
                        c1, c2, c3, c4 = st.columns(4)
                        with c1:
                            a1 = st.number_input("A1", value=0, min_value=0, key="s_a1")
                        with c2:
                            a2 = st.number_input("A2", value=0, min_value=0, key="s_a2")
                        with c3:
                            a3 = st.number_input("A3", value=0, min_value=0, key="s_a3")
                        with c4:
                            a4 = st.number_input("A4", value=0, min_value=0, key="s_a4")
                        atom_idxs = [i for i in [a1, a2, a3, a4] if i > 0]
                        if atom_idxs:
                            try:
                                xyz_for_incoord = get_xyz_content_from_h5(h5_path, internal_path, stride=1)
                                res = cal_incoord(xyz_for_incoord, atom_idxs, is_content=True)
                                if res:
                                    st.write(tr(f"计算结果: {res}", f"Result: {res}"))
                            except Exception as e:
                                st.warning(tr(f"内坐标计算失败：{e}", f"Internal-coordinate calculation failed: {e}"))

        # === Tab 3: 动量空间图 ===
        with tab3:
            render_analysis_note(
                tr('三体动量空间分析', 'Three-body momentum-space analysis'),
                tr('事件类别、碎片顺序和 analysis frame 会与生成的 Dalitz 图 / Newton 图一起保留，便于从图像反查轨迹。', 'Event class, fragment order, and analysis frame are kept with Dalitz plots / Newton diagrams for trajectory lookup.'),
                accent='purple',
            )
            source_scope = st.radio(
                tr("事件来源", "Event source"),
                options=["current_filtered_subset", "all_records"],
                index=0,
                horizontal=True,
                format_func=lambda x: tr("当前筛选集", "Current filtered subset") if x == "current_filtered_subset" else tr("全部记录", "All records"),
                key="dn_source_scope",
                help=tr("默认只从左侧当前筛选集生成动量空间图，保证统计、轨迹表和动量图共用同一事件集合。", "By default, momentum maps are generated only from the current filtered subset so statistics, the table, and maps share the same event set."),
            )
            df_dn_base = filtered.copy() if source_scope == "current_filtered_subset" else trajs_data.copy()
            all_3frag = sorted(set([f for f in df_dn_base.frag_string.astype(str) if len(f.split(" + ")) == 3]))
            state_options_dn = sorted(set(df_dn_base.state.astype(str)))
            if "dn_frag" in st.session_state:
                st.session_state["dn_frag"] = [x for x in st.session_state["dn_frag"] if x in all_3frag]
            if "dn_state" in st.session_state:
                st.session_state["dn_state"] = [x for x in st.session_state["dn_state"] if x in state_options_dn]

            with st.expander(tr("1. 选择事件类别", "1. Select event class"), expanded=True):
                d1, d2, d3 = st.columns(3)
                with d1:
                    frag_dn = st.multiselect(tr("三体产物", "Three-body products"), options=all_3frag, default=all_3frag, key="dn_frag")
                with d2:
                    state_dn = st.multiselect(tr("电子态", "Electronic state"), options=state_options_dn, default=state_options_dn, key="dn_state")
                with d3:
                    min_steps_dn = st.number_input(tr("最小步数", "Minimum steps"), value=1, min_value=1, key="dn_step")

            df_dn = df_dn_base[
                (df_dn_base.nstep >= min_steps_dn)
                & (df_dn_base.nfrag == 3)
                & (df_dn_base.state.isin(state_dn))
                & (df_dn_base.frag_string.isin(frag_dn))
            ].reset_index(drop=True)
            st.caption(tr(f"可用三体轨迹: {len(df_dn)}；来源: 当前筛选集" if source_scope == "current_filtered_subset" else f"可用三体轨迹: {len(df_dn)}；来源: 全部记录", f"Available three-body trajectories: {len(df_dn)}; source: current filtered subset" if source_scope == "current_filtered_subset" else f"Available three-body trajectories: {len(df_dn)}; source: all records"))

            with st.expander(tr("2. 定义动量映射", "2. Define momentum mapping"), expanded=True):
                st.caption(tr("A/B/C 对应 HDF5 中 frag_list 的三体碎片顺序。", "A/B/C refer to the three-fragment order stored in frag_list."))
                p1, p2, p3, p4 = st.columns([1, 1, 1, 1.2])
                frag_a = p1.selectbox("A", options=[1, 2, 3], index=2, key="dn_order_a")
                frag_b = p2.selectbox("B", options=[1, 2, 3], index=1, key="dn_order_b")
                frag_c = p3.selectbox("C", options=[1, 2, 3], index=0, key="dn_order_c")
                idx_frame = p4.number_input(tr("帧序号 (-1=末态)", "Frame index (-1 = final)"), value=-1, key="dn_frame")
                atom_orders = [frag_a, frag_b, frag_c]

            run_btn = st.button(tr("生成 Dalitz 图 / Newton 图", "Generate Dalitz plot / Newton diagram"), type="primary", use_container_width=True, key="dn_btn")

            if run_btn:
                if df_dn.empty:
                    st.error(tr("当前事件类别下没有可用三体轨迹。", "No three-body trajectories are available for the current event class."))
                elif len(set(atom_orders)) != 3:
                    st.error(tr("A、B、C 必须对应三个不同碎片。", "A, B, and C must correspond to three different fragments."))
                else:
                    with st.spinner(tr("计算中...", "Calculating...")):
                        try:
                            orders_adjusted = [x - 1 for x in atom_orders]
                            records_json = json.dumps(to_jsonable(df_dn.to_dict(orient="records")), ensure_ascii=False)
                            orders_json = json.dumps(orders_adjusted)
                            three_body_table = build_three_body_observable_table(h5_path, records_json, orders_json, int(idx_frame))
                            dn = dnplot(traj_ids=df_dn["file"].tolist(), orders=orders_adjusted, idx=int(idx_frame), hdf5_path=h5_path)
                            pc1, pc2 = st.columns(2, gap="large")
                            with pc1:
                                st.markdown(f"**{tr('Dalitz 图', 'Dalitz plot')}**")
                                fig_dalitz = normalize_matplotlib_figure(dn[0], size=(4.8, 4.8))
                                st.pyplot(fig_dalitz, use_container_width=True)
                            with pc2:
                                st.markdown(f"**{tr('Newton 图', 'Newton diagram')}**")
                                fig_newton = normalize_matplotlib_figure(dn[1], size=(4.8, 4.8))
                                st.pyplot(fig_newton, use_container_width=True)
                            st.download_button(
                                tr("下载三体 observable 表 CSV", "Download three-body observable table CSV"),
                                data=three_body_table.to_csv(index=False).encode("utf-8"),
                                file_name="molfragapp_three_body_observable_table.csv",
                                mime="text/csv",
                                use_container_width=True,
                                key="download_three_body_observable_table",
                            )
                            with st.expander(tr("三体 observable 表预览", "Three-body observable table preview"), expanded=False):
                                st.dataframe(three_body_table, use_container_width=True, hide_index=True)
                        except Exception as e:
                            st.error(tr(f"动量空间图生成失败：{e}", f"Momentum-space map generation failed: {e}"))

        # === Tab 4: 元数据 ===
        with tab4:
            render_analysis_note(
                tr('Metadata audit', 'Metadata audit'),
                tr('保留 HDF5 internal path、record identity、screening metadata 和 datasets 信息，用于结果复核和脚本复用。', 'Preserve HDF5 internal path, record identity, screening metadata, and dataset information for audit and reuse.'),
                accent='green',
            )
            if selected_row is None:
                st.info(tr("请选择一条轨迹。", "Please select one trajectory."))
            else:
                overview = get_h5_group_overview(h5_path, selected_row["file"])
                st.markdown(f"**{tr('File schema', 'File schema')}**")
                schema_df = pd.DataFrame(
                    [
                        [tr("schema status", "schema status"), schema_report.get("status", "unknown")],
                        [tr("trajectory groups", "trajectory groups"), schema_report.get("trajectory_group_count", 0)],
                        [tr("missing recommended file attrs", "missing recommended file attrs"), ", ".join(schema_report.get("missing_recommended_file_attributes", [])) or "None"],
                    ],
                    columns=[tr("字段", "Field"), tr("值", "Value")],
                )
                st.dataframe(schema_df, use_container_width=True, hide_index=True, height=150)

                st.markdown(f"**{tr('能量描述符', 'Energy descriptor')}**")
                energy_descriptor_df = pd.DataFrame(list(summarize_energy_descriptor(energy_meta).items()), columns=[tr("字段", "Field"), tr("值", "Value")])
                st.dataframe(energy_descriptor_df, use_container_width=True, hide_index=True, height=150)

                st.markdown(f"**{tr('轨迹摘要', 'Trajectory summary')}**")
                summary_df = pd.DataFrame(
                    [
                        [tr("轨迹 ID", "Trajectory ID"), selected_row.get("traj_id", "Unknown")],
                        [tr("HDF5 内部路径", "HDF5 internal path"), selected_row.get("file", "")],
                        [tr("电子态", "Electronic state"), selected_row.get("state", "Unknown")],
                        [tr("步数", "Number of steps"), selected_row.get("nstep", 0)],
                        [tr("片段数", "Fragment count"), selected_row.get("nfrag", 0)],
                        [tr("产物通道", "Product channel"), selected_row.get("frag_string", "Bound")],
                        [descriptor_axis, format_ker(selected_row.get("delta_energy"))],
                    ],
                    columns=[tr("字段", "Field"), tr("值", "Value")],
                )
                st.dataframe(summary_df, use_container_width=True, hide_index=True, height=260)

                st.markdown(f"**{tr('可用数据集', 'Available datasets')}**")
                if overview["datasets"]:
                    st.dataframe(pd.DataFrame(overview["datasets"]), use_container_width=True, hide_index=True)
                else:
                    st.caption(tr("没有检测到子数据集。", "No child datasets detected."))

                with st.expander(tr("高级：原始属性 JSON", "Advanced: raw attribute JSON"), expanded=False):
                    raw_dict = to_jsonable(selected_row.to_dict())
                    raw_dict["hdf5_attributes"] = overview.get("attributes", [])
                    raw_dict["schema_report"] = to_jsonable(schema_report)
                    raw_dict["energy_descriptor"] = to_jsonable(energy_meta)
                    st.json(raw_dict)
                    st.download_button(
                        tr("下载元数据 JSON", "Download metadata JSON"),
                        data=json.dumps(raw_dict, ensure_ascii=False, indent=2).encode("utf-8"),
                        file_name=f"{selected_row['traj_id']}_metadata.json",
                        mime="application/json",
                        use_container_width=True,
                        key="download_metadata_json",
                    )

    st.markdown("---")
    st.markdown(
        f"<div style='text-align: center; color: #999; font-size: 0.8rem;'>© 2026 MolFragApp | {tr('由 Streamlit 驱动', 'Powered by Streamlit')}</div>",
        unsafe_allow_html=True,
    )


if __name__ == "__main__":
    main()
