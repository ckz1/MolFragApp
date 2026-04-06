"""Web UI for MolFragApp
=========================

Streamlit 基础
---------------

* 基础：API与设置

    * `Streamlit API cheat sheet <https://docs.streamlit.io/develop/quick-reference/cheat-sheet>`_

        * `API reference <https://docs.streamlit.io/develop/api-reference>`_

        * `Configuration <https://docs.streamlit.io/develop/api-reference/configuration>`_

* 进阶：架构与执行流程

    * `Caching and state <https://docs.streamlit.io/develop/api-reference/caching-and-state>`_

        * `Working with Streamlit's execution model <https://docs.streamlit.io/develop/concepts/architecture>`_

        * `Execution flow <https://docs.streamlit.io/develop/api-reference/execution-flow>`_

        * 默认的执行流程是执行整个页面， 在 `st.fragment` 内交互时，只会重新运行 `st.fragment`，不会运行整个程序。可以按照指定的时间间隔 `rerun`

* 联合：

    * Plotly in Streamlit 可以返回选中的数据

        * 列表选择 `event.selection` 可以在 `st.fragment` 中使用

架构
-----

* `MolFragApp`
        #. 获取分析参数
        #. 根据分析参数对数据文件进行更新( `data` , `xyztraj` )
        #. 读取数据文件得到数据( `data` )
        #. 对数据进行不同角度的分析与可视化( `utils` , `py3bodyfrag` )
* `data`
        #. 轨迹数据的提取
        #. 数据文件的处理
* `xyztraj`
        #. xyz文件的读取与处理
* `py3bodyfrag`
        #. Dalitz Plot 和 Newton Plot 的绘制
* `utils`
        #. 存储 `st.fragment` 等 `Web UI` 部件
        #. 可视化与统计的实现

建议
-----

#. 尝试使用 `st.session_state` 处理共享的数据
#. 数据尽量集中，函数尽量封闭。尽可能明确函数所需要的参数，以及各个过程输入输出的参数。函数尽量使用显式声明所需的参数，没有声明但是调用的参数要在 `Note` 中声明
#. 不需要重新运行全部脚本的部分封装成函数，并使用 `st.fragment` 装饰器
#. 不需要立即更新的输入使用 `st.form()`
#. 适当添加提示信息 `st.toast()` 和 `st.spinner()`
#. 支持 `CMD` 和 `Web UI` 多种方式调用。

使用说明
---------

启动网页时会按照默认分析参数提取数据，如果默认分析与原来的参数不一致，原来的数据文件将被重命名(只有一次)，新的数据文件会被创建。如果需要保留原始的数据文件应该启用 `readonly` 模式，此时只会检查文件是否存在并不会对数据进行更新。

TODO

    多原子片段的Dalitz Plot 和Newton Plot

    模拟设置文件 和 模板文件需要修改为根据 数据文件的上上级目录检查

    文件路径  绝对 vs 相对
"""

import os
from glob import glob
import numpy as np
import pandas as pd
import streamlit as st

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
from glob import glob
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import colorsys
import matplotlib as mpl

# import xyztraj
import py3bodyfrag as ptf
import py3Dmol

# import h5py

# 自定义模块
import sharc
# import utils

import sys
import logging

logger = logging.getLogger(__package__)
logger.setLevel(logging.INFO)

stdout_handler = logging.StreamHandler(stream=sys.stdout)
stdout_handler.setLevel(
    logging.INFO
)  # CRITICAL > ERROR > WARNING > INFO > DEBUG > NOTSET
stdout_handler.setFormatter(
    logging.Formatter(
        fmt="[%(asctime)s] %(levelname)s: %(message)s",
        # datefmt='%H:%M:%S',
    )
)
logger.addHandler(stdout_handler)

__APP_NAME = "MolFragApp"
__AUTHORS = "Chenkai Zhang"
__VERSION = "2025.09.30_alpha"
__LOGO_PATH = "logo.png"
__GITHUB_LINK = "https://github.com/ckz1/MolFragApp"
__DOCS_LINK = "https://molfragapp.readthedocs.io/en/latest/"
__ISSUES_LINK = "https://github.com/ckz1/MolFragApp/issues"
__ABOUT_APP = f"**MolFragApp** is a project based on Python language that uses molecular structure files in `xyz` format for fragmentation analysis and visualization of trajectories.\n > Version: {__VERSION} [GitHub]({__GITHUB_LINK}) [Issues]({__ISSUES_LINK})"

## ==============================================
## 参数设置
## ==============================================
HDF5_FILE = "Data.hdf5"  # str: 数据文件名
DEFAULT_DAT_PATTERN = "*let_*/TRAJ_*/output.dat"  # str: 默认分析参数 -- dat文件路径模式
DEFAULT_MAX_BOND_LENGTH = 2.5  # float: 默认分析参数 -- 最大键长
READONLY_MODE = True  # 只读模式，用于在无原始数据的环境下读取数据文件，不会更新数据，同时支持的功能有限
## ==============================================


st.set_page_config(
    # 设置浏览器标签页
    page_icon="🔍",
    page_title=f"{__APP_NAME}.v{__VERSION}",
    layout="wide",
    # 右上角选项设置
    menu_items={
        "Get Help": __DOCS_LINK,
        "Report a bug": __ISSUES_LINK,
        "About": __ABOUT_APP,
    },
)

# logo与版权信息设置
# ===================

ASCII_LOGO = r"""
    __  ___           __    ______                             ___                   
   /  |/  /  ____    / /   / ____/   _____  ____ _   ____ _   /   |    ____     ____ 
  / /|_/ /  / __ \  / /   / /_      / ___/ / __ `/  / __ `/  / /| |   / __ \   / __ \
 / /  / /  / /_/ / / /   / __/     / /    / /_/ /  / /_/ /  / ___ |  / /_/ /  / /_/ /
/_/  /_/   \____/ /_/   /_/       /_/     \__,_/   \__, /  /_/  |_| / .___/  / .___/ 
                                                  /____/           /_/      /_/      

"""


@st.fragment
def display_logo():
    """显示logo

    Note:
        __LOGO_PATH (str): logo文件路径
        __GITHUB_LINK (str): GitHub 网址
    """
    try:
        file_path = os.path.realpath(__file__)
        logo_path = os.path.join(os.path.dirname(file_path), __LOGO_PATH)
        st.logo(logo_path, size="large", link=__GITHUB_LINK)
    except Exception as e:
        logger.error(f"Failed to display the logo.\n{e}")


@st.fragment
def display_copyright():
    """显示版权信息

    Note:
        __GITHUB_LINK (str): GitHub 网址
        __APP_NAME (str): 应用名称
        __AUTHORS (str): 作者名字
    """
    import datetime

    start_year = 2024
    end_year = datetime.datetime.today().year

    st.divider()

    st.components.v1.html(
        f"""
    <p align='center'>
        <a href="{__GITHUB_LINK}" target="_blank">
            <img border="0" src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png" alt="{__APP_NAME} in Github" width="30">
        </a>
    </p>
    <p align='center'> Copyright © {start_year}-{end_year} {__AUTHORS}. All rights reserved. </p>
    """
    )


# Define the rgbcolor class
class rgbcolor:
    def __init__(self, initlist):
        excluded = [[0.12, 0.22]]  # exclude yellow hues from the colorwheel
        excluded.sort(key=lambda x: x[0])
        temp1 = [
            [min(1.0, max(0.0, el[0])), max(0.0, min(1.0, el[1]))] for el in excluded
        ]
        temp2 = [[0.0, 0.0]]
        for el in temp1:
            if el[0] >= temp2[-1][1]:
                temp2.append(el)
            else:
                temp2[-1][1] = el[1]
        self.excluded = temp2[1:]
        self.initlist = [max(0, el) for el in initlist]
        self.n = sum(el > 0 for el in self.initlist)
        self.m = len(initlist)
        self.a = 1.0 - sum(el[1] - el[0] for el in self.excluded)
        self.startlist = [0.0] * self.m
        self.incrlist = [0.0] * self.m
        for i in range(1, self.m):
            self.startlist[i] = self.startlist[i - 1] + (
                self.a / self.n if self.initlist[i - 1] > 0 else 0
            )
        for i in range(self.m):
            if self.initlist[i] > 0:
                self.incrlist[i] = self.a / self.n / self.initlist[i]

    def rgb_to_hex(self, rgb):
        return "#{:02x}{:02x}{:02x}".format(
            int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
        )

    def hexcolor(self, index, el):
        if not (1 <= index <= self.m and 1 <= el <= self.initlist[index - 1]):
            return "#FFFFFF"
        hue = self.startlist[index - 1] + self.incrlist[index - 1] * (el - 1)
        for start, end in self.excluded:
            if hue > start:
                hue += end - start
        return self.rgb_to_hex(colorsys.hsv_to_rgb(hue, 1, 1))


def generate_colors(num: int = 10, alpha: float = 0.4, rgba_format: str = "plotly"):
    """

    Note:
        rgba_format
            * plotly 'rgba(255,255,255,0.5)'
            * matplotlib (0.1,0.5,0.9,0.5)
    """

    def hex2rgb(hex_str):
        # return tuple(int(hex_str.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
        return list(int(hex_str.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))

    def rgba2hex(rgb_tuple, alpha):
        rgba_tuple = (rgb_tuple[0] / 255, rgb_tuple[1] / 255, rgb_tuple[2] / 255, alpha)
        return mpl.colors.rgb2hex(rgba_tuple, keep_alpha=True)

    color_generator = rgbcolor([1] * num)
    colors_generated = [color_generator.hexcolor(i + 1, 1) for i in range(num)]
    # return [rgba2hex(hex2rgb(c),alpha) for c in colors_generated]
    if rgba_format == "plotly":
        return ["rgba" + str(tuple(hex2rgb(c) + [alpha])) for c in colors_generated]
    elif rgba_format == "matplotlib":
        return [
            tuple([value / 255 for value in hex2rgb(c)] + [alpha])
            for c in colors_generated
        ]
    else:
        raise ValueError("Wrong rgba_format")


# 获取分析参数相关
# =================


def read_template(filename="MOLCAS.template", output_type="str"):
    """读取 `SHARC` 的 `template` 文件

    Args:
        filename (str, optional): `template` 文件名. Defaults to "MOLCAS.template".
        output_type (str, optional): 输出类型: `str` or `dict`. Defaults to "str".

    Returns:
        object: `template` 文件内容，类型由 `output_type` 指定
    """
    assert os.path.exists(filename)

    if "MOLCAS.template" in filename:
        qc_software = "OpenMolcas"
    elif "MOLPRO.template" in filename:
        qc_software = "Molpro"
    else:
        qc_software = "Unknown"

    with open(filename, "r") as f:
        lines = f.readlines()
    f.close()

    template = {"qc_software": qc_software, "method": "CASSCF(default)"}
    for line in lines:
        temp = line.split()
        if len(temp) >= 2 and "#" not in line:
            template[temp[0]] = temp[1]

    if output_type == "str":
        if qc_software == "OpenMolcas":
            return f"SA({template['roots']})-{template['method']}({template['nactel']}e,{template['ras2']}o)/{template['basis']}, {template['qc_software']}"
        elif qc_software == "Molpro":
            return f"SA({template['roots']})-{template['method']}({int(template['nelec'])-int(template['closed'])*2}e,{int(template['occ'])-int(template['closed'])}o)/{template['basis']}, {template['qc_software']}"
    else:
        return template


def get_analysis_params() -> dict:
    """获取分析参数

    Return:
        dict: 分析参数 `dict(dat_pattern=dat_pattern, max_bond_length=max_bond_length, max_nstep=max_nstep)` ,
            其中最重要的 `dat_pattern` , `max_bond_length` 用于更新数据和筛选数据, `max_nstep` 用于筛选数据

    Note:
        两个默认的分析参数由 `DEFAULT_DAT_PATTERN` , `DEFAULT_MAX_BOND_LENGTH` 设定
    """
    st.write(f"当前路径:`{os.getcwd()}`")

    params = st.columns(4)

    # 分析参数: 文件名模式和最大键长
    with params[0]:
        dat_pattern = st.text_input(
            "数据文件(`dat`格式)路径:",
            # value=analysis_params['dat_pattern'],
            value=DEFAULT_DAT_PATTERN,
            help="使用提取到的最后一帧的结构分析解离片段的成分",
        )
        num_datfile = len(glob(dat_pattern))

        if num_datfile == 0:
            st.warning(
                (
                    f"ReadOnly mode on, `{HDF5_FILE}` used."
                    if READONLY_MODE
                    else "No dat files found!"
                ),
                icon="⚠️",
            )
        else:
            st.write("轨迹文件总数：", num_datfile)

    with params[1]:
        max_bond_length = st.number_input(
            "划分片段时的最大键长(Å)",
            # value=analysis_params['max_bond_length'],
            value=DEFAULT_MAX_BOND_LENGTH,
            min_value=0.1,
            step=0.1,
            help="若两原子之间距离大于最大键长，则被划分到不同片段",
        )

    with params[2]:
        # 辅助信息 -- 轨迹步数设置
        setup_traj_path = st.text_input(
            "`KEYSTROKES.setup_traj`文件路径:",
            value="KEYSTROKES.setup_traj",
        )
        max_nstep = 0
        try:
            max_nstep = get_max_nstep(setup_traj_path)
            st.write(f"轨迹步数设置：`{max_nstep}`")
        except Exception as e:
            st.warning(f"Can NOT read file `{setup_traj_path}` !\n{e}", icon="⚠️")

    with params[3]:
        # 辅助信息 -- 计算级别
        template_path = st.text_input(
            "`template`文件路径:",
            value="*.template",
        )
        # st.write(read_template(glob(template_path)[0]))
        try:
            cal_method = read_template(glob(template_path)[0])
            st.write(f"计算级别:`{cal_method}`")
        except Exception as e:
            st.warning(f"Can NOT read file `{template_path}` !\n{e}", icon="⚠️")

    return dict(
        dat_pattern=dat_pattern, max_bond_length=max_bond_length, max_nstep=max_nstep
    )


# 统计、计算、可视化功能
# ========================


@st.fragment
def counts_stats(count, name="steps", nbins=20):
    with st.popover("轨迹长度分布"):
        # fig = ff.create_distplot(
        # [count],
        # [name],
        # show_curve=False,
        # bin_size=(max(count) - min(count)) / nbins,
        # )
        fig = px.histogram(pd.DataFrame({name: count}), x=name, nbins=nbins)
        # return fig
        st.plotly_chart(fig)


def cal_incoord(traj, timestep, *idxs):
    """计算内坐标"""
    time = None
    incoord = None
    incoord_name = ""
    try:

        idxs = idxs[0]
        # st.write(idxs)  # 数据类似 ([1,2,3],)

        # traj = xyztraj.Traj(xyzfile)
        # traj.parse_title(lambda x: float(x.split()[1]))
        # time = traj.title_info
        time = np.arange(len(traj)) * timestep
        natom_selected = len(idxs)
        # st.write(natom_selected)
        # st.write(idxs)

        if natom_selected == 2:
            incoord = traj.bond(*idxs)
            incoord_name = "Bond %d-%d (Å)" % tuple(idxs)
        elif natom_selected == 3:
            incoord = traj.angle(*idxs)
            incoord_name = "Angle %d-%d-%d (°)" % tuple(idxs)
        elif natom_selected == 4:
            incoord = traj.dihedral(*idxs)
            incoord_name = "Dihedral %d-%d-%d-%d (°)" % tuple(idxs)
        else:
            pass

    except Exception as e:
        # print(f"Failed to calculate internal coordinate of {xyzfile}.")
        logger.error(f"Failed to calculate internal coordinate of {traj}. {e}")
    return time, incoord, incoord_name


def show_mol(xyz: str = "geom.xyz", width: int = 400, height: int = 400) -> None:
    """播放多帧 `xyz` 文件的动画

    Args:
        xyz (str, optional): `xyz` 文件名. Defaults to "geom.xyz".
        width (int, optional): 动画宽度(px). Defaults to 400.
        height (int, optional): 动画高度(px). Defaults to 400.
    """
    if os.path.exists(xyz):
        with open(xyz, "r") as f:
            mol_block = f.read()
        f.close()
    else:
        mol_block = xyz

    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModelsAsFrames(mol_block, "xyz")
    viewer.setStyle(
        {
            "stick": {"radius": 0.1, "colorscheme": "Jmol"},
            "sphere": {"radius": 0.3, "colorscheme": "Jmol"},
        }
    )
    viewer.setBackgroundColor("black")
    viewer.animate(
        {
            "loop": "forward",
        }
    )
    viewer.zoomTo()

    components.html(viewer._make_html(), width=width, height=height)
    # st.write('缩放: 滚动鼠标滑轮, 旋转: 拖动鼠标左键, 平移: 拖动鼠标滑轮')


def show_energy_curve(
    expec_out_df, incoord_data: tuple = None, width: int = 400, height: int = 400
):
    try:
        data = expec_out_df
        # 选择哪些列需要画图
        with st.form("data_to_plot"):
            options = st.multiselect(
                label="选择需要画图的数据",
                options=list(data.columns),
                # default = ['Trajectory'],
                placeholder="可选择多个",
                # help = 'X轴为Time, 默认画Trajectory',
            )
            submitted = st.form_submit_button("确认")

        # 对于内坐标曲线使用双y轴，否则使用单y轴
        if incoord_data is None:
            fig = go.Figure()
        else:
            # Create figure with secondary y-axis
            fig = make_subplots(specs=[[{"secondary_y": True}]])
            fig.add_trace(
                go.Scatter(x=incoord_data[0], y=incoord_data[1], name=incoord_data[2]),
                secondary_y=True,
            )
        # 绘制所模拟的轨迹的势能变化
        fig.add_trace(
            go.Scatter(
                x=data["Time"],
                y=data["Epot"],
                mode="markers",
                line=dict(color="black"),
                name="Trajectory",
            )
        )
        # 绘制所选数据的曲线
        for key in options:
            fig.add_trace(
                go.Scatter(x=data["Time"], y=data[key], mode="lines", name=key)
            )
        # 设置轴标题
        fig.update_layout(
            xaxis=dict(title="Time (fs)"),
            yaxis=dict(title="Energies in diagonal basis (eV)"),
        )
        st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.warning(f"Failed to show energy curve. {e}")


@st.fragment
def show_traj(traj: str, expec_out_df: str, frag_info: str = None) -> None:
    # 调用了 函数 show_mol show_energy_curve
    timestep = (expec_out_df["Time"].iloc[-1] - expec_out_df["Time"].iloc[0]) / (
        len(expec_out_df["Time"]) - 1
    )
    incoord_data = None
    # 选择需要计算内坐标的原子序号
    fraginfo_on = st.toggle("显示片段信息")
    if fraginfo_on:
        # st.write(st.session_state.df.iloc[idx_row]["frag"])
        st.write(frag_info)
        with st.form("atom_to_cal"):
            atom_idxs = st.multiselect(
                label="选择原子序号计算内坐标(原子序号从1开始) ",
                options=[i + 1 for i in range(traj.natom)],
                default=[],
            )
            submitted = st.form_submit_button("确认")
        if len(atom_idxs) > 0:
            # st.write(atom_idxs)
            incoord_data = cal_incoord(traj, timestep, atom_idxs)

    # show_mol(traj.xyz("traj.xyz"), width=850, height=350)
    show_mol(traj.xyz(), width=850, height=350)
    # st.write(incoord_data)
    # st.write(traj.bond(1,2))
    show_energy_curve(expec_out_df, incoord_data)


# f = h5py.File("Data.hdf5", "r")
# expec_out_df = pd.DataFrame(data = f["State_1/TRAJ_00074"]["expec_data"][:], columns = f["State_1/TRAJ_00074"].attrs["expec_columns"])
# frag_info = {0: {'atoms': ['H'], 'index': [0], 'formula': 'H'}, 1: {'atoms': ['H', 'O'], 'index': [1, 2], 'formula': 'HO'}}
# traj = xyztraj.Traj(natom=3, atoms=["H","H","O"], coord=f["State_1/TRAJ_00074"]["geometry"][:]*xyztraj.BOHR_IN_ANGSTROM)
# show_traj(traj, expec_out_df, frag_info)

# ================================================================
# KER
# def plot_KER(df):


def plot_KER_state(df, bin_count=30, kde_bandwidth=0.2):
    # 用户可配置参数
    # plot_range = (0, 16)
    # bin_count = 30
    # kde_bandwidth = 0.2
    opacity = 0.3

    # 加载数据
    # df = pd.read_csv('2024-11-15T08-01_export.csv')
    plot_range = (df["delta_energy"].min(), df["delta_energy"].max())

    color_generator = rgbcolor([1] * len(df["state"].unique()))
    state_colors = [
        color_generator.hexcolor(i + 1, 1) for i in range(len(df["state"].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [0.7, 0.3]}, sharex=True
    )

    y_kde_sum = np.zeros(1000)
    x_range = np.linspace(plot_range[0], plot_range[1], 1000)
    hist_max = 0

    for state, color in zip(df["state"].unique(), state_colors):
        subset = df[df["state"] == state]["delta_energy"]
        hist, edges = np.histogram(
            subset, bins=np.linspace(plot_range[0], plot_range[1], bin_count)
        )
        ax2.bar(
            edges[:-1],
            hist,
            width=np.diff(edges)[0],
            alpha=opacity,
            color=color,
            label=state,
        )
        if max(hist) > hist_max:
            hist_max = max(hist)

        kde = gaussian_kde(subset, bw_method=kde_bandwidth)
        y_kde = kde.evaluate(x_range)
        ax1.fill_between(x_range, 0, y_kde, alpha=opacity, color=color, label=state)
        y_kde_sum += y_kde

    # 归一化并绘制总和曲线
    y_kde_sum_normalized = y_kde_sum / y_kde_sum.max()
    ax1.plot(
        x_range,
        y_kde_sum_normalized,
        color="black",
        label="Normalized Total Sum",
        linewidth=2,
    )

    # Plot exp data for test
    # data_exp = pd.read_csv("data_exp.csv")
    # ax1.plot(
    # data_exp['X1'],
    # data_exp['Y1']/data_exp['Y1'].max(),
    # 'r.-'
    # )

    # ax1.plot(
    # data_exp['X2'],
    # data_exp['Y2']/data_exp['Y2'].max(),
    # 'b.-'
    # )

    # 设置图表的整体属性
    ax1.set_ylabel("Intensity (arb. units)")
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False, fontsize=8)

    ax2.set_ylabel("Count")
    ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel("KER (eV)")
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig


def plot_KER_state2(df, sigma, alpha=0.3, num_point=1000):

    def gaussian(x, sigma=1, mu=0):
        gx = np.exp(-((x - mu) ** 2) / sigma**2 / 2) / (sigma * np.sqrt(2 * np.pi))
        return gx

    plot_range = (df["delta_energy"].min(), df["delta_energy"].max())

    color_generator = rgbcolor([1] * len(df["state"].unique()))
    generated_colors = [
        color_generator.hexcolor(i + 1, 1) for i in range(len(df["state"].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [0.7, 0.3]}, sharex=True
    )

    x_range = np.linspace(plot_range[0], plot_range[1], num_point)
    y_sum = np.zeros(num_point)
    plot_data = dict()

    for plot_type in df["state"].unique():
        y = np.zeros(num_point)
        temp_delta_energy = df[df["state"] == plot_type].delta_energy
        for tde in temp_delta_energy:
            y += gaussian(x_range, sigma=sigma, mu=tde)
        y_sum += y
        plot_data[plot_type] = dict(curve_data=y, point_data=temp_delta_energy)

    counter = 0
    ratio = max(y_sum)
    for plot_type in plot_data:
        ax1.fill_between(
            x_range,
            0,
            plot_data[plot_type]["curve_data"] / ratio,
            alpha=alpha,
            color=generated_colors[counter],
            label=plot_type,
        )
        ax2.plot(
            plot_data[plot_type]["point_data"],
            np.zeros(len(plot_data[plot_type]["point_data"])) - counter,
            marker="|",
            lw=0,
            color=generated_colors[counter],
            alpha=alpha,
        )
        counter += 1
    # 归一化并绘制总和曲线

    ax1.plot(
        x_range, y_sum / ratio, color="black", label="Normalized Total Sum", linewidth=2
    )

    # 设置图表的整体属性
    ax1.set_ylabel("Intensity (arb. units)")
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False, fontsize=8)

    # ax2.set_ylabel('Count')
    ax2.set_yticks([])
    # ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel("KER (eV)")
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig


def plot_KER_frag(df, bin_count=30, kde_bandwidth=0.2):
    # 用户可配置参数
    # plot_range = (0, 16)
    # bin_count = 30
    # kde_bandwidth = 0.2
    opacity = 0.3

    # 加载数据
    # df = pd.read_csv('2024-11-15T08-01_export.csv')
    plot_range = (df["delta_energy"].min(), df["delta_energy"].max())

    color_generator = rgbcolor([1] * len(df["frag_string"].unique()))
    frag_string_colors = [
        color_generator.hexcolor(i + 1, 1)
        for i in range(len(df["frag_string"].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [0.7, 0.3]}, sharex=True
    )

    y_kde_sum = np.zeros(1000)
    x_range = np.linspace(plot_range[0], plot_range[1], 1000)
    hist_max = 0

    for frag_string, color in zip(df["frag_string"].unique(), frag_string_colors):
        subset = df[df["frag_string"] == frag_string]["delta_energy"]
        hist, edges = np.histogram(
            subset, bins=np.linspace(plot_range[0], plot_range[1], bin_count)
        )
        ax2.bar(
            edges[:-1],
            hist,
            width=np.diff(edges)[0],
            alpha=opacity,
            color=color,
            label=frag_string,
        )
        if max(hist) > hist_max:
            hist_max = max(hist)

        kde = gaussian_kde(subset, bw_method=kde_bandwidth)
        y_kde = kde.evaluate(x_range)
        ax1.fill_between(
            x_range, 0, y_kde, alpha=opacity, color=color, label=frag_string
        )
        y_kde_sum += y_kde

    # 归一化并绘制总和曲线
    y_kde_sum_normalized = y_kde_sum / y_kde_sum.max()
    ax1.plot(
        x_range,
        y_kde_sum_normalized,
        color="black",
        label="Normalized Total Sum",
        linewidth=2,
    )

    # 设置图表的整体属性
    ax1.set_ylabel("Intensity (arb. units)")
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False, fontsize=8)

    ax2.set_ylabel("Count")
    ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel("KER (eV)")
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig


def plot_KER_state2(df, sigma, alpha=0.3, num_point=1000):

    def gaussian(x, sigma=1, mu=0):
        gx = np.exp(-((x - mu) ** 2) / sigma**2 / 2) / (sigma * np.sqrt(2 * np.pi))
        return gx

    plot_range = (df["delta_energy"].min(), df["delta_energy"].max())

    color_generator = rgbcolor([1] * len(df["frag_string"].unique()))
    generated_colors = [
        color_generator.hexcolor(i + 1, 1)
        for i in range(len(df["frag_string"].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [0.7, 0.3]}, sharex=True
    )

    x_range = np.linspace(plot_range[0], plot_range[1], num_point)
    y_sum = np.zeros(num_point)
    plot_data = dict()

    for plot_type in df["frag_string"].unique():
        y = np.zeros(num_point)
        temp_delta_energy = df[df["frag_string"] == plot_type].delta_energy
        for tde in temp_delta_energy:
            y += gaussian(x_range, sigma=sigma, mu=tde)
        y_sum += y
        plot_data[plot_type] = dict(curve_data=y, point_data=temp_delta_energy)

    counter = 0
    ratio = max(y_sum)
    for plot_type in plot_data:
        ax1.fill_between(
            x_range,
            0,
            plot_data[plot_type]["curve_data"] / ratio,
            alpha=alpha,
            color=generated_colors[counter],
            label=plot_type,
        )
        ax2.plot(
            plot_data[plot_type]["point_data"],
            np.zeros(len(plot_data[plot_type]["point_data"])) - counter,
            marker="|",
            lw=0,
            color=generated_colors[counter],
            alpha=alpha,
        )
        counter += 1
    # 归一化并绘制总和曲线

    ax1.plot(
        x_range, y_sum / ratio, color="black", label="Normalized Total Sum", linewidth=2
    )

    # 设置图表的整体属性
    ax1.set_ylabel("Intensity (arb. units)")
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False, fontsize=8)

    # ax2.set_ylabel('Count')
    ax2.set_yticks([])
    # ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel("KER (eV)")
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig


def plot_KER(df, col_name="state", sigma=1, alpha=0.3, num_point=1000):
    # FWHM
    # sigma = FWHM / (2*np.sqrt(2*np.log(2)))

    def gaussian(x, sigma=1, mu=0):
        gx = np.exp(-((x - mu) ** 2) / sigma**2 / 2) / (sigma * np.sqrt(2 * np.pi))
        return gx

    plot_range = (df["delta_energy"].min(), df["delta_energy"].max())

    color_generator = rgbcolor([1] * len(df[col_name].unique()))
    generated_colors = [
        color_generator.hexcolor(i + 1, 1) for i in range(len(df[col_name].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(6, 4), gridspec_kw={"height_ratios": [0.7, 0.3]}, sharex=True
    )

    x_range = np.linspace(plot_range[0], plot_range[1], num_point)
    y_sum = np.zeros(num_point)
    plot_data = dict()

    for plot_type in df[col_name].unique():
        y = np.zeros(num_point)
        temp_delta_energy = df[df[col_name] == plot_type].delta_energy
        for tde in temp_delta_energy:
            y += gaussian(x_range, sigma=sigma, mu=tde)
        y_sum += y
        plot_data[plot_type] = dict(curve_data=y, point_data=temp_delta_energy)

    counter = 0
    ratio = max(y_sum)
    for plot_type in plot_data:
        ax1.fill_between(
            x_range,
            0,
            plot_data[plot_type]["curve_data"] / ratio,
            alpha=alpha,
            color=generated_colors[counter],
            label=plot_type,
        )
        ax2.plot(
            plot_data[plot_type]["point_data"],
            np.zeros(len(plot_data[plot_type]["point_data"])) - counter,
            marker="|",
            ms=6,
            lw=0,
            color=generated_colors[counter],
            alpha=alpha,
        )
        counter += 1
    # 归一化并绘制总和曲线

    ax1.plot(
        x_range, y_sum / ratio, color="black", label="Normalized Total Sum", linewidth=2
    )

    # 设置图表的整体属性
    ax1.set_ylabel("Intensity (arb. units)")
    ax1.set_ylim(0.0, 1.05)
    ax1.legend(frameon=False, fontsize=8)

    # ax2.set_ylabel('Count')
    ax2.set_yticks([])
    ax2.set_ylim(-counter, 0.5)
    ax2.set_xlabel("KER (eV)")
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig


@st.fragment
def show_trajstats(trajs_data, trajs_stats):
    # plot_KER_state
    params_filter = st.columns(5)

    with params_filter[0]:
        min_steps = st.number_input(
            "轨迹的最小步数",
            value=1,
            min_value=1,
        )

    with params_filter[1]:
        min_nfrag = st.number_input(
            "片段数目最小值",
            value=0,
            min_value=0,
            max_value=trajs_stats["nfrag_max"],
        )

    with params_filter[2]:
        max_nfrag = st.number_input(
            "片段数目最大值",
            value=trajs_stats["nfrag_max"],
            min_value=trajs_stats["nfrag_min"],
            max_value=trajs_stats["nfrag_max"],
        )

    with params_filter[3]:
        min_delta_energy = st.number_input(
            "$\Delta E_{\min}$ (eV)",
            value=trajs_stats["delta_energy_min"],
            min_value=trajs_stats["delta_energy_min"] - 0.1,
            max_value=trajs_stats["delta_energy_max"] + 0.1,
        )

    with params_filter[4]:
        max_delta_energy = st.number_input(
            "$\Delta E_{\max}$ (eV)",
            value=trajs_stats["delta_energy_max"],
            min_value=trajs_stats["delta_energy_min"] - 0.1,
            max_value=trajs_stats["delta_energy_max"] + 0.1,
        )

    frag_string_selected = st.multiselect(
        label="片段类型",
        options=trajs_stats["all_frag_string"],
        default=trajs_stats["all_frag_string"],
    )

    states_selected = st.multiselect(
        label="轨迹的初始态",
        options=trajs_stats["all_states"],
        default=trajs_stats["all_states"],
    )

    trajs_data_filtered = trajs_data[
        (trajs_data.nstep >= min_steps)
        & (trajs_data.nfrag >= min_nfrag)
        & (trajs_data.nfrag <= max_nfrag)
        & (trajs_data.state.isin(states_selected))
        & (trajs_data.delta_energy >= min_delta_energy)
        & (trajs_data.delta_energy <= max_delta_energy)
        & (trajs_data.frag_string.isin(frag_string_selected))
    ]

    # 手动选择
    st.write("手动选择")

    # if "dff" not in st.session_state:
    # st.session_state.dff = trajs_data_filtered

    st.session_state.dff = trajs_data_filtered

    event_filtered = st.dataframe(
        st.session_state.dff,
        key="data_filtered",
        use_container_width=True,
        # height=1000,
        hide_index=True,
        column_order=("dir_path", "nstep", "nfrag", "frag_list", "delta_energy"),
        column_config={
            "dir_path": st.column_config.TextColumn(label="文件路径"),
            "nstep": st.column_config.NumberColumn(label="轨迹长度", width="small"),
            "nfrag": st.column_config.NumberColumn(label="片段数量", width="small"),
            "frag_list": st.column_config.ListColumn(label="片段成分", width="small"),
            "delta_energy": st.column_config.NumberColumn(
                label="初末能量差 (eV)", format="%.2f", width="small"
            ),
        },
        on_select="rerun",
        selection_mode=["multi-row"],
    )

    trajs_data_filtered_selected = trajs_data_filtered.iloc[
        event_filtered.selection["rows"], :
    ]
    # st.write(trajs_data_filtered_selected)

    # st.write(trajs_data_filtered_selected[["state","delta_energy"]])

    # 显示筛选后各个态的数目
    states_data = trajs_data_filtered_selected.state.value_counts()
    states_data = dict(state_type=list(states_data.index), state_num=list(states_data))
    # st.write(states_data)
    states_data_string_list = []
    for i in range(len(states_data["state_type"])):
        states_data_string_list.append(
            f"{states_data['state_num'][i]} `{states_data['state_type'][i]}`"
        )

    st.write(", ".join(states_data_string_list))
    # st.write(trajs_data_filtered)

    # 显示各个片段类型的数目
    frag_string_data = trajs_data_filtered_selected.frag_string.value_counts()
    frag_string_data = dict(
        frag_type=list(frag_string_data.index), frag_num=list(frag_string_data)
    )

    bar_data = frag_string_data
    for s in states_data["state_type"]:
        temp_df = trajs_data_filtered_selected[
            trajs_data_filtered_selected["state"] == s
        ]["frag_string"]
        temp_nstate = []
        for ft in frag_string_data["frag_type"]:
            temp_nstate.append(temp_df[temp_df == ft].count())
        bar_data[s] = temp_nstate

    bar_data = pd.DataFrame(bar_data)
    # bar_colors =

    # st.write(bar_data)
    # 显示各个片段类型的数目
    # , color_discrete_sequence=bar_colors
    fig = px.bar(
        bar_data,
        x="frag_type",
        y=sorted(states_data["state_type"]),
        color_discrete_sequence=generate_colors(len(states_data["state_type"])),
    )
    # fig = px.bar(frag_string_data, x='frag_type',y='frag_num',color='frag_num')  # 不对state分类
    fig.update_layout(
        title=dict(
            text=f"{len(trajs_data_filtered_selected)}条轨迹, {len(frag_string_data['frag_type'])}种解离路径"
        ),
        xaxis=dict(title="片段类型"),
        yaxis=dict(title="片段数目"),
    )
    st.plotly_chart(fig, use_container_width=True)

    # 片段中分子temp_frag
    all_frag_mol = []
    for temp_frag in frag_string_data["frag_type"]:
        all_frag_mol += [t.strip() for t in temp_frag.split("+")]
    all_frag_mol = sorted(set(all_frag_mol))

    # df_pie = {"frag_string": trajs_data_filtered['frag_string']}
    # for temp_frag in all_frag_mol:
    # df_pie[temp_frag] = trajs_data_filtered['frag_list'].apply(lambda x: temp_frag if temp_frag in x else 'No'+temp_frag)
    # st.write(df_pie)

    # 显示各个片段类型的比例
    df_pie = {
        "frag_type": frag_string_data["frag_type"],
        "frag_num": frag_string_data["frag_num"],
    }
    for temp_frag in all_frag_mol:
        # 分隔符' + '
        df_pie[temp_frag] = [
            (temp_frag if temp_frag in f.split(" + ") else "No " + temp_frag)
            for f in df_pie["frag_type"]
        ]
    # st.write(df_pie)

    df_pie = pd.DataFrame(df_pie)
    frag_mol_selected = st.multiselect(
        label="轨迹的片段所包含的分子",
        options=all_frag_mol,
        default=None,
        placeholder="先后顺序会影响结果",
        help="先后顺序会影响结果",
    )

    fig = px.sunburst(
        df_pie,
        path=frag_mol_selected + ["frag_type"],
        values="frag_num",
        color_discrete_sequence=generate_colors(len(frag_string_data["frag_type"])),
    )
    st.plotly_chart(fig, use_container_width=True)

    # state_type 类别大于1
    if len(trajs_data_filtered_selected["state"].unique()) > 1:
        params_hist = st.columns(2)
        with params_hist[0]:
            bin_count = st.number_input(
                "bin count",
                value=30,
                min_value=1,
                help="调整分箱数影响直方图精细度",
            )
        with params_hist[1]:
            kde_bandwidth = st.number_input(
                "kde bandwidth",
                value=0.2,
                min_value=0.1,
                help="带宽较大，曲线更平滑，可能合并邻近的峰；带宽较小，则曲线细节更丰富，可以显示更多的峰",
            )

        try:
            st.pyplot(
                plot_KER_state(
                    trajs_data_filtered_selected[["state", "delta_energy"]],
                    bin_count,
                    kde_bandwidth,
                ),
                use_container_width=True,
            )
            # # # 测试
            # st.pyplot(plot_KER( trajs_data_filtered_selected[["state", "delta_energy"]], col_name = "state", sigma=kde_bandwidth),
            # use_container_width=True)
        except Exception as e:
            st.warning(
                f"Failed to plot hist of `state`. The min of number of each kind should be larger than 1. {e}"
            )

    else:
        st.info("所选轨迹的态的种类不足")


# Dalitz Plot & Newton Plot
# ==========================

# def dnplot(
# xyzfiles: list[str] = glob("*let_*/TRAJ_*"),
# frag_orders: list[int] = [0, 1, 2],
# xylim: tuple[float] = (-1.0, 1.0, -1.0, 1.0),
# idx: int = 0,
# hdf5_file: str = "Data.hdf5",
# ) -> tuple:
# """

# TODO
# 检查片段的数目是否是3

# 计算片段的动能和（读取attrs 的frag 和atomic masses）
# """
# f = h5py.File(hdf5_file, "r")
# elements = [str(i) for i in range(3)]
# kins = []
# monts = []

# for path in (os.path.dirname(xf) for xf in xyzfiles):
# try:
# elements = f[path].attrs["elements"]
# atomic_masses = f[path].attrs["atomic_masses"]
# velocities = f[path]["velocity"][idx, :]

# # elements,atomic_masses,velocities = ptf.read_output_dat(f, idx)
# # print(elements,atomic_masses,velocities)

# p = np.array([m * v for m, v in zip(atomic_masses, velocities)])
# e = np.array(
# [0.5 * m * np.sum(v**2) for m, v in zip(atomic_masses, velocities)]
# )
# # print(p,e)
# monts.append(p)
# kins.append(e)
# except Exception as e:
# print(f"Failed to read file {hdf5_file} or calculate Ekin and p. {e}")

# kins = np.array(kins)
# monts = np.array(monts)

# # print(kins.shape)

# if len(kins) > 1:

# np.save("ps.npy", monts)
# np.save("es.npy", kins)

# elements = ["%s_{%d}" % (e, i + 1) for i, e in enumerate(elements)]

# return ptf.dalitzplot(kins, elements, orders, "top"), ptf.newtonplot(
# monts, elements, orders, xylim, "top"
# )
# else:
# print(f"No data to plot.")
# return None, None


@st.fragment
def show_dnplot(trajs_data, trajs_stats, analysis_params, trajs_dataset):
    # note 调用函数 dnplot
    params_dnplot = st.columns(3)

    with params_dnplot[0]:
        min_steps_dnplot = st.number_input(
            "轨迹的最小步数 ",
            value=1,
            min_value=1,
        )

    with params_dnplot[1]:
        min_delta_energy_dnplot = st.number_input(
            "$\Delta E_{\min}$ (eV) ",
            value=trajs_stats["delta_energy_min"],
            min_value=trajs_stats["delta_energy_min"] - 0.1,
            max_value=trajs_stats["delta_energy_max"] + 0.1,
        )

    with params_dnplot[2]:
        max_delta_energy_dnplot = st.number_input(
            "$\Delta E_{\max}$ (eV) ",
            value=trajs_stats["delta_energy_max"],
            min_value=trajs_stats["delta_energy_min"] - 0.1,
            max_value=trajs_stats["delta_energy_max"] + 0.1,
        )

    all_frag_string_3frag = [
        f for f in trajs_stats["all_frag_string"] if len(f.split(" + ")) == 3
    ]
    frag_string_selected_dnplot = st.multiselect(
        label="片段类型(只列出三个片段的解离) ",
        options=all_frag_string_3frag,
        default=all_frag_string_3frag,
    )

    states_selected_dnplot = st.multiselect(
        label="轨迹的初始态 ",
        options=trajs_stats["all_states"],
        default=trajs_stats["all_states"],
    )

    trajs_data_dnplot = trajs_data[
        (trajs_data.nstep >= min_steps_dnplot)
        & (trajs_data.nfrag == 3)
        & (trajs_data.state.isin(states_selected_dnplot))
        & (trajs_data.delta_energy >= min_delta_energy_dnplot)
        & (trajs_data.delta_energy <= max_delta_energy_dnplot)
        & (trajs_data.frag_string.isin(frag_string_selected_dnplot))
    ]

    # 手动选择
    st.write("选择筛选后的数据")

    st.session_state.df_dnplot = trajs_data_dnplot

    event_dnplot = st.dataframe(
        st.session_state.df_dnplot,
        key="data_dnplot",
        use_container_width=True,
        # height=1000,
        hide_index=True,
        column_order=("dir_path", "nstep", "nfrag", "frag_list", "delta_energy"),
        column_config={
            "dir_path": st.column_config.TextColumn(label="文件路径"),
            "nstep": st.column_config.NumberColumn(label="轨迹长度", width="small"),
            "nfrag": st.column_config.NumberColumn(label="片段数量", width="small"),
            "frag_list": st.column_config.ListColumn(label="片段成分", width="small"),
            "delta_energy": st.column_config.NumberColumn(
                label="初末能量差 (eV)", format="%.2f", width="small"
            ),
        },
        on_select="rerun",
        selection_mode=["multi-row"],
    )

    trajs_data_dnplot_selected = trajs_data_dnplot.iloc[
        event_dnplot.selection["rows"], :
    ]

    # st.write(trajs_data_dnplot_selected)
    st.write(f"已选中{trajs_data_dnplot_selected['dir_path'].count()}条轨迹")

    if len(trajs_data_dnplot_selected) > 2:
        frag_info_first = eval(trajs_data_dnplot_selected["frag"].iloc[0])
        # st.write([f"{f['formula']} | {f['index']}" for f in  frag_info_first.values()])
        frag_orders_options = [
            f"{f['formula']} | {sorted(f['index'])}" for f in frag_info_first.values()
        ]
        frag_orders = st.multiselect(
            label="片段顺序",
            options=frag_orders_options,
            default=frag_orders_options,
        )

        # st.write(frag_orders)

        st.write("Newton Plot 坐标轴范围")
        xylim_nplot = st.columns(4)
        with xylim_nplot[0]:
            xmin = st.number_input("$x_{\min}$", value=-1.0)
        with xylim_nplot[1]:
            xmax = st.number_input("$x_{\max}$", value=1.0)
        with xylim_nplot[2]:
            ymin = st.number_input("$y_{\min}$", value=-1.0)
        with xylim_nplot[3]:
            ymax = st.number_input("$y_{\max}$", value=1.0)

        idx = st.number_input(
            "提取速度所用帧的序号",
            value=-1 if READONLY_MODE else analysis_params["max_nstep"],
        )

        if len(frag_orders) == 3 and len(trajs_data_dnplot_selected) > 2:
            st.write(len(frag_orders), len(trajs_data_dnplot_selected))

            try:
                ekin_frag, momentum_frag = trajs_dataset.frag_ek_p(
                    trajs_data_dnplot_selected["dir_path"], idx, frag_orders
                )
                st.pyplot(
                    ptf.dalitzplot(
                        ekin_frag,
                        frag_name=[
                            foo.split("|")[0].strip() for foo in frag_orders_options
                        ],
                        location="top",
                    ),
                    use_container_width=True,
                )
                st.pyplot(
                    ptf.newtonplot(
                        momentum_frag,
                        frag_name=[
                            foo.split("|")[0].strip() for foo in frag_orders_options
                        ],
                        xylim=(xmin, xmax, ymin, ymax),
                        location="top",
                    ),
                    use_container_width=True,
                )
            except Exception as e:
                st.warning(f"Failed to plot. {e}")
    else:
        st.info("选择轨迹后画图")


# idx = -1
# frag_orders = ["H","H","O"]
# paths = ["State_1/TRAJ_00021","State_1/TRAJ_00042", "State_1/TRAJ_00045", "State_1/TRAJ_00094", "State_1/TRAJ_00096", "State_1/TRAJ_00100"]

# ekin_frag, momentum_frag = [], []

# for path in paths:
# try:
# f = h5py.File(HDF5_FILE, "r")
# frag_info, mass, veloc = eval(f[path].attrs["frag"]), f[path].attrs["atomic_masses"], f[path]["velocity"][idx]
# frag_info = dict((fv["formula"],fv["index"]) for fv in frag_info.values())
# fragidx_list = [frag_info.get(formula) for formula in frag_orders]
# momentum_atom = np.einsum("ij,i->ij", veloc, mass)
# momentum_frag.append([
# momentum_atom[fragidx].sum(axis=0)
# for fragidx in fragidx_list
# ])
# ekin_frag.append([
# np.sum((momentum_atom[fragidx].sum(axis=0)) ** 2)
# / (2 * mass[fragidx].sum())
# for fragidx in fragidx_list
# ])
# f.close()
# except Exception as e:
# print(f"{path} {e}")

# st.write(np.array(ekin_frag).shape)
# st.write(ptf.dalitzplot(np.array(ekin_frag)))
# st.write(np.array(momentum_frag))


# 查看同路径下的其它文件
@st.fragment
def show_file(dir_path):
    filename_display = st.text_input(
        "文件名:", value="input", help="文本文件一般包含 `geom`, `input`, `veloc`"
    )
    # path_temp = os.path.dirname(xyzfile) + f"/{filename_display}"
    path_temp = os.path.join(
        dir_path,
    )
    # 文件名非空且存在
    if filename_display and os.path.exists(path_temp):
        suffix = filename_display.split(".")[-1]
        # 打开图片文件
        if suffix in ["jpg", "png"]:
            st.image(path_temp)
        # 打开文本文件
        else:
            with open(path_temp, "r") as f:
                content = f.read()
            f.close()

            st.text(content)
    else:
        st.warning(f"文件{path_temp}不存在!")


# ================================================================
# 获取最大步数
def get_max_nstep(filename="KEYSTROKES.setup_traj"):
    with open(filename, "r") as f:
        lines = f.readlines()
    f.close()

    sim_time = 0
    sim_step = 1

    for line in lines:
        if "Simulation time (fs)" in line:
            temp = line.split()[0]
            sim_time = float(temp) if temp[0] != "#" else 1000
        elif "Simulation timestep (fs)" in line:
            temp = line.split()[0]
            sim_step = float(temp) if temp[0] != "#" else 0.5
        else:
            continue

    max_nstep = int(sim_time / sim_step)
    return max_nstep


def test_get_max_nstep():
    print(get_max_nstep())


# 主程序
# =======


def molfragapp():

    logger.info("\n%s", ASCII_LOGO)

    display_logo()
    st.title(f"分子解离片段分析 v{__VERSION}")
    # -------------------------------------------

    analysis_params = get_analysis_params()
    # st.write(analysis_params)

    trajs_dataset = sharc.Dataset(HDF5_FILE)

    ## 检查、更新文件
    ## ------------------------------------------
    if os.path.exists(HDF5_FILE):
        # 文件存在，则对比分析参数
        #   - 相同 -> 更新文件
        #   - 不同 -> 删除文件，从头创建
        st.toast(f"文件`{HDF5_FILE}`加载中...")
        if READONLY_MODE:
            st.toast("只读模式下不会更新数据.")
        # 非只读模式下，先检查分析参数，如果分析参数相同且启动更新数据选项，则对比修改时间并更新数据
        elif trajs_dataset.analysis_params() == analysis_params:
            # 默认为不更新数据，此时查看数据时页面执行速度更快
            update_data_on = st.toggle("更新数据(可能会消耗较长时间)")
            if update_data_on:
                st.toast("分析参数不变，开始更新文件...")

                mtime_old = trajs_dataset.mtime()
                mtime_new = dict(
                    (p, int(os.path.getmtime(p)))
                    for p in sorted(glob(analysis_params["dat_pattern"]))
                )
                # st.write('旧的修改时间',mtime_old,'新的修改时间',mtime_new)
                file_path_to_update = []
                for dir_path, mtime in mtime_new.items():
                    if mtime_old.get(dir_path) != mtime:
                        # print(dir_path,'需要更新')
                        # print(mtime_old.get(dir_path) , mtime)
                        file_path_to_update.append(dir_path)
                # st.write('需要更新数据的xyz文件',file_path_to_update)
                trajs_dataset.update(
                    sharc.parallel_trajdata(
                        sorted(file_path_to_update), analysis_params["max_bond_length"]
                    ),
                    analysis_params,
                )

                st.toast("文件更新结束")
        else:
            # 分析参数不同，则备份原始文件，按照新的参数提取数据
            st.toast("分析参数改变，重新提取数据...")
            os.rename(HDF5_FILE, HDF5_FILE + ".bk")
            # 测试数据
            trajs_dataset.create(
                sharc.parallel_trajdata(
                    sorted(glob(analysis_params["dat_pattern"])),
                    analysis_params["max_bond_length"],
                ),
                analysis_params,
            )
    else:
        # 文件不存在，则按照分析参数提取数据并写入数据文件
        st.toast(f"未发现文件`{HDF5_FILE}`, 正在提取数据,请耐心等待")

        # with st.spinner("Wait for it...", show_time=True):  # streamlit >= 1.42
        with st.spinner("正在提取数据,请耐心等待"):
            trajs_dataset.create(
                sharc.parallel_trajdata(
                    sorted(glob(analysis_params["dat_pattern"])),
                    analysis_params["max_bond_length"],
                ),
                analysis_params,
            )
        st.toast("数据提取结束")

    ## 加载、统计数据
    ## ------------------------------------------
    trajs_data = pd.DataFrame(trajs_dataset.fragdata())
    # st.write(trajs_data)

    cols_button = st.columns(2)
    with cols_button[0]:
        if st.button("导出所有数据"):
            csv_export = "data_export.csv"
            trajs_data.to_csv(csv_export)
            st.toast(f"数据已导出到文件`{csv_export}`")
    with cols_button[1]:
        counts_stats(trajs_data.nstep, "nstep", 20)

    # 统计数据
    trajs_stats = {}
    trajs_stats["nfrag_max"], trajs_stats["nfrag_min"] = (
        trajs_data.nfrag.max(),
        trajs_data.nfrag.min(),
    )
    trajs_stats["all_states"], trajs_stats["all_frag_string"] = sorted(
        set(trajs_data.state)
    ), sorted(set(trajs_data.frag_string))
    trajs_stats["delta_energy_max"], trajs_stats["delta_energy_min"] = (
        trajs_data.delta_energy.max(),
        trajs_data.delta_energy.min(),
    )

    # st.write(trajs_stats)

    ## 分析、展示数据
    ## ------------------------------------------
    st.header("解离片段分析结果")

    cols = st.columns(2)
    # 左栏
    with cols[0]:
        st.session_state.df = trajs_data

        event = st.dataframe(
            st.session_state.df,
            key="data",
            use_container_width=True,
            height=1000,
            hide_index=True,
            # 设置需要显示的列
            column_order=("dir_path", "nstep", "nfrag", "frag_list", "delta_energy"),
            column_config={
                "dir_path": st.column_config.TextColumn(label="文件路径"),
                "nstep": st.column_config.NumberColumn(label="轨迹长度", width="small"),
                "nfrag": st.column_config.NumberColumn(label="片段数量", width="small"),
                "frag_list": st.column_config.ListColumn(label="片段成分"),
                "delta_energy": st.column_config.NumberColumn(
                    label="初末动能差 (eV)", format="%.2f", width="small"
                ),
            },
            on_select="rerun",
            selection_mode=["single-row"],
        )

    # 右栏
    with cols[1]:
        # st.session_state
        tabs = st.tabs(
            [
                "轨迹动画&能量曲线",
                "轨迹筛选与统计",
                "模拟Dalitz图和Newton图",
                "查看其它文件",
            ]
        )
        # 第一个标签
        with tabs[0]:
            # 查看轨迹
            # 需要使用 event.selection
            if len(event.selection["rows"]) == 1:
                idx_row = event.selection["rows"][0]
                # st.write(st.session_state.df.iloc[idx_row])
                # st.write(trajs_dataset.expec_out(st.session_state.df.iloc[idx_row]["dir_path"]))
                # show_traj(
                # # xyzfile=st.session_state.df.iloc[idx_row]["xyz_path"],
                # expec_out_df = trajs_dataset.expec_out(st.session_state.df.iloc[idx_row]["dir_path"]) ,
                # frag_info=st.session_state.df.iloc[idx_row]["frag"],
                # )
                dir_path = st.session_state.df.iloc[idx_row]["dir_path"]
                st.write(dir_path)
                try:
                    show_traj(*trajs_dataset.traj4show(dir_path))
                except Exception as e:
                    st.warning(e)
            else:
                st.info("选择一个轨迹查看动画与能量")

        # 第二个标签
        with tabs[1]:
            # 需要 trajs_stats
            try:
                show_trajstats(trajs_data=trajs_data, trajs_stats=trajs_stats)
            except Exception as e:
                st.warning(e)

        # 第三个标签
        with tabs[2]:
            # 需要trajs_stats trajs_data
            try:
                # show_dnplot(
                # trajs_data=trajs_data,
                # trajs_stats=trajs_stats,
                # analysis_params=analysis_params,
                # hdf5_file=HDF5_FILE,
                # )

                show_dnplot(
                    trajs_data=trajs_data,
                    trajs_stats=trajs_stats,
                    analysis_params=analysis_params,
                    trajs_dataset=trajs_dataset,
                )

                # st.write( trajs_dataset.frag_ek_p(["State_1/TRAJ_00001","State_1/TRAJ_00002"],-1,["H","H","O"]) )

            except Exception as e:
                st.warning(e)

        # 第四个标签
        with tabs[3]:
            if READONLY_MODE:
                st.info("只读模式下不可用")
            else:
                # 检查输入
                if len(event.selection["rows"]) == 1:
                    idx_row = event.selection["rows"][0]
                    try:
                        show_file(
                            dir_path=st.session_state.df.iloc[idx_row]["dir_path"]
                        )
                    except Exception as e:
                        st.warning(e)
                else:
                    st.info("选择一个轨迹查看同一路径下的其它文件")

    display_copyright()
    # -------------------------------------------


if __name__ == "__main__":
    molfragapp()
