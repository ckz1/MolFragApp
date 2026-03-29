r"""MolFragApp
===============

分子动力学数据处理

todo
-----

changelog
----------
2024.12.27  Feat: Location of colorbar

"""

import datetime
import multiprocessing as mp
import os
import time
from glob import glob

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
import psutil
import streamlit as st
import streamlit.components.v1 as components
import xyztraj
from plotly.subplots import make_subplots
from utils import *

VERSION = '2024.12.27'
st.set_page_config(page_icon="🔍", page_title='MolFragApp', layout='wide')


def display_logo():
    """显示logo
    """
    try:
        file_path = os.path.realpath(__file__)
        logo_path = os.path.join(os.path.dirname(file_path),
                                 "MolFragApp_logo.jpg")
        st.logo(logo_path,
                size='large',
                link='https://github.com/ckz1/MolFragApp')
    except:
        print("Failed to display the logo")


class TrajDir:
    """使用SHARC进行分子动力学模拟时一条轨迹的目录

    Attributes:
        path (str): `xyz` 文件路径
        dirpath (str): 目录路径
        traj (object): 轨迹对象
        nstep (int): 轨迹步数
        frag (dict): 片段
        frag_list (list[str]): 片段化学式列表
        state (str): 态名称，即多重度与序号
        multiplicity (str): 多重度
        expec_out (str): `output_data/expec.out` 文件路径
        delta_energy (float): 初末能量差
    """

    def __init__(self, path='PATH/of/output.xyz'):
        # 路径必须存在, 使用相对路径
        # self.path = os.path.abspath(path)
        self.path = path
        self.dirpath = os.path.dirname(self.path)

        self.traj = None
        self.nstep = 0
        self.frag = None
        self.frag_list = []
        self.nfrag = 0
        self.state = ''
        self.multiplicity = ''
        self.expec_out = f"{self.dirpath}/output_data/expec.out"
        self.delta_energy = np.nan

    def update_data(self, max_bond_length: float = 2.5) -> None:
        """提取数据

        Args:
            max_bond_length (float, optional): 划分片段的最大键长. Defaults to 2.5.
        """
        try:
            # 轨迹信息
            self.traj = xyztraj.Traj(self.path)
            self.nstep = len(self.traj)
            self.frag, self.frag_list = self.traj[-1].get_frag(max_bond_length)
            self.frag = str(self.frag)
            self.nfrag = len(self.frag_list)
            # 激发态与多重度信息
            self.state = self.path.split('/')[-3]
            self.multiplicity = self.state.split('_')[0]
        except:
            print(f'Failed to update data for `{self.path}`')

    def get_delta_energy(self) -> None:
        r"""读取 `output.lis` 计算初末能量差 :math:`\Delta E = E_0 - \bar{E}_{last\ 100}`
        """
        try:
            output_lis_file = f"{self.dirpath}/output.lis"
            output_lis = np.loadtxt(output_lis_file)
            # self.delta_energy = output_lis[0,5] - output_lis[-1,5]
            # total energy of first frame - mean of potential energy of last 100 frames
            self.delta_energy = output_lis[0, 6] - np.mean(
                output_lis[-100:, 5])  # eV
        except:
            print(f'Failed to get delta_energy from `{self.dirpath}`')

    def generate_expec_out(self) -> None:
        """生成文件 `output_data/expec.out`
        """
        try:
            res = os.system(
                f"cd {self.dirpath}; $SHARC/data_extractor.x output.dat > /dev/null 2>&1 &"
            )
        except:
            print(f"Failed to generate expec.out in {self.dirpath}")

    def __repr__(self):
        return f"Path {self.path}\nSteps {self.nstep}\nFrag {self.frag_list}\nState {self.state} {self.multiplicity}\nDeltaE {self.delta_energy}"


def get_single_traj_info(params: tuple[str, float, int]) -> list:
    """获取后续分析所需的轨迹信息

    Args:
        params (tuple[str,float,int]): 文件名, 最大键长, 最大步数

    Returns:
        list: 获取的轨迹信息
    """
    xyzfile, max_bond_length, max_nstep = params
    temp_traj = TrajDir(xyzfile)
    temp_traj.update_data(max_bond_length)
    temp_traj.get_delta_energy()
    if temp_traj.nstep >= max_nstep + 1 > 1 and os.path.exists(
            temp_traj.expec_out):
        print(f"Trajectory {temp_traj.path} finished.")
    else:
        temp_traj.generate_expec_out()
    # 此处输出的顺序应该与`frag_analysis` 函数的`columns`一致
    return [
        temp_traj.path, temp_traj.state, temp_traj.nstep, temp_traj.nfrag,
        temp_traj.frag_list, ' + '.join(temp_traj.frag_list), temp_traj.frag,
        temp_traj.delta_energy, temp_traj.expec_out
    ]


# get_all_traj_info
@st.cache_data
def frag_analysis(xyz_path: str,
                  max_bond_length: float = 2.5,
                  max_nstep: int = 0) -> pd.DataFrame:
    """分析 `xyz_path` 匹配的所有轨迹

    Args:
        xyz_path (str): `xyz` 文件路径的模式
        max_bond_length (float, optional): 最大键长. Defaults to 2.5.
        max_nstep (int, optional): 最大步数. Defaults to 0.

    Returns:
        pd.DataFrame: 所有轨迹信息 DataFrame
    """

    xyzfiles = sorted(glob(xyz_path))

    ncpu = psutil.cpu_percent(interval=1, percpu=True).count(0)
    ncpu = max(1, ncpu)

    with mp.Pool(ncpu) as p:
        # results = p.map(get_single_traj_info,((os.path.abspath(f),max_bond_length, max_nstep) for f in xyzfiles))
        results = p.map(get_single_traj_info,
                        ((f, max_bond_length, max_nstep) for f in xyzfiles))

    return pd.DataFrame(results,
                        columns=[
                            "file", "state", "nstep", "nfrag", "frag_list",
                            "frag_string", "frag", "delta_energy", "expec_out"
                        ])


# 布局
#==============
def main():
    display_logo()

    st.title(f'分子解离片段分析 v{VERSION}')
    st.write(f'当前路径:`{os.getcwd()}`')

    #path1, path2, para3, para4 = st.columns(3)
    path1, path2, para3, para4 = st.columns(4)

    with path1:
        xyz_path = st.text_input(
            "轨迹文件(`xyz`格式)路径:",
            value="*let_*/TRAJ_*/output.xyz",
            help='使用`xyz`文件最后一帧的结构分析解离片段的成分',
        )
        num_xyzfile = len(glob(xyz_path))
        st.write("轨迹文件总数：", num_xyzfile)

    if num_xyzfile == 0:
        st.warning('No xyz file found!', icon="⚠️")

    with path2:
        template_path = st.text_input(
            "`template`文件路径:",
            value="*.template",
        )
        # st.write(read_template(glob(template_path)[0]))
        try:
            cal_method = read_template(glob(template_path)[0])
        except:
            cal_method = f"Warning: Can NOT read template from {template_path}!"

        st.write(f"计算级别:`{cal_method}`")

    with para3:
        setup_traj_path = st.text_input(
            "`KEYSTROKES.setup_traj`文件路径:",
            value="KEYSTROKES.setup_traj",
        )
        try:
            max_nstep = get_max_nstep(setup_traj_path)
        except:
            max_nstep = 0

        st.write(f"轨迹步数设置：`{max_nstep}`")

    with para4:
        max_bond_length = st.number_input(
            "划分片段时的最大键长(Å)",
            value=2.5,
            min_value=0.1,
            help='若两原子之间距离大于最大键长，则被划分到不同片段',
        )

    # 结果展示
    st.header('解离片段分析结果')
    start_time = time.time()
    trajs_data = frag_analysis(xyz_path, max_bond_length, max_nstep)
    stop_time = time.time()
    st.write(f"耗时{stop_time - start_time:.2f} s")
    # 保存数据
    csv_export = 'data_export.csv'
    trajs_data.to_csv(csv_export)
    # st.toast(f"片段分析结果已导出为`{csv_export}`")

    nfrag_max, nfrag_min = trajs_data.nfrag.max(), trajs_data.nfrag.min()
    all_states = sorted(set(trajs_data.state))
    all_frag_string = sorted(set(trajs_data.frag_string))
    delta_energy_max, delta_energy_min = trajs_data.delta_energy.max(
    ), trajs_data.delta_energy.min()

    # DEFAULT_COLORS = [
    # '#FF0000', '#FF7F00', '#FFFF00', '#7FFF00', '#00FF00', '#00FF7F',
    # '#00FFFF', '#0000FF', '#FF00FF'
    # ]
    # # ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
    # DEFAULT_COLORS_LIGHT = [
    # "#FFB2B2",
    # "#D8FFB2",
    # "#B2FFD0",
    # "#B2E7FF",
    # "#C1B2FF",
    # "#FFB2F7",

    # # "#FF8080",  # 浅红色
    # # "#FFB84D",  # 浅橙色
    # # "#FFFF80",  # 浅黄色
    # # "#B8FF66",  # 浅黄绿色
    # # "#80FF80",  # 浅绿色
    # # "#80FFB3",  # 浅青绿色
    # # "#80FFFF",  # 浅蓝绿色
    # # "#8080FF",  # 浅蓝色
    # # "#FF80FF"  # 浅品红色
    # ]

    col_data, col_analysis = st.columns(2)

    with col_data:

        st.session_state.df = trajs_data

        event = st.dataframe(
            st.session_state.df,
            key="data",
            use_container_width=True,
            height=1000,
            hide_index=True,
            column_order=("file", "nstep", "nfrag", "frag_list",
                          "delta_energy"),
            column_config={
                "file":
                st.column_config.TextColumn(label='文件路径'),
                "nstep":
                st.column_config.NumberColumn(label='轨迹长度', width='small'),
                "nfrag":
                st.column_config.NumberColumn(label='片段数量', width='small'),
                "frag_list":
                st.column_config.ListColumn(label='片段成分'),
                "delta_energy":
                st.column_config.NumberColumn(label='初末能量差 (eV)',
                                              format="%.2f",
                                              width='small'),
            },
            on_select="rerun",
            selection_mode=["single-row"],
        )

        # event.selection

    with col_analysis:
        tab1, tab2, tab3, tab4 = st.tabs(
            ["轨迹筛选与统计", "轨迹动画&能量曲线", "模拟Dalitz图和Newton图", "查看其它文件"])

        with tab1:
            para1_filter, para2_filter, para3_filter, para4_filter, para5_filter = st.columns(
                5)

            with para1_filter:
                min_steps = st.number_input(
                    "轨迹的最小步数",
                    value=1,
                    min_value=1,
                )

            with para2_filter:
                min_nfrag = st.number_input(
                    "片段数目最小值",
                    value=0,
                    min_value=0,
                    max_value=nfrag_max,
                )

            with para3_filter:
                max_nfrag = st.number_input(
                    "片段数目最大值",
                    value=nfrag_max,
                    min_value=nfrag_min,
                    max_value=nfrag_max,
                )

            with para4_filter:
                min_delta_energy = st.number_input(
                    "$\Delta E_{\min}$ (eV)",
                    value=delta_energy_min,
                    min_value=delta_energy_min - .1,
                    max_value=delta_energy_max + .1,
                )

            with para5_filter:
                max_delta_energy = st.number_input(
                    "$\Delta E_{\max}$ (eV)",
                    value=delta_energy_max,
                    min_value=delta_energy_min - .1,
                    max_value=delta_energy_max + .1,
                )

            frag_string_selected = st.multiselect(
                label='片段类型',
                options=all_frag_string,
                default=all_frag_string,
            )

            states_selected = st.multiselect(
                label='轨迹的初始态',
                options=all_states,
                default=all_states,
            )

            trajs_data_filtered = trajs_data[
                (trajs_data.nstep >= min_steps)
                & (trajs_data.nfrag >= min_nfrag) &
                (trajs_data.nfrag <= max_nfrag) &
                (trajs_data.state.isin(states_selected)) &
                (trajs_data.delta_energy >= min_delta_energy) &
                (trajs_data.delta_energy <= max_delta_energy) &
                (trajs_data.frag_string.isin(frag_string_selected))]

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
                column_order=("file", "nstep", "nfrag", "frag_list",
                              "delta_energy"),
                column_config={
                    "file":
                    st.column_config.TextColumn(label='文件路径'),
                    "nstep":
                    st.column_config.NumberColumn(label='轨迹长度', width='small'),
                    "nfrag":
                    st.column_config.NumberColumn(label='片段数量', width='small'),
                    "frag_list":
                    st.column_config.ListColumn(label='片段成分', width='small'),
                    "delta_energy":
                    st.column_config.NumberColumn(label='初末能量差 (eV)',
                                                  format="%.2f",
                                                  width='small'),
                },
                on_select="rerun",
                selection_mode=["multi-row"],
            )

            trajs_data_filtered_selected = trajs_data_filtered.iloc[
                event_filtered.selection["rows"], :]
            # st.write(trajs_data_filtered_selected)

            # st.write(trajs_data_filtered_selected[["state","delta_energy"]])

            # 显示筛选后各个态的数目
            states_data = trajs_data_filtered_selected.state.value_counts()
            states_data = dict(state_type=list(states_data.index),
                               state_num=list(states_data))
            # st.write(states_data)
            states_data_string_list = []
            for i in range(len(states_data['state_type'])):
                states_data_string_list.append(
                    f"{states_data['state_num'][i]} `{states_data['state_type'][i]}`"
                )

            st.write(', '.join(states_data_string_list))
            # st.write(trajs_data_filtered)

            # 显示各个片段类型的数目
            frag_string_data = trajs_data_filtered_selected.frag_string.value_counts(
            )
            frag_string_data = dict(frag_type=list(frag_string_data.index),
                                    frag_num=list(frag_string_data))

            bar_data = frag_string_data
            for s in states_data['state_type']:
                temp_df = trajs_data_filtered_selected[
                    trajs_data_filtered_selected['state'] == s]['frag_string']
                temp_nstate = []
                for ft in frag_string_data['frag_type']:
                    temp_nstate.append(temp_df[temp_df == ft].count())
                bar_data[s] = temp_nstate

            bar_data = pd.DataFrame(bar_data)
            # bar_colors =

            # st.write(bar_data)
            # 显示各个片段类型的数目
            # , color_discrete_sequence=bar_colors
            fig = px.bar(bar_data,
                         x='frag_type',
                         y=sorted(states_data['state_type']),
                         color_discrete_sequence=generate_colors(
                             len(states_data['state_type'])))
            # fig = px.bar(frag_string_data, x='frag_type',y='frag_num',color='frag_num')  # 不对state分类
            fig.update_layout(title=dict(
                text=
                f"{len(trajs_data_filtered_selected)}条轨迹, {len(frag_string_data['frag_type'])}种解离路径"
            ),
                              xaxis=dict(title='片段类型'),
                              yaxis=dict(title='片段数目'))
            st.plotly_chart(fig, use_container_width=True)

            # 片段中分子temp_frag
            all_frag_mol = []
            for temp_frag in frag_string_data['frag_type']:
                all_frag_mol += [t.strip() for t in temp_frag.split('+')]
            all_frag_mol = sorted(set(all_frag_mol))

            # df_pie = {"frag_string": trajs_data_filtered['frag_string']}
            # for temp_frag in all_frag_mol:
            # df_pie[temp_frag] = trajs_data_filtered['frag_list'].apply(lambda x: temp_frag if temp_frag in x else 'No'+temp_frag)
            # st.write(df_pie)

            # 显示各个片段类型的比例
            df_pie = {
                "frag_type": frag_string_data['frag_type'],
                "frag_num": frag_string_data['frag_num'],
            }
            for temp_frag in all_frag_mol:
                # 分隔符' + '
                df_pie[temp_frag] = [
                    (temp_frag if temp_frag in f.split(' + ') else 'No ' +
                     temp_frag) for f in df_pie["frag_type"]
                ]
            # st.write(df_pie)

            df_pie = pd.DataFrame(df_pie)
            frag_mol_selected = st.multiselect(
                label='轨迹的片段所包含的分子',
                options=all_frag_mol,
                default=None,
                placeholder='先后顺序会影响结果',
                help='先后顺序会影响结果',
            )

            fig = px.sunburst(df_pie,
                              path=frag_mol_selected + ["frag_type"],
                              values="frag_num",
                              color_discrete_sequence=generate_colors(
                                  len(frag_string_data['frag_type'])))
            st.plotly_chart(fig, use_container_width=True)

            # state_type 类别大于1
            if len(trajs_data_filtered_selected["state"].unique()) > 1:
                param1_hist, param2_hist = st.columns(2)
                with param1_hist:
                    bin_count = st.number_input("bin count",
                                                value=30,
                                                min_value=1,
                                                help="调整分箱数影响直方图精细度")
                with param2_hist:
                    kde_bandwidth = st.number_input(
                        "kde bandwidth",
                        value=0.2,
                        min_value=0.1,
                        help="带宽较大，曲线更平滑，可能合并邻近的峰；带宽较小，则曲线细节更丰富，可以显示更多的峰")

                try:
                    st.pyplot(plot_KER_state(
                        trajs_data_filtered_selected[["state",
                                                      "delta_energy"]],
                        bin_count, kde_bandwidth),
                              use_container_width=True)
                    # # # 测试
                    # st.pyplot(plot_KER( trajs_data_filtered_selected[["state", "delta_energy"]], col_name = "state", sigma=kde_bandwidth),
                    # use_container_width=True)
                except:
                    st.warning(
                        "Failed to plot hist of `state`. The min of number of each kind should be larger than 1."
                    )

            else:
                st.info("所选轨迹的态的种类不足")

            # if len(trajs_data_filtered_selected["frag_string"].unique()) > 1:
            # param1_hist, param2_hist = st.columns(2)
            # with param1_hist:
            # bin_count = st.number_input("bin count ",
            # value=30,
            # min_value=1,
            # help="调整分箱数影响直方图精细度")
            # with param2_hist:
            # kde_bandwidth = st.number_input(
            # "kde bandwidth ",
            # value=0.2,
            # min_value=0.1,
            # help="带宽较大，曲线更平滑，可能合并邻近的峰；带宽较小，则曲线细节更丰富，可以显示更多的峰")
            # # st.pyplot(plot_KER_frag(trajs_data_filtered_selected[["frag_string", "delta_energy"]],
            # # bin_count, kde_bandwidth),
            # # use_container_width=True)
            # try:
            # # # 测试
            # st.pyplot(plot_KER(trajs_data_filtered_selected[["frag_string", "delta_energy"]], col_name="frag_string",
            # sigma=kde_bandwidth),
            # use_container_width=True)

            # st.pyplot(plot_KER_frag(trajs_data_filtered_selected[["frag_string", "delta_energy"]],
            # bin_count, kde_bandwidth),
            # use_container_width=True)
            # except:
            # st.warning("Failed to plot hist of `frag_string`. The min of number of each kind should be larger than 1.")
            # else:
            # st.info("所选轨迹的片段类型的种类不足")

            # # 显示能量差的直方图
            # hist_data = [
            # trajs_data_filtered_selected[trajs_data_filtered_selected['state']
            # == st]["delta_energy"]
            # for st in states_data['state_type']
            # ]

            # # state_type 类别大于1
            # if len(trajs_data_filtered_selected["state"].unique()) > 1:
            # bin_size = st.number_input(
            # "bin size",
            # value=0.25,
            # min_value=0.1,
            # )
            # fig = ff.create_distplot(hist_data,
            # states_data['state_type'],
            # bin_size=bin_size,
            # colors=DEFAULT_COLORS_LIGHT)
            # # , curve_type ='normal'
            # fig.update_layout(xaxis_title='ΔE (eV)')
            # st.plotly_chart(fig, use_container_width=True)
            # else:
            # st.info("所选轨迹的态的种类不足")

        with tab2:
            if len(event.selection["rows"]) == 1:
                idx_row = event.selection["rows"][0]
                xyzfile = st.session_state.df.iloc[idx_row]['file']
                expec_out_file = st.session_state.df.iloc[idx_row]['expec_out']
                incoord_data = None

                fraginfo_on = st.toggle("显示片段信息")
                if fraginfo_on:
                    st.write(st.session_state.df.iloc[idx_row]['frag'])

                    # TODO multiselect
                    st.write("选择原子序号计算内坐标(原子序号从1开始)")
                    idx1_selected, idx2_selected, idx3_selected, idx4_selected = st.columns(
                        4)
                    with idx1_selected:
                        atom_idx1 = st.number_input("选择第一个原子",
                                                    value=0,
                                                    min_value=0)
                    with idx2_selected:
                        atom_idx2 = st.number_input("选择第二个原子",
                                                    value=0,
                                                    min_value=0)
                    with idx3_selected:
                        atom_idx3 = st.number_input("选择第三个原子",
                                                    value=0,
                                                    min_value=0)
                    with idx4_selected:
                        atom_idx4 = st.number_input("选择第四个原子",
                                                    value=0,
                                                    min_value=0)

                    atom_idxs = [atom_idx1, atom_idx2, atom_idx3, atom_idx4]
                    atom_idxs = [i for i in atom_idxs if i != 0]

                    if len(atom_idxs) > 0:
                        incoord_data = cal_incoord(xyzfile, atom_idxs)
                    # st.write(atom_idxs)
                    # st.write(incoord_data)

                show_mol(xyzfile, width=850, height=350)
                show_energy_curve(expec_out_file, incoord_data)
            else:
                st.info("选择一个轨迹查看动画与能量")

        with tab3:
            pass
            para1_dnplot, para2_dnplot, para3_dnplot = st.columns(3)

            with para1_dnplot:
                min_steps_dnplot = st.number_input(
                    "轨迹的最小步数 ",
                    value=1,
                    min_value=1,
                )

            with para2_dnplot:
                min_delta_energy_dnplot = st.number_input(
                    "$\Delta E_{\min}$ (eV) ",
                    value=delta_energy_min,
                    min_value=delta_energy_min - .1,
                    max_value=delta_energy_max + .1,
                )

            with para3_dnplot:
                max_delta_energy_dnplot = st.number_input(
                    "$\Delta E_{\max}$ (eV) ",
                    value=delta_energy_max,
                    min_value=delta_energy_min - .1,
                    max_value=delta_energy_max + .1,
                )

            all_frag_string_3frag = [
                f for f in all_frag_string if len(f.split(' + ')) == 3
            ]
            frag_string_selected_dnplot = st.multiselect(
                label='片段类型 ',
                options=all_frag_string_3frag,
                default=all_frag_string_3frag,
            )

            states_selected_dnplot = st.multiselect(
                label='轨迹的初始态 ',
                options=all_states,
                default=all_states,
            )

            trajs_data_dnplot = trajs_data[
                (trajs_data.nstep >= min_steps_dnplot)
                & (trajs_data.nfrag == 3) &
                (trajs_data.state.isin(states_selected_dnplot)) &
                (trajs_data.delta_energy >= min_delta_energy_dnplot) &
                (trajs_data.delta_energy <= max_delta_energy_dnplot) &
                (trajs_data.frag_string.isin(frag_string_selected_dnplot))]

            # 手动选择
            st.write("手动选择")

            st.session_state.df_dnplot = trajs_data_dnplot

            event_dnplot = st.dataframe(
                st.session_state.df_dnplot,
                key="data_dnplot",
                use_container_width=True,
                # height=1000,
                hide_index=True,
                column_order=("file", "nstep", "nfrag", "frag_list",
                              "delta_energy"),
                column_config={
                    "file":
                    st.column_config.TextColumn(label='文件路径'),
                    "nstep":
                    st.column_config.NumberColumn(label='轨迹长度', width='small'),
                    "nfrag":
                    st.column_config.NumberColumn(label='片段数量', width='small'),
                    "frag_list":
                    st.column_config.ListColumn(label='片段成分', width='small'),
                    "delta_energy":
                    st.column_config.NumberColumn(label='初末能量差 (eV)',
                                                  format="%.2f",
                                                  width='small'),
                },
                on_select="rerun",
                selection_mode=["multi-row"],
            )

            trajs_data_dnplot_selected = trajs_data_dnplot.iloc[
                event_dnplot.selection["rows"], :]

            # st.write(trajs_data_dnplot_selected)
            st.write(f"{trajs_data_dnplot_selected['file'].count()}条轨迹")

            atom_orders = st.multiselect(
                label='原子顺序',
                options=[1, 2, 3],
                default=[1, 2, 3],
            )

            # st.write(atom_orders)

            st.write("Newton Plot 坐标轴范围")
            xylim_nplot1, xylim_nplot2, xylim_nplot3, xylim_nplot4 = st.columns(
                4)
            with xylim_nplot1:
                xmin = st.number_input("$x_{\min}$", value=-1.)
            with xylim_nplot2:
                xmax = st.number_input("$x_{\max}$", value=1.)
            with xylim_nplot3:
                ymin = st.number_input("$y_{\min}$", value=-1.)
            with xylim_nplot4:
                ymax = st.number_input("$y_{\max}$", value=1.)

            idx = st.number_input("提取速度所用帧的序号", value=max_nstep)

            if len(atom_orders) == 3:
                atom_orders = [ao - 1 for ao in atom_orders]

                # try:
                dalitzplot, newtonplot = dnplot(
                        trajs_data_dnplot_selected['file'], atom_orders,
                        (xmin, xmax, ymin, ymax), idx)

                # dalitzplot, newtonplot = dnplot(trajs_data_dnplot_selected['file'],atom_orders)
                st.pyplot(dalitzplot, use_container_width=True)
                st.pyplot(newtonplot, use_container_width=True)
                # except:
                    # st.warning("Failed to plot.")

        with tab4:
            if len(event.selection["rows"]) == 1:
                idx_row = event.selection["rows"][0]
                xyzfile = st.session_state.df.iloc[idx_row]['file']
                show_file(xyzfile)
            else:
                st.info("选择一个轨迹查看同一路径下的其它文件")

    st.divider()
    components.html('''
    <p align='center'> © Copyright
    <a href="https://github.com/ckz1/MolFragApp" target="_blank">
    <img border="0" src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png" alt="MolFragApp in Github" width="30">
    </a>
    </p>
    ''')


if __name__ == '__main__':
    main()