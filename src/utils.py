import os
import psutil
import json
import numpy as np
import multiprocessing as mp 
import py3Dmol
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import colorsys
import matplotlib as mpl
import xyztraj
import py3bodyfrag as ptf
import logger
import h5py
# =======================================================================



# test_get_single_traj_info()



# =======================================================================

def read_template(filename='MOLCAS.template', output_type='str'):
    assert os.path.exists(filename)

    if 'MOLCAS.template' in filename:
        qc_software = 'OpenMolcas'
    elif 'MOLPRO.template' in filename:
        qc_software = 'Molpro'

    with open(filename, 'r') as f:
        lines = f.readlines()
    f.close()

    template = {'qc_software': qc_software, 'method' : 'CASSCF(default)'}
    for line in lines:
        temp = line.split()
        if len(temp) >= 2 and '#' not in line:
            template[temp[0]] = temp[1]

    if output_type == 'str':
        if qc_software == 'OpenMolcas':
            return f"SA({template['roots']})-{template['method']}({template['nactel']}e,{template['ras2']}o)/{template['basis']}, {template['qc_software']}"
        elif qc_software == 'Molpro':
            return f"SA({template['roots']})-{template['method']}({int(template['nelec'])-int(template['closed'])*2}e,{int(template['occ'])-int(template['closed'])}o)/{template['basis']}, {template['qc_software']}"

    return template
    
def read_output_dat(
        output_data_file: str = 'output.dat') -> tuple[dict, dict, list[dict]]:
    """读取 `output.dat` 中的所有信息

    Args:
        output_data_file (str, optional): `output.dat` 文件路径. Defaults to 'output.dat'.

    Returns:
        tuple[dict,dict,list[dict]]: (SHARC 设置, 头数据, 轨迹数据) , 例如 

        .. code-block:: 

            sharc_setting = {
            'SHARC_version': 3.0,
            'method': 0,
            'integrator': 2,
            'maxmult': 2,
            'nstates_m': [0, 5],
            'natom': 3,
            'dtstep': 20.6706868947804,
            'nsteps': 1000,
            'nsubsteps': 25,
            'ezero': -185.2605928,
            'write_overlap': 1,
            'write_grad': 0,
            'write_nacdr': 0,
            'write_property1d': 0,
            'write_property2d': 0,
            'n_property1d': 1,
            'n_property2d': 1,
            'laser': 0
            }

            traj_data = {
             'Atomic numbers': ['0.6000000000000E+001',
              '0.8000000000000E+001',
              '0.8000000000000E+001'],
             'Elements': ['C', 'O', 'O'],
             'Atomic masses': ['0.2187466181995E+005',
              '0.2915694637199E+005',
              '0.2915694637199E+005'],
             '0 Step': '3',
             '1 Hamiltonian (MCH) in a.u.': [['0.1458710000000E-001',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               '0.0000000000000E+000',
               ...
    """
    is_settings = True
    is_header_array_data = False
    # 原子质量等不会更新，其它内容会更新
    sharc_setting = dict()
    header_data = dict()
    data = dict()
    traj_data = []

    name = ''
    value = []

    with open(output_data_file, 'r') as f:
        line = f.readline()
        # 读取设置
        while line:
            if 'End of settings' in line:
                is_settings = False
                is_header_array_data = True
                # 跳过注释行
                line = f.readline()
                break

            if is_settings:
                line_split = line.split()
                setting_name = line_split[0]
                setting_value = [eval(v) for v in line_split[1:]]
                sharc_setting[setting_name] = setting_value[0] if len(
                    setting_value) == 1 else setting_value

            line = f.readline()

        while line:
            if 'End of header array data' in line:
                # 最后一段数据
                header_data[name] = value
                is_header_array_data = False
                # 跳过注释行
                line = f.readline()
                break

            if is_header_array_data:
                if line[0] == '!':
                    # 跳过第一个空白
                    if name:
                        header_data[name] = value
                    name = line[1:].strip()
                    value = []
                else:
                    line_split = line.split()
                    if 'Elements' not in name:
                        line_split = list(map(eval, line_split))
                    value.append(line_split[0] if len(line_split) ==
                                 1 else line_split)

            line = f.readline()

        # 重置
        name = ''
        value = []
        # 读取以 `!` 标记名称的数组数据
        while line:
            if line[0] == '!':
                if name:
                    data[name] = value[0] if len(value) == 1 else value

                    if "Step" in line:
                        traj_data.append(data)
                        data = dict()

                # 读取下一个数据的名称并初始化值
                name = line[1:].strip()
                value = []
            # 数据行
            else:
                line_split = line.split()
                # line_split = list(map(eval, line_split))
                try:
                    line_split = [eval(v) for v in line_split]
                except NameError:
                    # 'NaN' 会转换失败
                    pass
                value.append(line_split[0] if len(line_split) ==
                             1 else line_split)

            line = f.readline()

    # 把最后一步的数据加上
    traj_data.append(data)
    return sharc_setting, header_data, traj_data


def show_mol(xyzfile, width:int=400,height:int=400):
    with open(xyzfile,'r') as f:
        mol_block = f.read()
    f.close()
    
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModelsAsFrames(mol_block, 'xyz')
    viewer.setStyle({'stick': {'radius': 0.1,'colorscheme':'Jmol'}, 'sphere': {'radius': 0.3,'colorscheme':'Jmol'}})
    viewer.setBackgroundColor('black')
    viewer.animate({'loop': 'forward',})
    viewer.zoomTo()
    
    components.html(viewer._make_html(), width=width, height=height)
    # st.write('缩放: 滚动鼠标滑轮, 旋转: 拖动鼠标左键, 平移: 拖动鼠标滑轮')
    
def show_energy_curve(expec_out_file, incoord_data = None, width:int=400,height:int=400):
    try:
        with open(expec_out_file,'r') as f:
            f.readline()
            title_row = f.readline()
        f.close()
        names=[col_header.strip().replace(' ','') for col_header in title_row[1:].split('|')[:-1]]
        
        data = pd.read_csv(expec_out_file,delimiter='\s+',skiprows=3,names=names)
        
        options = st.multiselect(
            label = '选择需要画图的数据',
            options = names,
            # default = ['Trajectory'],
            placeholder = '可选择多个',
            # help = 'X轴为Time, 默认画Trajectory',
            )
        if incoord_data is None:
            fig = go.Figure()
        else:
            # Create figure with secondary y-axis
            fig = make_subplots(specs=[[{"secondary_y": True}]])
            fig.add_trace(
                go.Scatter(x=incoord_data[0], y=incoord_data[1], name=incoord_data[2]),
                secondary_y=True,
            )
            
        fig.add_trace(go.Scatter(x=data['Time'], y=data['Epot'],
                mode='markers',
                line=dict(color='black'),
                name='Trajectory'))
        for key in options:
            fig.add_trace(go.Scatter(x=data['Time'], y=data[key],
                mode='lines',
                name=key))

            
        fig.update_layout(xaxis=dict(title='Time (fs)'),yaxis=dict(title='Energies in diagonal basis (eV)'))
        st.plotly_chart(fig, use_container_width=True)
    except:
        st.warning(f'Failed to extract info from {os.path.dirname(expec_out_file)}')
        
        
def show_file(xyzfile):
    file_type_display = st.text_input('文件名:',value='input')
    path_temp = os.path.dirname(xyzfile)+f"/{file_type_display}"
    if os.path.exists(path_temp):
        suffix = file_type_display.split('.')[-1]
        if suffix in ['jpg','png']:
            st.image(path_temp)
        else:
            with open(path_temp,'r') as f:
                content = f.read()
            f.close()
            
            st.text(content)
    else:
        st.warning(f'文件{path_temp}不存在!')
    

def cal_incoord(xyzfile, *idxs):
    time = None 
    incoord = None 
    incoord_name = ''
    try:
        idxs = idxs[0]
        traj = xyztraj.Traj(xyzfile)
        traj.parse_title(lambda x: float(x.split()[1]))
        time = traj.title_info
        natom = len(idxs)
        # st.write(natom)
        # st.write(idxs)
        
        if natom == 2:
            incoord = traj.bond(*idxs) 
            incoord_name = "Bond %d-%d (Å)"%tuple(idxs)
        elif natom == 3:
            incoord = traj.angle(*idxs) 
            incoord_name = "Angle %d-%d-%d (°)"%tuple(idxs)
        elif natom == 4:
            incoord = traj.dihedral(*idxs)  
            incoord_name = "Dihedral %d-%d-%d-%d (°)"%tuple(idxs)
        else:
            pass
        return time, incoord, incoord_name
    except:
        # print(f"Failed to calculate internal coordinate of {xyzfile}.")
        logger.logger.error((f"Failed to calculate internal coordinate of {xyzfile}."))
        return None
        
    
# ================================================================
# 获取最大步数
def get_max_nstep(filename="KEYSTROKES.setup_traj"):
    with open(filename, 'r') as f:
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

# ================================================================
# KER
# def plot_KER(df):
# Define the rgbcolor class
class rgbcolor:
    def __init__(self, initlist):
        excluded = [[0.12, 0.22]]  # exclude yellow hues from the colorwheel
        excluded.sort(key=lambda x: x[0])
        temp1 = [[min(1., max(0., el[0])), max(0., min(1., el[1]))] for el in excluded]
        temp2 = [[0., 0.]]
        for el in temp1:
            if el[0] >= temp2[-1][1]:
                temp2.append(el)
            else:
                temp2[-1][1] = el[1]
        self.excluded = temp2[1:]
        self.initlist = [max(0, el) for el in initlist]
        self.n = sum(el > 0 for el in self.initlist)
        self.m = len(initlist)
        self.a = 1. - sum(el[1] - el[0] for el in self.excluded)
        self.startlist = [0.] * self.m
        self.incrlist = [0.] * self.m
        for i in range(1, self.m):
            self.startlist[i] = self.startlist[i - 1] + (self.a / self.n if self.initlist[i - 1] > 0 else 0)
        for i in range(self.m):
            if self.initlist[i] > 0:
                self.incrlist[i] = self.a / self.n / self.initlist[i]

    def rgb_to_hex(self, rgb):
        return '#{:02x}{:02x}{:02x}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))

    def hexcolor(self, index, el):
        if not (1 <= index <= self.m and 1 <= el <= self.initlist[index - 1]):
            return '#FFFFFF'
        hue = self.startlist[index - 1] + self.incrlist[index - 1] * (el - 1)
        for start, end in self.excluded:
            if hue > start:
                hue += end - start
        return self.rgb_to_hex(colorsys.hsv_to_rgb(hue, 1, 1))


def generate_colors(num:int=10,alpha:float=0.4,rgba_format:str='plotly'):
    """rgba_format
        plotly 'rgba(255,255,255,0.5)'
        matplotlib (0.1,0.5,0.9,0.5)
    """
    def hex2rgb(hex_str):
        # return tuple(int(hex_str.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
        return list(int(hex_str.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
    
    def rgba2hex(rgb_tuple, alpha):
        rgba_tuple = (rgb_tuple[0]/255,rgb_tuple[1]/255,rgb_tuple[2]/255,alpha)
        return mpl.colors.rgb2hex(rgba_tuple, keep_alpha=True)
        
    color_generator = rgbcolor([1] * num)
    colors_generated = [color_generator.hexcolor(i + 1, 1) for i in range(num)]
    # return [rgba2hex(hex2rgb(c),alpha) for c in colors_generated]
    if rgba_format == 'plotly':
        return [ 'rgba'+str(tuple(hex2rgb(c) + [alpha])) for c in colors_generated]
    elif rgba_format == 'matplotlib':
        return [ tuple([value/255 for value in hex2rgb(c)] + [alpha]) for c in colors_generated]
    else:
        raise ValueError("Wrong rgba_format")

def plot_KER_state(df,bin_count = 30, kde_bandwidth = 0.2 ):
    # 用户可配置参数
    # plot_range = (0, 16)
    # bin_count = 30
    # kde_bandwidth = 0.2
    opacity = 0.3
    
    # 加载数据
    # df = pd.read_csv('2024-11-15T08-01_export.csv')
    plot_range = (df['delta_energy'].min(), df['delta_energy'].max()) 
    
    color_generator = rgbcolor([1] * len(df['state'].unique()))
    state_colors = [color_generator.hexcolor(i + 1, 1) for i in range(len(df['state'].unique()))]
    
    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4), gridspec_kw={'height_ratios': [0.7, 0.3]}, sharex=True)
    
    y_kde_sum = np.zeros(1000)
    x_range = np.linspace(plot_range[0], plot_range[1], 1000)
    hist_max = 0
    
    for state, color in zip(df['state'].unique(), state_colors):
        subset = df[df['state'] == state]['delta_energy']
        hist, edges = np.histogram(subset, bins=np.linspace(plot_range[0], plot_range[1], bin_count))
        ax2.bar(edges[:-1], hist, width=np.diff(edges)[0], alpha=opacity, color=color, label=state)
        if max(hist) > hist_max:
            hist_max = max(hist) 
        
        kde = gaussian_kde(subset, bw_method=kde_bandwidth)
        y_kde = kde.evaluate(x_range)
        ax1.fill_between(x_range, 0, y_kde, alpha=opacity, color=color, label=state)
        y_kde_sum += y_kde
    
    # 归一化并绘制总和曲线
    y_kde_sum_normalized = y_kde_sum / y_kde_sum.max()
    ax1.plot(x_range, y_kde_sum_normalized, color='black', label='Normalized Total Sum', linewidth=2)
    
    # 设置图表的整体属性
    ax1.set_ylabel('Intensity (arb. units)')
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False,fontsize=8)
    
    ax2.set_ylabel('Count')
    ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel('KER (eV)')
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白
    
    # plt.show()
    return fig
    
def plot_KER_state2(df, sigma, alpha=0.3, num_point=1000):

    def gaussian(x, sigma=1, mu=0):
        gx = np.exp(-(x - mu)**2 / sigma**2 / 2) / (sigma * np.sqrt(2 * np.pi))
        return gx

    plot_range = (df['delta_energy'].min(), df['delta_energy'].max())

    color_generator = rgbcolor([1] * len(df['state'].unique()))
    generated_colors = [
        color_generator.hexcolor(i + 1, 1)
        for i in range(len(df['state'].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(2,
                                   1,
                                   figsize=(6, 4),
                                   gridspec_kw={'height_ratios': [0.7, 0.3]},
                                   sharex=True)

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
        ax1.fill_between(x_range,
                         0,
                         plot_data[plot_type]["curve_data"] / ratio,
                         alpha=alpha,
                         color=generated_colors[counter],
                         label=plot_type)
        ax2.plot(plot_data[plot_type]["point_data"],
                 np.zeros(len(plot_data[plot_type]["point_data"])) - counter,
                 marker='|',
                 lw=0,
                 color=generated_colors[counter],
                 alpha=alpha)
        counter += 1
    # 归一化并绘制总和曲线

    ax1.plot(x_range,
             y_sum / ratio,
             color='black',
             label='Normalized Total Sum',
             linewidth=2)

    # 设置图表的整体属性
    ax1.set_ylabel('Intensity (arb. units)')
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False, fontsize=8)

    # ax2.set_ylabel('Count')
    ax2.set_yticks([])
    # ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel('KER (eV)')
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig

    
def plot_KER_frag(df,bin_count = 30, kde_bandwidth = 0.2 ):
    # 用户可配置参数
    # plot_range = (0, 16)
    # bin_count = 30
    # kde_bandwidth = 0.2
    opacity = 0.3
    
    # 加载数据
    # df = pd.read_csv('2024-11-15T08-01_export.csv')
    plot_range = (df['delta_energy'].min(), df['delta_energy'].max()) 
    
    color_generator = rgbcolor([1] * len(df['frag_string'].unique()))
    frag_string_colors = [color_generator.hexcolor(i + 1, 1) for i in range(len(df['frag_string'].unique()))]
    
    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4), gridspec_kw={'height_ratios': [0.7, 0.3]}, sharex=True)
    
    y_kde_sum = np.zeros(1000)
    x_range = np.linspace(plot_range[0], plot_range[1], 1000)
    hist_max = 0
    
    for frag_string, color in zip(df['frag_string'].unique(), frag_string_colors):
        subset = df[df['frag_string'] == frag_string]['delta_energy']
        hist, edges = np.histogram(subset, bins=np.linspace(plot_range[0], plot_range[1], bin_count))
        ax2.bar(edges[:-1], hist, width=np.diff(edges)[0], alpha=opacity, color=color, label=frag_string)
        if max(hist) > hist_max:
            hist_max = max(hist) 
        
        kde = gaussian_kde(subset, bw_method=kde_bandwidth)
        y_kde = kde.evaluate(x_range)
        ax1.fill_between(x_range, 0, y_kde, alpha=opacity, color=color, label=frag_string)
        y_kde_sum += y_kde
    
    # 归一化并绘制总和曲线
    y_kde_sum_normalized = y_kde_sum / y_kde_sum.max()
    ax1.plot(x_range, y_kde_sum_normalized, color='black', label='Normalized Total Sum', linewidth=2)
    
    # 设置图表的整体属性
    ax1.set_ylabel('Intensity (arb. units)')
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False,fontsize=8)
    
    ax2.set_ylabel('Count')
    ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel('KER (eV)')
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白
    
    # plt.show()
    return fig
    

def plot_KER_state2(df, sigma, alpha=0.3, num_point=1000):

    def gaussian(x, sigma=1, mu=0):
        gx = np.exp(-(x - mu)**2 / sigma**2 / 2) / (sigma * np.sqrt(2 * np.pi))
        return gx

    plot_range = (df['delta_energy'].min(), df['delta_energy'].max())

    color_generator = rgbcolor([1] * len(df['frag_string'].unique()))
    generated_colors = [
        color_generator.hexcolor(i + 1, 1)
        for i in range(len(df['frag_string'].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(2,
                                   1,
                                   figsize=(6, 4),
                                   gridspec_kw={'height_ratios': [0.7, 0.3]},
                                   sharex=True)

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
        ax1.fill_between(x_range,
                         0,
                         plot_data[plot_type]["curve_data"] / ratio,
                         alpha=alpha,
                         color=generated_colors[counter],
                         label=plot_type)
        ax2.plot(plot_data[plot_type]["point_data"],
                 np.zeros(len(plot_data[plot_type]["point_data"])) - counter,
                 marker='|',
                 lw=0,
                 color=generated_colors[counter],
                 alpha=alpha)
        counter += 1
    # 归一化并绘制总和曲线

    ax1.plot(x_range,
             y_sum / ratio,
             color='black',
             label='Normalized Total Sum',
             linewidth=2)

    # 设置图表的整体属性
    ax1.set_ylabel('Intensity (arb. units)')
    ax1.set_ylim(0.01, 1.05)
    ax1.legend(frameon=False, fontsize=8)

    # ax2.set_ylabel('Count')
    ax2.set_yticks([])
    # ax2.set_ylim(0, hist_max + 2)
    ax2.set_xlabel('KER (eV)')
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig
    
    
def plot_KER(df, col_name = "state", sigma = 1, alpha=0.3, num_point=1000):
    # FWHM
    # sigma = FWHM / (2*np.sqrt(2*np.log(2)))

    def gaussian(x, sigma=1, mu=0):
        gx = np.exp(-(x - mu)**2 / sigma**2 / 2) / (sigma * np.sqrt(2 * np.pi))
        return gx

    plot_range = (df['delta_energy'].min(), df['delta_energy'].max())

    color_generator = rgbcolor([1] * len(df[col_name].unique()))
    generated_colors = [
        color_generator.hexcolor(i + 1, 1)
        for i in range(len(df[col_name].unique()))
    ]

    # 创建图形并调整子图比例，设置sharex=True以共享x轴
    fig, (ax1, ax2) = plt.subplots(2,
                                   1,
                                   figsize=(6, 4),
                                   gridspec_kw={'height_ratios': [0.7, 0.3]},
                                   sharex=True)

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
        ax1.fill_between(x_range,
                         0,
                         plot_data[plot_type]["curve_data"] / ratio,
                         alpha=alpha,
                         color=generated_colors[counter],
                         label=plot_type)
        ax2.plot(plot_data[plot_type]["point_data"],
                 np.zeros(len(plot_data[plot_type]["point_data"])) - counter,
                 marker='|',
                 ms = 6,
                 lw=0,
                 color=generated_colors[counter],
                 alpha=alpha)
        counter += 1
    # 归一化并绘制总和曲线

    ax1.plot(x_range,
             y_sum / ratio,
             color='black',
             label='Normalized Total Sum',
             linewidth=2)

    # 设置图表的整体属性
    ax1.set_ylabel('Intensity (arb. units)')
    ax1.set_ylim(0., 1.05)
    ax1.legend(frameon=False, fontsize=8)

    # ax2.set_ylabel('Count')
    ax2.set_yticks([])
    ax2.set_ylim(-counter , 0.5)
    ax2.set_xlabel('KER (eV)')
    fig.subplots_adjust(hspace=0)  # 移除子图之间的空白

    # plt.show()
    return fig
    

# def dnplot(xyzfiles:list[str], orders:list[int]=[0,1,2], xylim: tuple[float]=(-1.,1.,-1.,1.)) -> tuple:
    # datfiles = [xf.replace('output.xyz','output.dat') for xf in xyzfiles]
    # try:
        # elements = ptf.read_elements(datfiles[0])
        # atomic_masses = ptf.read_atomic_masses(datfiles[0])
    # except:
        # elements = ptf.read_elements(datfiles[-1])
        # atomic_masses = ptf.read_atomic_masses(datfiles[-1])
        
    # kins = []
    # monts = []
    
    # for f in datfiles:
        # try:
            # velocities = ptf.read_velocities(f)
            
            # p = np.array([m*v for m,v in zip(atomic_masses,velocities)])
            # e = np.array([0.5*m*np.sum(v**2) for m,v in zip(atomic_masses,velocities)])
            
            # monts.append(p)
            # kins.append(e)
        # except:
            # print(f"Failed to read file {f}")
            
    # np.save('ps.npy',monts)
    # np.save('es.npy',kins)
    
    # kins = np.array(kins)
    # monts = np.array(monts)
    
    # elements = ["%s_{%d}"%(e,i+1) for i,e in enumerate(elements)]
    
    # return ptf.dalitzplot_rect(kins, elements, orders), ptf.newtonplot(monts, elements, orders)
    
def dnplot(xyzfiles:list[str], orders:list[int]=[0,1,2], xylim: tuple[float]=(-1.,1.,-1.,1.), idx : int = 0) -> tuple:
    datfiles = [xf.replace('output.xyz','output.dat') for xf in xyzfiles]
    elements = [str(i) for i in range(3)]
    kins = []
    monts = []
    
    for f in datfiles:
        try:
            elements,atomic_masses,velocities = ptf.read_output_dat(f, idx)
            # print(elements,atomic_masses,velocities)
            
            p = np.array([m*v for m,v in zip(atomic_masses,velocities)])
            e = np.array([0.5*m*np.sum(v**2) for m,v in zip(atomic_masses,velocities)])
            # print(p,e)
            monts.append(p)
            kins.append(e)
        except:
            print(f"Failed to read file {f} or calculate Ekin and p")
            
    kins = np.array(kins)
    monts = np.array(monts)
    
    np.save('ps.npy',monts)
    np.save('es.npy',kins)
    
    elements = ["%s_{%d}"%(e,i+1) for i,e in enumerate(elements)]
    
    return ptf.dalitzplot_rect(kins, elements, orders, 'top'), ptf.newtonplot(monts, elements, orders, xylim, 'top')
    
# def dnplot(xyzfiles:list[str]=glob("*let_*/TRAJ_*/output.xyz"), orders:list[int]=[0,1,2], xylim: tuple[float]=(-1.,1.,-1.,1.), idx : int = 0, hdf5_file:str='Data.hdf5') -> tuple:
    # f = h5py.File(hdf5_file, "r")
    # elements = [str(i) for i in range(3)]
    # kins = []
    # monts = []
    
    # for path in (os.path.dirname(xf) for xf in xyzfiles):
        # try:
            # elements = f[path].attrs['elements']
            # atomic_masses = f[path].attrs['atomic_masses']
            # velocities = f[path]['velocity'][idx, :]
            
            # # elements,atomic_masses,velocities = ptf.read_output_dat(f, idx)
            # # print(elements,atomic_masses,velocities)
            
            # p = np.array([m*v for m,v in zip(atomic_masses,velocities)])
            # e = np.array([0.5*m*np.sum(v**2) for m,v in zip(atomic_masses,velocities)])
            # # print(p,e)
            # monts.append(p)
            # kins.append(e)
        # except:
            # print(f"Failed to read file {f} or calculate Ekin and p")
            
    # kins = np.array(kins)
    # monts = np.array(monts)
    
    # # print(kins.shape)
    
    # if len(kins) > 1:
        
    
        # np.save('ps.npy',monts)
        # np.save('es.npy',kins)
        
        # elements = ["%s_{%d}"%(e,i+1) for i,e in enumerate(elements)]
    
        # return ptf.dalitzplot_rect(kins, elements, orders, 'top'), ptf.newtonplot(monts, elements, orders, xylim, 'top')
    # else:
        # print(f"No data to plot.")
        # return None, None
    
# dnplot()