
'''
version: 2024-10-21  Feat : Filter of Multiplicity.
version: 2024-10-16  Feat : Statistics on fragmentation information
version: 2024-10-10  Fix : find and read template file
'''

VERSION = '2024.10.21'

import os
import json
import numpy as np
import py3Dmol
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from glob import glob

st.set_page_config(page_title='MolFragApp', layout='wide')

def read_template(filename='MOLCAS.template', output_type='str'):
    assert os.path.exists(filename)

    if filename == 'MOLCAS.template':
        qc_software = 'OpenMolcas'

    with open(filename, 'r') as f:
        lines = f.readlines()
    f.close()

    template = {'qc_software': qc_software}
    for line in lines:
        temp = line.split()
        if len(temp) >= 2 and '#' not in line:
            template[temp[0]] = temp[1]

    if output_type == 'str':
        return f"SA({template['roots']})-{template['method']}({template['nactel']}e,{template['ras2']}o)/{template['basis']}, {template['qc_software']}"

    return template

def get_delta_kin(lisfile):
    # lis只有一步会报错
    lis_data = np.loadtxt(lisfile)
    delta_kin = lis_data[-1,4] - lis_data[0,4] # in eV 
    return delta_kin


def parse_coord(coord_str_list: list):
    coord_arr = np.array([row.split()[1:] for row in coord_str_list],
                         dtype='float64')
    return coord_arr


def read_xyz(xyzfile):
    with open(xyzfile, 'r') as f:
        lines = f.readlines()
    f.close()

    num_atom = int(lines[0])
    coords = []
    for i in range(2, len(lines), num_atom + 2):
        coords.append(parse_coord(lines[i:i + num_atom]))
    return np.array(coords)


def gen_chemical_formula(atom_list: list = ['H', 'H', 'O']):
    ele = []
    num = []

    for atom in atom_list:
        if atom not in ele:
            ele.append(atom)
            num.append(1)
        else:
            num[ele.index(atom)] += 1

    chemical_formula = []
    for i in range(len(ele)):
        chemical_formula.append(ele[i] + (str(num[i]) if num[i] > 1 else ''))

    return ''.join(sorted(chemical_formula))


class XYZ:

    def __init__(self, xyzfile) -> None:
        self.xyzfile = xyzfile

    def read_xyz(self):
        with open(self.xyzfile, 'r') as f:
            lines = f.readlines()
        f.close()

        try:
            num_atom = int(lines[0])
            self.num_atom = num_atom
        except:
            raise ValueError(f"Can NOT read file {self.xyzfile}!")

        atom_list = [lines[i].split()[0] for i in range(2, 2 + num_atom)]
        self.atom_list = atom_list

        coords = []
        for i in range(2, len(lines), num_atom + 2):
            coords.append(parse_coord(lines[i:i + num_atom]))

        self.coords = np.array(coords)

    def get_frag(self,
                 frame_idx: int = -1,
                 max_bond_length: float = 2.5,
                 min_atom_num: int = 1):

        from sklearn.cluster import DBSCAN
        clustering = DBSCAN(eps=max_bond_length, min_samples=min_atom_num).fit(
            self.coords[frame_idx])
        n_clusters = max(clustering.labels_) + 1
        frag = dict()
        for i in range(n_clusters):
            frag[i] = {'atom': [], 'number': []}
        for i, idx_cluster in enumerate(clustering.labels_):
            frag[idx_cluster]['atom'].append(self.atom_list[i])
            frag[idx_cluster]['number'].append(i + 1)

        self.frag = frag

    def get_formula(self):
        for idx_frag in self.frag:
            self.frag[idx_frag]['formula'] = gen_chemical_formula(
                self.frag[idx_frag]['atom'])

   
def show_mol(xyzfile,width:int=400,height:int=400,show_mol=True,show_energy=True):
    with open(xyzfile,'r') as f:
        mol_block = f.read()
    f.close()
    
    col1, col2 = st.columns(2)
    
    # with col1:
    if show_mol:
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModelsAsFrames(mol_block, 'xyz')
        viewer.setStyle({'stick': {'radius': 0.1,'colorscheme':'Jmol'}, 'sphere': {'radius': 0.3,'colorscheme':'Jmol'}})
        viewer.setBackgroundColor('black')
        viewer.animate({'loop': 'forward',})
        viewer.zoomTo()
        
        components.html(
            viewer._make_html(),
            width=width,
            height=height
            )
        
        # st.write('缩放: 滚动鼠标滑轮, 旋转: 拖动鼠标左键, 平移: 拖动鼠标滑轮')
            
    # with col2:
    if show_energy:
        res = os.system(f"cd {os.path.dirname(xyzfile)}; $SHARC/data_extractor.x output.dat; ")
    
        if res == 0:
            # st.write(os.path.dirname(xyzfile))
            with open(f"{os.path.dirname(xyzfile)}/output_data/expec.out",'r') as f:
                f.readline()
                title_row = f.readline()
            f.close()
            names=[col_header.strip().replace(' ','') for col_header in title_row[1:].split('|')[:-1]]
            
            # 调试
            # st.write(names)
            
            data = pd.read_csv(f"{os.path.dirname(xyzfile)}/output_data/expec.out",delimiter='\s+',skiprows=3,names=names)
            
            options = st.multiselect(
                label = '选择需要画图的数据',
                options = names,
                # default = ['Trajectory'],
                placeholder = '可选择多个',
                # help = 'X轴为Time, 默认画Trajectory',
                )
                
            # data_selected = {"Time":data['Time'],'Trajectory':data['Epot'],'Etot':data['Etot']}
            # for key in options:
            #     data_selected[key] = data[key]
            # data_plot = pd.DataFrame(data_selected)
            # st.line_chart(data_plot,x="Time",y=options+["Trajectory"],x_label = 'Time [fs]',y_label= 'Energies in diagonal basis [eV]',height=350,width=800)

            fig = go.Figure()
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



@st.cache_data
def frag_analysis(xyz_path, max_bond_length: float = 2.5, min_atom_num: int = 1):
    xyzfiles = sorted(glob(xyz_path))
    xyzs = []
    steps = []
    frags = []

    delta_kin = []

    # st.write(max_bond_length,min_atom_num)

    progress_text = "Operation in progress. Please wait."
    my_bar = st.progress(0, text=progress_text)
    counter = 0
    num_xyzfile = len(xyzfiles)

    for xyzfile in xyzfiles:
        try:            
            xyz = XYZ(xyzfile)
            xyz.read_xyz()
            
            xyz.get_frag(frame_idx=-1,
                         max_bond_length=max_bond_length,
                         min_atom_num=min_atom_num)
            xyz.get_formula()
            temp = []
            for idx_frag in xyz.frag:
                temp.append(xyz.frag[idx_frag]['formula'])
                
            # frags.append(sorted(temp))
            for i in range(1,len(temp)):
                for j in range(0,len(temp)-i):
                    if len(temp[j]) < len(temp[j+1]) :
                        temp[j] , temp[j+1] = temp[j+1] , temp[j]
                        
            frags.append(temp)

            
            xyzs.append(xyzfile)
            steps.append(len(xyz.coords))

            try:
                dk = get_delta_kin(f"{os.path.dirname(xyzfile)}/output.lis")
                delta_kin.append(dk)
            except:
                print(f"Can NOT get delta kin from file {xyzfile}!")
                delta_kin.append(np.nan)

        except KeyError:
            print(f"`max_bond_length` and `min_atom_num` of {xyzfile} are inappropriate")
            # st.warning(f'max_bond_length and min_atom_num of {xyzfile} are inappropriate')
            # frags.append([None])
            continue
        
        except ValueError:
            print(f'Content of {xyzfile} is inappropriate')
            continue

        counter += 1
        my_bar.progress(counter / num_xyzfile, text=progress_text)

    my_bar.empty()

    # st.write(frags)

    return xyzs, steps, frags, delta_kin


st.title(f'分子解离片段分析 v{VERSION}')

st.write(f'当前路径:`{os.getcwd()}`')

path1, path2 = st.columns(2)

with path1:
    xyz_path = st.text_input(
        "轨迹文件(`xyz`格式)路径:",
        value = "Singlet_*/TRAJ_*/output.xyz",
        help = '使用`xyz`文件最后一帧的结构分析解离片段的成分',
    )
    num_xyzfile = len(glob(xyz_path))
    st.write("轨迹文件总数：", num_xyzfile)

if num_xyzfile == 0:
    st.warning('No xyz file found!', icon="⚠️")

with path2:
    template_path = st.text_input(
        "`template`文件路径:",
        value = "*.template",
    )
    try:
        cal_method = read_template(glob(template_path)[0])
    except:
        cal_method = f"Warning: Can NOT read template from {template_path}!"
        
    st.write(f"计算级别:`{cal_method}`")

# st.write('**片段划分参数设置:**')

para1, para2 = st.columns(2)


with para1:
    max_bond_length = st.number_input(
        "划分片段时的最大键长(Å)",
        value=2.5,min_value=0.1,
        help='若两原子之间距离大于最大键长，则被划分到不同片段',
    )

with para2:
    delta_kin_shift = st.number_input(
        "$\Delta E_{kin}$ (eV)",
        value=0., min_value = 0.,
        help='',
    )

    # st.write(delta_kin_shift)


# st.divider()

st.header('解离片段分析结果')

xyzs, steps, frags, delta_kin = frag_analysis(xyz_path, max_bond_length)
num_frag = [len(f) for f in frags]
max_num_frag = max(num_frag)

# st.write(delta_kin)
# st.write(np.array(delta_kin))

# # np.save('frag.npy',(steps,num_frag,frags))
# with open('frag.json','w') as f:
    # json.dump(dict(steps=steps,nfrag=num_frag,frags=frags),f)
# f.close()

data = {
    "plot": [False]*len(xyzs),
    "xyz": xyzs,
    "steps": steps,
    "nfrag": num_frag,
    "delta_kin": delta_kin,
}
# data = {"xyz file": [os.path.abspath(f) for f in xyzfiles],"steps": steps,}
for idx_frag in range(max_num_frag):
    data[f'frag-{idx_frag+1}'] = []
    for frag_list in frags:
        if idx_frag < len(frag_list):
            data[f'frag-{idx_frag+1}'].append(frag_list[idx_frag])
        else:
            data[f'frag-{idx_frag+1}'].append(None)

data["delta_kin_shifted"] = np.array(delta_kin) - delta_kin_shift
data["delta_kin_distribution"] = [False]*len(xyzs)

df = pd.DataFrame(data)

st.markdown('可按照某列的值进行排序、搜索；选中某一行观看动画')


col1, col2 = st.columns(2)

with col1:
    edited_df = st.data_editor(
        data = df,
        height=1000,
        hide_index = True,
        column_config = {
            "plot" : st.column_config.CheckboxColumn(
                width = 'small',
                help = '选择轨迹查看动画和能量曲线',
            ),
            "delta_kin_shifted" : st.column_config.NumberColumn(
                width = 'small',
            ),
            "delta_kin_distribution" : st.column_config.CheckboxColumn(
                width = 'small',
                help = '选择多条轨迹画直方图',
            ),
        },
        disabled = ['xyz','steps','nfrag','delta_kin','delta_kin_shifted']+[f'frag-{i+1}' for i in range(max_num_frag)],
    )


with col2:
    tab1, tab2, tab3, tab4 = st.tabs(["原始轨迹筛选与统计","轨迹动画&能量曲线", "能量差直方图", "查看其它文件"])

    with tab1:
        frags_info = dict(filepath=xyzs, steps=steps,nfrag=num_frag,frags=frags)
        
        frags_type = []
        frags_num = []

        para1_filter, para2_filter, para3_filter = st.columns(3)
        
        with para1_filter:
            min_steps = st.number_input(
            "轨迹的最小步数",
            value=1000,min_value=1,
        )
        
        with para2_filter:
            min_nfrag = st.number_input(
            "片段数目最小值",
            value=2,min_value=1,max_value=max(frags_info['nfrag'])-1
        )
        
        with para3_filter:
            max_nfrag = st.number_input(
            "片段数目最大值",
            value=3,min_value=2,max_value=max(frags_info['nfrag'])
        )
        
        all_multiplicity = []
        for file_path in frags_info['filepath']:
            multi = file_path.split('/')[0]
            if multi not in all_multiplicity:
                all_multiplicity.append(multi)                
            
        
        multiplicity = st.multiselect(
                label = '轨迹的初始多重度',
                options = all_multiplicity,
                default = all_multiplicity,
            )
        
        # frag_chosen = st.text_input(
        # "是否包含片段",
        # )
        
        multiplicity_filtered = []
        for i in range(len(frags_info['steps'])):
            # 根据设定的条件筛选轨迹
            multi_temp = frags_info['filepath'][i].split('/')[0]
            if frags_info['steps'][i] >= min_steps and min_nfrag <= frags_info['nfrag'][i] <= max_nfrag and multi_temp in multiplicity: 
                multiplicity_filtered.append(multi_temp)
                frags_filter = frags_info['frags'][i]
                if frags_filter in frags_type:
                    idx = frags_type.index(frags_filter)
                    frags_num[idx] += 1
                else:
                    frags_type.append(frags_filter)
                    frags_num.append(1)
        
        # frags_type,frags_num,sum(frags_num)
        # 降序排序
        for i in range(1,len(frags_num)):
            for j in range(0,len(frags_num)-i):
                if frags_num[j] < frags_num[j+1]:
                    frags_num[j], frags_num[j+1] = frags_num[j+1], frags_num[j]
                    frags_type[j], frags_type[j+1] = frags_type[j+1], frags_type[j]
                    
        multiplicity_filtered_info = {'type':[],'num':[]}
        for mf in multiplicity_filtered:
            if mf not in multiplicity_filtered_info['type']:
                multiplicity_filtered_info['type'].append(mf)
                multiplicity_filtered_info['num'].append(1)
            else:
                idx_mf = multiplicity_filtered_info['type'].index(mf)
                multiplicity_filtered_info['num'][idx_mf] += 1

        for m in range(1,len(multiplicity_filtered_info['num'])):
            for n in range(0,len(multiplicity_filtered_info['num'])-m):
                if multiplicity_filtered_info['num'][n] < multiplicity_filtered_info['num'][n+1]:
                    multiplicity_filtered_info['num'][n], multiplicity_filtered_info['num'][n+1] = multiplicity_filtered_info['num'][n+1], multiplicity_filtered_info['num'][n]
                    multiplicity_filtered_info['type'][n], multiplicity_filtered_info['type'][n+1] = multiplicity_filtered_info['type'][n+1], multiplicity_filtered_info['type'][n]

        # st.write(multiplicity_filtered_info)

        multiplicity_filtered_info_strlist = [f"{multiplicity_filtered_info['num'][i]} {multiplicity_filtered_info['type'][i]}" for i in range(len(multiplicity_filtered_info['num']))]
        # st.write(multiplicity_filtered_info_string)
        st.write(f"{sum(frags_num)}条轨迹({', '.join(multiplicity_filtered_info_strlist)}), {len(frags_type)}种解离路径 ")

        fig = px.bar(dict(frags_type=[" + ".join(f) for f in frags_type],frags_num=frags_num), x='frags_type',y='frags_num')
        fig.update_layout(xaxis=dict(title='片段类型'),yaxis=dict(title='片段数目'))
        st.plotly_chart(fig, use_container_width=True)
        
        
        df_pie = pd.DataFrame(
            {
                "frags_type" : [" + ".join(f) for f in frags_type],
                "frags_num" : frags_num,
                "H" : [("H" if "H" in f else "No H") for f in frags_type],
                "H2" : [("H2" if "H2" in f else "No H2") for f in frags_type],
                "H channel" : [("H channel" if "H2" in f or "H" in f else "No H channel") for f in frags_type],
            }
        )
        fig = px.sunburst(df_pie, path=['H2', 'H', 'frags_type'], values='frags_num',color_discrete_sequence=['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52'])
        # fig = px.sunburst(
            # df_pie, 
            # path=['H channel', 'H2', 'frags_type'], 
            # values='frags_num',
            # color = 'frags_num',
            # # color_discrete_map={'(?)':'black', 'Lunch':'gold', 'Dinner':'darkblue'},
        # )
        st.plotly_chart(fig, use_container_width=True)
        
        

        frag_chosen1 = st.text_input('是否包含片段(1/3)',value=None)
        if frag_chosen1 is not None:
            frag_chosen1 = frag_chosen1.strip()
            df_pie[frag_chosen1] = [(frag_chosen1 if frag_chosen1 in f else "Other") for f in frags_type]
            
            fig1 = px.pie(df_pie,values='frags_num',names=frag_chosen1)
            
            st.plotly_chart(fig1, use_container_width=True)
            
        frag_chosen2 = st.text_input('⋙ 是否包含片段(2/3)',value=None)
        if frag_chosen2 is not None:
            frag_chosen2 = frag_chosen2.strip()
            df_pie[frag_chosen2] = [(frag_chosen2 if frag_chosen2 in f else "Other") for f in frags_type]
            
            option1 = st.radio(
                label = '选择片段(1/3)',
                options = [frag_chosen1,'Other'],
                horizontal = True,
                )
                
            df_pie2 = df_pie[ df_pie[frag_chosen1] == option1 ]
            
            fig2 = px.pie(df_pie2,values='frags_num',names=frag_chosen2,title=f"{frag_chosen1} = {option1}")
            
            st.plotly_chart(fig2, use_container_width=True)
            
        frag_chosen3 = st.text_input('⋙⋙ 是否包含片段(3/3)',value=None)
        if frag_chosen3 is not None:
            frag_chosen3 = frag_chosen3.strip()
            df_pie[frag_chosen3] = [(frag_chosen3 if frag_chosen3 in f else "Other") for f in frags_type]
            
            option2 = st.radio(
                label = '选择片段(2/3)',
                options = [frag_chosen2,'Other'],
                horizontal = True,
                )
                
            df_pie3 = df_pie[ df_pie[frag_chosen1] == option1 ]
            df_pie3 = df_pie3[ df_pie3[frag_chosen2] == option2 ]
            
            fig3 = px.pie(df_pie3,values='frags_num',names=frag_chosen3,title=f"{frag_chosen1} = {option1} / {frag_chosen2} = {option2}")
            
            st.plotly_chart(fig3, use_container_width=True)
                
            

    with tab2:
        xyzfiles_selected = edited_df[edited_df['plot']==True]['xyz']
        num_xyzfiles_selected = len(xyzfiles_selected)
        xyzfile_display = None
        if num_xyzfiles_selected == 0:
            st.info("未选择轨迹")
        elif num_xyzfiles_selected > 1:
            xyzfile_display = st.radio(
                label = '选择一个轨迹查看动画与能量曲线',
                options = xyzfiles_selected,
                horizontal = True,
            )
            xyzfile_display = xyzfile_display
        else:
            xyzfile_display = xyzfiles_selected.iloc[0]

        if xyzfile_display is not None:
            show_mol(xyzfile_display,width=850,height=350)
            
    with tab3:
        # st.write(edited_df)
        delta_kin_shifted = edited_df[edited_df['delta_kin_distribution']==True]['delta_kin_shifted']
        if len(delta_kin_shifted) > 1:
            bin_size = st.number_input("bin size", value=1., min_value = 0.1,)
            fig = ff.create_distplot(
                [delta_kin_shifted],
                ['KER'],
                bin_size=[bin_size],
            )
            fig.update_layout(xaxis=dict(title='Energy (eV)'))
            st.plotly_chart(fig, use_container_width=True)

            st.write('选中的轨迹:')
            st.dataframe(
                data = edited_df[edited_df['delta_kin_distribution']==True],
                use_container_width = True,
                hide_index = True,
                )
        else:
            st.info("未选择轨迹或选择轨迹数不足")
            
    with tab4:
        xyz_path_selected = edited_df[edited_df['plot']==True]['xyz']
        num_xyz_path_selected = len(xyz_path_selected)
        xyz_path_display = None
        if num_xyz_path_selected == 0:
            st.info("未选择轨迹")
        elif num_xyz_path_selected > 1:
            xyz_path_display = st.radio(
                label = '选择一个轨迹查看统一路径下的其它文件',
                options = xyz_path_selected,
                horizontal = True,
            )
            xyz_path_display = xyz_path_display
        else:
            xyz_path_display = xyz_path_selected.iloc[0]
        
        if xyz_path_display is not None:
            file_type_display = st.text_input('文件名:',value='input')
            path_temp = os.path.dirname(xyz_path_display)+f"/{file_type_display}"
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


st.divider()
components.html('''
<p align='center'>
<a href="https://github.com/ckz1/MolFragApp" target="_blank">
<img border="0" src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png" alt="MolFragApp in Github" width="30">
</a>
</p>
''')
