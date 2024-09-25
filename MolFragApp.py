import os
import numpy as np
import py3Dmol
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
import plotly.figure_factory as ff
import plotly.graph_objects as go
from glob import glob

st.set_page_config(layout='wide')

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
        if len(temp) == 2:
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


st.title(f'分子解离片段分析(当前路径:`{os.getcwd()}`)')

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
        value = "MOLCAS.template",
    )
    st.write(f"计算级别:`{read_template(template_path)}`")

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
    tab1, tab2 = st.tabs(["轨迹动画&能量曲线", "能量差直方图"])

    with tab1:
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


    with tab2:
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


        
        
components.html('''
<p align='center'>
<a href="https://github.com/ckz1/MolFragApp" target="_blank">
<img border="0" src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png" alt="MolFragApp in Github" width="30">
</a>
</p>
''')