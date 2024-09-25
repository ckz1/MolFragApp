import os
import numpy as np
import py3Dmol
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from glob import glob

st.set_page_config(layout='wide')


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

        num_atom = int(lines[0])
        self.num_atom = num_atom

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

@st.cache_data
def main(xyz_path, max_bond_length: float = 2.5, min_atom_num: int = 1):
    xyzfiles = sorted(glob(xyz_path))
    steps = []
    frags = []
    
    # st.write(max_bond_length,min_atom_num)
    
    progress_text = "Operation in progress. Please wait."
    my_bar = st.progress(0, text=progress_text)
    counter = 0
    num_xyzfile = len(xyzfiles)
    
    for xyzfile in xyzfiles:
        xyz = XYZ(xyzfile)
        xyz.read_xyz()
        steps.append(len(xyz.coords))
        try:
            xyz.get_frag(frame_idx=-1,
                         max_bond_length=max_bond_length,
                         min_atom_num=min_atom_num)
            xyz.get_formula()
            temp = []
            for idx_frag in xyz.frag:
                temp.append(xyz.frag[idx_frag]['formula'])
            frags.append(temp)

        except KeyError:
            print(f'max_bond_length and min_atom_num of {xyzfile} are inappropriate')
            # st.warning(f'max_bond_length and min_atom_num of {xyzfile} are inappropriate')
            frags.append([None])
            
        counter += 1
        my_bar.progress(counter/num_xyzfile, text=progress_text)
    
    my_bar.empty()
    
    # st.write(frags)

    return xyzfiles, steps, frags


st.title('分子解离片段分析')


xyz_path = st.text_input(
    f"当前路径:`{os.getcwd()}`, 轨迹文件(`xyz`格式)路径:",
    value = "Singlet_*/TRAJ_*/output.xyz",
    help = '使用`xyz`文件最后一帧的结构分析解离片段的成分',
)
num_xyzfile = len(glob(xyz_path))
st.write("轨迹文件总数：", num_xyzfile)

if num_xyzfile == 0:
    st.warning('No xyz file found!', icon="⚠️")
    
    
# st.write('**片段划分参数设置:**')

para1, para2 = st.columns(2)


with para1:
    max_bond_length = st.number_input(
        "划分片段时的最大键长(Å)",
        value=2.5,min_value=0.1,
        help='若两原子之间距离大于最大键长，则被划分到不同片段',
    )

with para2:
    min_atom_num = st.number_input(
        "每个片段中原子数的下限",
        value=1, min_value = 1,
        help='取值为1表示划分后的每个片段中至少包含1个原子，即允许存在由单个原子构成的片段',
    )


# st.divider()

st.header('解离片段分析结果')

xyzfiles, steps, frags = main(xyz_path, max_bond_length, min_atom_num)
num_frag = [len(f) for f in frags]
max_num_frag = max(num_frag)

data = {
    "xyz": xyzfiles,
    "steps": steps,
    "nfrag": num_frag,
}
# data = {"xyz file": [os.path.abspath(f) for f in xyzfiles],"steps": steps,}
for idx_frag in range(max_num_frag):
    data[f'frag-{idx_frag+1}'] = []
    for frag_list in frags:
        if idx_frag < len(frag_list):
            data[f'frag-{idx_frag+1}'].append(frag_list[idx_frag])
        else:
            data[f'frag-{idx_frag+1}'].append(None)

st.markdown('可按照某列的值进行排序、搜索；选中某一行观看动画')

# st.dataframe(pd.DataFrame(data), width=800, height=800)

st.session_state.df = pd.DataFrame(data)
    
# event = st.dataframe(
    # st.session_state.df,
    # width=800,
    # key="data",
    # hide_index = True,
    # on_select="rerun",
    # selection_mode="single-row",
# )

# event.selection

# idx_row = event.selection["rows"][0]
# st.session_state.df.iloc[idx_row]['xyz']



# def show_mol(xyzfile,width:int=400,height:int=400):
    # with open(xyzfile,'r') as f:
        # mol_block = f.read()
    # f.close()
    
    # viewer = py3Dmol.view(width=400, height=400)
    # viewer.addModelsAsFrames(mol_block, 'xyz')
    # viewer.setStyle({'stick': {'radius': 0.1,'colorscheme':'Jmol'}, 'sphere': {'radius': 0.3,'colorscheme':'Jmol'}})
    # viewer.setBackgroundColor('black')
    # viewer.animate({'loop': 'forward',})
    # viewer.zoomTo()
    
    # components.html(
        # viewer._make_html(),
        # width=width,
        # height=height
        # )
        
# def show_mol(xyzfile,width:int=400,height:int=400):
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
            names=[col_header.strip() for col_header in title_row[1:].split('|')[:-1]]
            
            data = pd.read_csv(f"{os.path.dirname(xyzfile)}/output_data/expec.out",delimiter='\s+',skiprows=3,names=names)
            
            data_plot = pd.DataFrame(
                {
                    "Time": data.iloc[:,0],
                    "State1": data.iloc[:,4],
                    "State2": data.iloc[:,5],
                    "State3": data.iloc[:,6],
                    "State4": data.iloc[:,7],
                    "Trajectory": data.iloc[:,2],
                }
            )
            
            st.line_chart(data_plot,x="Time",y=["State1","State2","State3","State4","Trajectory"],x_label = 'Time [fs]',y_label= 'Energy [eV]',height=350,width=800)
            
            # st.write(data_plot)


# if len(event.selection["rows"]) != 0:
    # idx_row = event.selection["rows"][0]
    # xyzfile = st.session_state.df.iloc[idx_row]['xyz']
    
    # show_mol(xyzfile,700)
    
    # st.write('缩放: 滚动鼠标滑轮, 旋转: 拖动鼠标左键, 平移: 拖动鼠标滑轮')
    
    
col1, col2 = st.columns(2)

with col1:
    event = st.dataframe(
        st.session_state.df,
        height=700,
        width=800,
        key="data",
        hide_index = True,
        on_select="rerun",
        selection_mode="single-row",
    )
    
    st.write('缩放: 滚动鼠标滑轮, 旋转: 拖动鼠标左键, 平移: 拖动鼠标滑轮')
    
with col2:
    if len(event.selection["rows"]) != 0:
        idx_row = event.selection["rows"][0]
        xyzfile = st.session_state.df.iloc[idx_row]['xyz']
        
        show_mol(xyzfile,width=850,height=350)
        

components.html('''
<p align='center'>
<a href="https://github.com/ckz1/MolFragApp" target="_blank">
<img border="0" src="https://github.githubassets.com/assets/GitHub-Mark-ea2971cee799.png" alt="MolFragApp in Github" width="30">
</a>
</p>
''')