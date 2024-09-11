import os
import numpy as np
import py3Dmol
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from glob import glob


def parse_coord(coord_str_list: list):
    """将`element coord_x coord_y coord_z`的文本转换为一帧的原子坐标数组

    Args:
        coord_str_list (list): 包含元素和原子坐标的字符串数组

    Returns:
        NDarray: 原子坐标数组(一帧)
    """
    coord_arr = np.array([row.split()[1:] for row in coord_str_list],
                         dtype='float64')
    return coord_arr


def read_xyz(xyzfile):
    """读取多帧xyz文件，输出原子坐标数组

    Args:
        xyzfile (str): xyz文件路径

    Returns:
        NDarray: 原子坐标数组
    """
    with open(xyzfile, 'r') as f:
        lines = f.readlines()
    f.close()

    num_atom = int(lines[0])
    coords = []
    for i in range(2, len(lines), num_atom + 2):
        coords.append(parse_coord(lines[i:i + num_atom]))
    return np.array(coords)


def gen_chemical_formula(atom_list: list = ['H', 'H', 'O']):
    """将原子列表转化为化学式

    Args:
        atom_list (list, optional): 原子列表. Defaults to ['H', 'H', 'O'], i.e. H2O.

    Returns:
        str: 化学式
    """
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
    """传入xyz文件路径初始化

    Attribute:
        read_xyz : 读取坐标
        get_frag : 碎片分析
        get_formula : 将碎片结果转换为化学式

    Note:
        碎片分析结果包含原子的元素名以及索引,还有更多信息可以挖掘
    """

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
            print(
                f'max_bond_length and min_atom_num of {xyzfile} are inappropriate'
            )
            # st.warning(f'max_bond_length and min_atom_num of {xyzfile} are inappropriate')
            frags.append([None])

        counter += 1
        my_bar.progress(counter / num_xyzfile, text=progress_text)

    my_bar.empty()

    # st.write(frags)

    return xyzfiles, steps, frags


st.title('分子解离片段分析')

xyz_path = st.text_input(
    f"当前路径:`{os.getcwd()}`, 轨迹文件(`xyz`格式)路径:",
    value="Singlet_*/TRAJ_*/output.xyz",
    help='使用`xyz`文件最后一帧的结构分析解离片段的成分',
)
num_xyzfile = len(glob(xyz_path))
st.write("轨迹文件总数：", num_xyzfile)

if num_xyzfile == 0:
    st.warning('No xyz file found!', icon="⚠️")

# st.write('**片段划分参数设置:**')

max_bond_length = st.number_input(
    "划分片段时的最大键长(Å)",
    value=2.5,
    min_value=0.1,
    help='若两原子之间距离大于最大键长，则被划分到不同片段',
)

min_atom_num = st.number_input(
    "每个片段中原子数的下限",
    value=1,
    min_value=1,
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

event = st.dataframe(
    st.session_state.df,
    width=800,
    key="data",
    hide_index=True,
    on_select="rerun",
    selection_mode="single-row",
)

# st.markdown('没有碎片结果说明最少原子数设置不合适')

# event.selection

# idx_row = event.selection["rows"][0]
# st.session_state.df.iloc[idx_row]['xyz']


def show_mol(xyzfile, width: int = 400, height: int = 400):
    """通过py3dmol显示轨迹动画，并将其嵌入到html

    Args:
        xyzfile (_type_): _description_
        width (int, optional): _description_. Defaults to 400.
        height (int, optional): _description_. Defaults to 400.
    """
    with open(xyzfile, 'r') as f:
        mol_block = f.read()
    f.close()
    viewer = py3Dmol.view(width=800, height=400)
    viewer.addModelsAsFrames(mol_block, 'xyz')
    viewer.setStyle({
        'stick': {
            'radius': 0.1,
            'colorscheme': 'Jmol'
        },
        'sphere': {
            'radius': 0.3,
            'colorscheme': 'Jmol'
        }
    })
    viewer.setBackgroundColor('black')
    viewer.animate({
        'loop': 'forward',
    })
    viewer.zoomTo()

    components.html(viewer._make_html(), width=width, height=height)


if len(event.selection["rows"]) != 0:
    idx_row = event.selection["rows"][0]
    xyzfile = st.session_state.df.iloc[idx_row]['xyz']

    show_mol(xyzfile, 700)

    st.write('缩放: 滚动鼠标滑轮, 旋转: 拖动鼠标左键, 平移: 拖动鼠标滑轮')