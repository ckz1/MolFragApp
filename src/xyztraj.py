"""`xyztraj` -- 处理 `xyz` 格式的轨迹文件
===========================================

`xyztraj` 主要由 `Traj` , `Mol` 两个对象构成。 `Traj` 在读取 `xyz` 文件后,可获取原子数,
原子列表,以及每帧的标题和坐标,并可进一步计算键长,键角,二面角。 `Traj` 可通过索引得到
对应帧的 `Mol` , `Mol` 可输出 `xyz` 格式的文本内容,也可进行片段划分。`Traj` 侧重于以
不同属性处理整条轨迹的数据,具体来说,整条轨迹的坐标存储为 `[nframe, natom, 3]` 形状的
数组,便于进行矢量计算 , `Mol` 则侧重于通过不同帧处理单个分子的数据。 `Traj` 既支持
通过类方法计算键长,键角,二面角,也支持通过实例方法计算

注意：

* 支持格式: `xyz` , 文件可以只有一帧也可以包含多帧, 包含多帧时每一帧的原子数和原子的排序必须是相同的。

* 读取原子坐标时由于 `float64` 精度有限,可能存在舍入误差

* `Traj` 和 `Mol` 只包含基础的功能其它的功能应该通过继承或者组合添加。

TODO
    COM 质心与质心距离

    距离矩阵与连接矩阵
"""

import os
import numpy as np
import py3Dmol
from sklearn.cluster import DBSCAN
from typing import Any, TypeVar

mol = TypeVar('mol')
"""单分子类变量注解
"""

traj = TypeVar('traj')
"""轨迹类注解
"""


class Mol:
    """单分子类, 表示一个分子或者一帧。可以输出 `xyz` 格式的文本内容； 进行片段分析

    Args:
        natom (int): 原子数 
        atoms (list[str]): 原子列表
        coord (list[list[float]]): 每一帧的原子坐标
        title (str, optional): 每一帧的标题 

    Attributes:
        frag (dict): 片段划分信息,片段索引从0开始。 例如`{0:{'atoms':['H','H','O'],'index':[1,3,5]}}`
        frag_list (list[str]): 排序后的各个片段化学式。

    Example:

        >>> mol = Mol(2, ['H', 'O'], np.array([[1, 2, 3], [4, 5, 6]]), 'molecular OH')
        >>> mol.xyzcontent()
        >>> frag, frag_list = mol.get_frag()
        >>> mol.frag
        >>> mol.gen_chemical_formula(['S', 'C', 'O'])
    """

    def __init__(self,
                 natom: int,
                 atoms: list[str],
                 coord: list[list[float]],
                 title: str = '') -> None:
        self.natom = natom
        self.atoms = atoms
        self.coord = coord
        self.title = title
        self.frag = None
        self.frag_list = None
        # 检查数据
        assert len(self.atoms) == len(
            self.coord
        ), f"原子数目({len(self.atoms)})与原子坐标数目({len(self.coord)})不匹配"

    def xyzcontent(self, savefile=False, filename='output.xyz') -> str:
        """返回 `xyz` 格式的文本内容

        Args:
            savefile (bool, optional): 是否将输出保存为文件. Defaults to False.
            filename (str, optional): 输出文件名. Defaults to 'output.xyz'.

        Returns:
            str: `xyz` 格式的文本内容
        """
        content = f"{self.natom}\n{self.title.rstrip()}\n"
        for i in range(len(self.atoms)):
            content += "{: <4}{:15.8f}{:15.8f}{:15.8f}\n".format(
                self.atoms[i], *self.coord[i])

        if savefile:
            with open(filename, 'w') as f:
                f.write(content)
            f.close()
        return content

    def get_frag(self,
                 max_bond_length: float = 2.5,
                 min_atom_number: int = 1,
                 get_formula=True) -> tuple:
        """根据原子距离使用聚类算法划分片段

        Args:
            max_bond_length (float, optional): 划分不同片段的键长标准. Defaults to 2.5.
            min_atom_number (int, optional): 每个片段的最少的原子数. Defaults to 1.
            get_formula (bool, optional): 是否为输出每个片段的分子式. Defaults to False.

        Note:
            `frag_list` 是按照分子式排序后的列表, 而 `frag` 没有进行排序
        """
        clustering = DBSCAN(eps=max_bond_length,
                            min_samples=min_atom_number).fit(self.coord)
        n_clusters = max(clustering.labels_) + 1
        frag = dict()
        frag_list = []
        # 片段数目
        for i in range(n_clusters):
            frag[i] = {'atoms': [], 'index': []}
        # 按照每个原子的标签将其划分到不同的片段, 索引从1开始与可视化软件一致,使用此索引调用原子坐标时需要减一
        for i, idx_cluster in enumerate(clustering.labels_):
            frag[idx_cluster]['atoms'].append(self.atoms[i])
            frag[idx_cluster]['index'].append(i + 1)
        # 添加片段的化学式
        if get_formula:
            # key 即为片段的索引,从0开始
            for idx_frag in frag:
                frag[idx_frag]['formula'] = Mol.gen_chemical_formula(
                    frag[idx_frag]['atoms'])
                frag_list.append(frag[idx_frag]['formula'])

        frag_list = sorted(frag_list)

        self.frag = frag
        self.frag_list = frag_list

        return frag, frag_list

    def __repr__(self) -> str:
        return f"Mol({self.natom} atoms)"

    @classmethod
    def gen_chemical_formula(cls, atoms: list[str] = ['H', 'H', 'O']) -> str:
        """根据原子列表生成分子式

        Args:
            atoms (list[str], optional): 原子列表. Defaults to ['H', 'H', 'O'].

        Returns:
            str: 分子式, 按照元素的字母序号排序(此处的处理比较粗糙,分子写法有更复杂的规范)
        """
        ele = []
        num = []

        for atom in atoms:
            if atom not in ele:
                ele.append(atom)
                num.append(1)
            else:
                num[ele.index(atom)] += 1

        chemical_formula = []
        for i in range(len(ele)):
            chemical_formula.append(ele[i] +
                                    (str(num[i]) if num[i] > 1 else ''))

        return ''.join(sorted(chemical_formula))


class Traj:
    """轨迹类

    Args:
        trajfile (str): 轨迹文件的路径

    Attributes:
        natom (int): 原子数 `read_trajfile`
        title (list[str]): 每一帧的标题 `read_trajfile`
        title_info (list[str]): 每一帧的标题信息 `parse_title`
        atoms (list[str]): 原子列表 `read_trajfile`
        coord (list[list[float]]): 每一帧的原子坐标, `read_trajfile`

    Example:

        >>> traj = Traj('test_files/traj.xyz')
        >>> traj.natom == 84
        >>> traj.bond(1, 2)
        >>> traj.angle(8, 9, 25)
        >>> traj.dihedral(32, 13, 12, 29)
    """

    def __init__(self, trajfile: str) -> None:
        self.trajfile = trajfile
        assert os.path.exists(
            self.trajfile), f"Can NOT find file {self.trajfile}!"
        # 初始值
        self.natom = None
        self.title = None
        self.title_info = None
        self.atoms = None
        self.coord = None

        # 加载数据
        self.read_trajfile()

    def read_trajfile(self) -> None:
        """从轨迹文件中读取如下信息: 原子数 `natom` , 标题 `title` list(str), 
        标题信息 `title_info` list(obj) 原子列表 原子坐标 化学式 字符串  片段

        Raises:
            ValueError: 轨迹文件读取失败。检查轨迹文件是否满足 `xyz` 格式: 
            1. 第一行是原子数 2. 第二行是标题 3. 第三行到最后一行是元素符号和原子坐标
        """

        # 读取文件
        try:
            with open(self.trajfile, 'r') as f:
                lines = f.readlines()
            f.close()

            # 相同的信息仅读取第一帧
            self.natom = int(lines[0])
            self.atoms = [
                lines[i].split()[0] for i in range(2, 2 + self.natom)
            ]

            # 读取每一帧的标题与原子坐标
            self.title = [
                lines[i] for i in range(1, len(lines), self.natom + 2)
            ]

            coord = []
            for irow in range(2, len(lines), self.natom + 2):
                coord.append(Traj.parse_coord(lines[irow:irow + self.natom]))

            self.coord = np.array(coord, dtype=np.float64)
        except:
            raise IOError(f"Failed to read file {self.trajfile}, please check it!")

    def parse_title(self, parse_func: object = None) -> None:
        """解析xyz文件的标题

        Args:
            parse_func (object): 解析标题内容的函数,要求输入为文本,输出不限

        Example:
            标题为 `i =        5, time =        5.000, E =      -542.5433274126` 解析函数为::

                parse_func = lambda title: [eval(item.split("=")[-1]) for item in title.split(",")]

        Returns:
            object: 具体格式与内容由 `parse_func` 决定
        """
        if parse_func is not None:
            self.title_info = [parse_func(t) for t in self.title]

    def bond(self, atom_idx1: int, atom_idx2: int, index_start: int = 1):
        """计算两个原子之间的键长,原子的索引从1开始,与可视化程序一致。

        Args:
            atom_idx1 (int): 第一个原子的索引
            atom_idx2 (int): 第二个原子的索引
            index_start (int, optional): 原子索引的起始值,从1开始与可视化程序一致,从零开始与 `self.coord` 索引一致,索引起始值与索引应该一致. Defaults to 1.

        Returns:
            list[float]: 键长数组

        Example:

            >>> Traj.bond(1, 2)

            .. image:: ../../src/img/bond.jpg
                :width: 200 px
                :align: center
        """
        return Traj.cal_bond(self.coord, atom_idx1 - index_start,
                             atom_idx2 - index_start)

    def angle(self,
              side_atom_idx1: int,
              vertex_atom_idx: int,
              side_atom_idx2: int,
              index_start: int = 1):
        """计算三个原子之间的键角,原子的索引从1开始,与可视化程序一致。

        Args:
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx (int): 顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引
            index_start (int, optional): 原子索引的起始值. Defaults to 1.

        Returns:
            list[float]: 键角,单位为 `degree`

        Example:

            >>> Traj.angle(3, 1, 2)

            .. image:: ../../src/img/angle.jpg
                :width: 200 px
                :align: center
        """
        return Traj.cal_angle(self.coord, side_atom_idx1 - index_start,
                              vertex_atom_idx - index_start,
                              side_atom_idx2 - index_start)

    def dihedral(self,
                 side_atom_idx1: int,
                 vertex_atom_idx1: int,
                 vertex_atom_idx2: int,
                 side_atom_idx2: int,
                 index_start: int = 1):
        """计算二面角(side_atom_idx1)-(vertex_atom_idx1)-(vertex_atom_idx2)-(side_atom_idx2)。二面角公式为两个平面法向量的夹角。

        Args:
            coord (np.ndarray): 原子坐标 `[nframe,natom,3]`
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx1 (int): 第一个顶点处原子的索引
            vertex_atom_idx2 (int): 第二顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引
            index_start (int): 原子索引的起始值. Defaults to 1.

        Returns:
            np.ndarray: 二面角,形状为 `[natom,1]` , 单位为 `degree`

        Example:

            >>> Traj.dihedral(4, 1, 2, 3)

            .. image:: ../../src/img/dihedral.jpg
                :width: 400 px
                :align: center
        """
        return Traj.cal_dihedral(self.coord, side_atom_idx1 - index_start,
                                 vertex_atom_idx1 - index_start,
                                 vertex_atom_idx2 - index_start,
                                 side_atom_idx2 - index_start)

    def __len__(self):
        return len(self.coord)

    def __getitem__(self, key: int) -> mol:
        frame = Mol(self.natom, self.atoms, self.coord[key], self.title[key])
        return frame

    def __repr__(self) -> str:
        return f"Trajectory with {len(self.coord)} frames({self.natom} atoms) from file {os.path.abspath(self.trajfile)}."

    @classmethod
    def parse_coord(cls, coord_str_list: list) -> np.ndarray:
        """将xyz文件坐标部分的一行转换为原子坐标

        Args:
            coord_str_list (str): xyz文件坐标部分的一行

        Returns:
            np.ndarray: 原子坐标
        """
        coord_arr = np.array([row.split()[1:] for row in coord_str_list],
                             dtype='float64')
        return coord_arr

    @classmethod
    def cal_bond(cls, coord: np.ndarray, atom_idx1: int,
                 atom_idx2: int) -> np.ndarray:
        """计算键长,索引从0开始

        Args:
            coord (np.ndarray): 原子坐标 `[nframe,natom,3]`
            atom_idx1 (int): 第一个原子的索引
            atom_idx2 (int): 第二个原子的索引

        Returns:
            np.ndarray: 键长, 形状为 `[natom,1]` ,单位为 `Angstrom`, 与原子坐标的单位相同
        """
        bond = np.sum((coord[:, atom_idx1] - coord[:, atom_idx2])**2,
                      axis=1)**0.5
        return bond

    @classmethod
    def cal_angle(cls, coord: np.ndarray, side_atom_idx1: int,
                  vertex_atom_idx: int, side_atom_idx2: int) -> float:
        """计算键角(side_atom_idx1)-(vertex_atom_idx)-(side_atom_idx2)。二面角公式为两个键的向量的夹角。原子索引从0开始。

        Args:
            coord (np.ndarray): 原子坐标, 形状为`[nframe,natom,3]`
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx (int): 顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引

        Returns:
            np.ndarray: 键角,形状为 `[natom,1]` ,单位为 `degree`
        """
        side1 = coord[:, side_atom_idx1] - coord[:, vertex_atom_idx]
        side2 = coord[:, side_atom_idx2] - coord[:, vertex_atom_idx]
        angle = np.arccos(
            np.sum(side1 * side2, axis=1) / np.sum(side1**2, axis=1)**0.5 /
            np.sum(side2**2, axis=1)**0.5) * (180 / np.pi)
        return angle

    @classmethod
    def cal_dihedral(cls, coord: np.ndarray, side_atom_idx1: int,
                     vertex_atom_idx1: int, vertex_atom_idx2: int,
                     side_atom_idx2: int) -> float:
        """计算二面角(side_atom_idx1)-(vertex_atom_idx1)-(vertex_atom_idx2)-(side_atom_idx2)。二面角公式为两个平面法向量的夹角。原子索引从0开始。

        Args:
            coord (np.ndarray): 原子坐标 `[nframe,natom,3]`
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx1 (int): 第一个顶点处原子的索引
            vertex_atom_idx2 (int): 第二顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引

        Returns:
            np.ndarray: 二面角,形状为 `[natom,1]` , 单位为 `degree`
        """
        side1 = coord[:, side_atom_idx1] - coord[:, vertex_atom_idx1]
        side2 = coord[:, side_atom_idx2] - coord[:, vertex_atom_idx2]
        vertex = coord[:, vertex_atom_idx1] - coord[:, vertex_atom_idx2]
        # 两个平面的法向量
        normal1 = np.cross(vertex, side1)
        normal2 = np.cross(vertex, side2)
        dihedral = np.arccos(
            np.sum(normal1 * normal2, axis=1) / np.sum(normal1**2, axis=1)**0.5
            / np.sum(normal2**2, axis=1)**0.5) * (180 / np.pi)
        return dihedral

    @classmethod
    def cal_com(cls, atoms: list[str], coord: list[list[float]],
                idxs_1: list[int], idxs_2: list[int]) -> float:
        """_summary_

        Args:
            atoms (list[str]): 所有原子的元素符号
            coord (list[list[float]]): 所有原子的坐标
            idxs_1 (list[int]): 第一个片段的bool索引 
            idxs_2 (list[int]): 第二个片段的bool索引
        """
        assert len(atoms) == len(coord)
        # 计算质心 center of mass
        pass


def display_atoms_by_frag(Traj,
                          idx_frag: int,
                          idx_display: int,
                          width: int = 500,
                          height: int = 500):
    """使用不同颜色标记不同解离片段的原子(在Jupyter Notebook中显示)

    Args:
        Traj (object): `xyztraj.Traj` 对象
        idx_frag (int): 轨迹中划分片段所用帧的索引
        idx_display (int): 轨迹中显示分子结构所有帧的索引
        width (int): `py3Dmol viewer` 宽度,单位 `px`
        height (int): `py3Dmol viewer` 高度,单位 `px`

    Returns:
        object: `py3Dmol viewer`

    Example:
        >>> traj = xyztraj.Traj("output.xyz")
        >>> display_atoms_by_frag(traj, -1, 0)

        .. image:: ../../src/img/frag_label.png
            :width: 400 px
            :align: center

    """
    frag_colors = [
        '#FF0000', '#FFFF00', '#00FF00', '#00FFFF', '#FF00FF', '#FF7F50',
        '#CCCCFF'
    ]

    xyzcontent = Traj[idx_display].xyzcontent()
    coord = Traj[idx_display].coord
    frag, _ = Traj[idx_frag].get_frag()

    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(xyzcontent, 'xyz')

    for idx_frag in frag.keys():
        for idx_atom, atom in zip(frag[idx_frag]['index'],
                                  frag[idx_frag]['atoms']):
            viewer.addLabel(
                f"{atom}{idx_atom}", {
                    'position': {
                        'x': coord[idx_atom - 1][0],
                        'y': coord[idx_atom - 1][1],
                        'z': coord[idx_atom - 1][2]
                    },
                    'backgroundColor': frag_colors[idx_frag],
                    'backgroundOpacity': 0.8
                })

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
    viewer.zoomTo()
    viewer.render()

    return viewer


if __name__ == '__main__':
    pass