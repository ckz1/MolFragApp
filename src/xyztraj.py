"""`xyztraj` -- 处理 `xyz` 格式的轨迹文件
===========================================

`xyztraj` 主要由 `Traj` , `Mol` 两个对象构成。 `Traj` 在读取 `xyz` 文件后,可获取原子数,
原子列表,以及每帧的标题和坐标,并可进一步计算键长,键角,二面角。 `Traj` 可通过索引得到
对应帧的 `Mol` , `Mol` 可输出 `xyz` 格式的文本内容,也可进行片段划分。`Traj` 侧重于以
不同属性处理整条轨迹的数据,具体来说,整条轨迹的坐标存储为 `(nframe, natom, 3)` 形状的
数组,便于进行矢量计算 , `Mol` 则侧重于通过不同帧处理单个分子的数据。 `Traj` 既支持
通过类方法计算键长,键角,二面角,也支持通过实例方法计算

注意：

* 支持格式: `xyz` , 文件可以只有一帧也可以包含多帧, 包含多帧时每一帧的原子数和原子的排序必须是相同的。

* 读取原子坐标时由于 `float64` 精度有限,可能存在舍入误差

* `Traj` 和 `Mol` 只包含基础的功能其它的功能应该通过继承或者组合添加。

TODO
    -[x] COM 质心与质心距离

    -[x] 距离矩阵与连接矩阵

    -[x] 应该也能处理 `xyz` 格式的速度文件,此时 `self.coord` 表示速度

Note:
    原子索引从 `1` 开始,与可视化程序保持一致

Changelog:

- 2025年10月2日
"""

import os
from typing import Any, TypeVar, Type, Self, Callable, Literal

import numpy as np
import py3Dmol
from sklearn.cluster import DBSCAN
from periodictable import elements, constants


# 常数
# Ref: https://docs.scipy.org/doc/scipy-1.16.0/reference/constants.html
FEMTO = 1e-15
"""float: femto 代表的数量级
"""
BOHR_IN_METER = 5.29177210544e-11
"""float: Bohr 半径(原子单位下的距离单位) -> 米
"""
ANGSTROM_IN_METER = 1.00001495e-10
"""float: 埃米 ( `xyz` 文件中的坐标单位) -> 米
"""
BOHR_IN_ANGSTROM = BOHR_IN_METER / ANGSTROM_IN_METER
"""float: Bohr 半径 -> 埃米
"""
AMU_IN_KG = 1.66053906892e-27
"""float: 原子质量单位 `AMU` -> 千克
"""
HARTREE_IN_EV = 27.211386245981
"""float: Hartree (原子单位下的能量单位) -> 电子伏特
"""
EV_IN_JOULE = 1.602176634e-19
"""float: 电子伏特 -> 焦耳
"""
TIME_AU_IN_SECOND = 2.4188843265864e-17
"""float: atomic unit of time in sec
"""
MASS_AU_IN_KG = 9.1093837139e-31
"""float: atomic unit of mass in kg
"""
AMU_IN_AU = AMU_IN_KG / MASS_AU_IN_KG
"""float: 原子质量单位转原子单位
"""

# 类型
form = TypeVar("form", bound="Formula")
"""化学式
"""

mol = TypeVar("mol", bound="Mol")
"""单分子类变量注解
"""

traj = TypeVar("traj", bound="Traj")
"""轨迹类注解
"""


class Formula:
    """化学式。参考 [Chemical formula type — ASE documentation](https://wiki.fysik.dtu.dk/ase/ase/formula.html)

    Example:

        >>> Formula(formula = 'H2O', fmt = '')  # 直接创建

        >>> Formula.from_list(['H', 'H', 'O'], fmt = '')  # 使用方法创建

        >>> W = Formula.from_dict({"H": 2, "O": 1}, fmt='latex')
        >>> W.formula, W.fmt

        >>> Formula('H2O',''), f"{Formula('H2O','')}", str(Formula('H2O',''))
        >>> (<__main__.F at 0x270b3c05e80>, 'H2O', 'H2O')

        >>> print(Formula('H2O',''))
        >>> H2O

    Important:
        化学式的格式 `fmt` 设定后不要改变
    """

    # 用于上下标的特殊符号
    # [下标/上标符号 - ²](https://cn.piliapp.com/symbol/subscript-superscript/)
    # [上标符号 - 特殊符号大全](https://www.shubang.net/fuhao/shangbiao/)
    # [Cool Symbols & Cool Fonts - Symbols, Emoji & Fonts ✮✢❂✶✧](https://coolsymbol.com/)
    # [Unicode 符号表 - 所有 Unicode 字符及其代码都在一页上 (◕‿◕) SYMBL](https://symbl.cc/cn/unicode-table/#malayalam)
    # [Unicode字符表](https://www.rapidtables.org/zh-CN/code/text/unicode-characters.html)
    SPECIAL_SYMBOL: dict[str, dict[str, str]] = {
        "sup": {
            "0": "⁰",
            "1": "¹",
            "2": "²",
            "3": "³",
            "4": "⁴",
            "5": "⁵",
            "6": "⁶",
            "7": "⁷",
            "8": "⁸",
            "9": "⁹",
            "+": "⁺",
            "-": "⁻",
            ".": "˙",
        },
        "sub": {
            "0": "₀",
            "1": "₁",
            "2": "₂",
            "3": "₃",
            "4": "₄",
            "5": "₅",
            "6": "₆",
            "7": "₇",
            "8": "₈",
            "9": "₉",
            "+": "₊",
            "-": "₋",
            ".": "̣",
        },
    }
    # print(''.join(SPECIAL_SYMBOL['sup'].get(s,"␣") for s in str('0123456789.+-')))

    FORMULA_FORMAT = Literal["", "latex", "html", "unicode"]

    def __init__(self, formula: str = "H2O", fmt: FORMULA_FORMAT = ""):
        """初始化化学式

        Args:
            formula (str, optional): 化学式字符串. Defaults to "H2O".
            fmt (str, optional): 化学式格式. Defaults to "". `fmt = ''` 表示无格式, `fmt = 'latex'` 表示使用 LaTeX 格式的上下标; `fmt = 'html'` 表示使用 HTML 格式表示上下标; `fmt = 'unicode'` 表示使用 Unicode 特殊字符表示上下标, 此选项只支持 `FORMULA.SPECIAL_SYMBOL` 中包含的字符映射。
        """
        self.formula = formula
        self.fmt = fmt

    @classmethod
    def hill_order(cls, symbol: str = "H", contain_carbon: bool = False) -> int:
        """获取元素符号的Hill顺序

        Args:
            symbol (str, optional): 元素符号. Defaults to "H".
            contain_carbon (bool, optional): 化学式中是否包含"C". Defaults to False.

        Returns:
            int: Hill顺序

        Note:
            Hill order, 结构式的系统,首先表示分子中的碳原子数,随后是氢原子数,然后是其他所有化学元素的原子数,按元素符号的字母顺序排列。当化学式不含碳时,所有元素（包括氢）按字母顺序列出。使用希尔系统的化学式列表按字母顺序排列如上。当符号以相同字母开头时,单字母元素位于两个字母符号之前(“B”在“Be”之前,“Be”在“B”之前”)。

            参考:

            - [有机物分子式中元素的排列顺序应该是怎么样的？ - 知乎](https://www.zhihu.com/question/312415825?sort=created)
            - [Hill Notation and Hill Order - CHEM.2600: Information Retrieval for Chemists - LibGuides at University of Massachusetts Lowell](https://libguides.uml.edu/c.php?g=110997&p=719419)

            通过 `ord('A') * 1000` 将元素符号映射为数值,每个字符设置三位数字表示,数值为ascii十进制表示

            ASCII 表

            - [ASCII 表 | 菜鸟教程](https://www.runoob.com/w3cnote/ascii.html)
            - [ASCII码对照表,ASCII码一览表（非常详细） - C语言中文网](https://c.biancheng.net/c/ascii/)
        """
        # 包含碳时,C在前,H其次,其它元素按照字母表顺序排列,单字符元素放在双字符元素前
        C_H_map = {"C": 0, "H": 1}
        if len(symbol) == 1:
            if contain_carbon and symbol in ["C", "H"]:
                return C_H_map[symbol]
            # 保证单个字符在双字符之前
            return ord(symbol) * 1000
        else:
            return ord(symbol[0]) * 1000 + ord(symbol[1])

    def add_charge(self, charge: float = 0.0) -> Self:
        """为化学式字符串添加电荷表示

        Args:
            charge (float, optional): 电荷数. Defaults to 0.0.

        Returns:
            Self@Formula: `self` 对象
        """
        charge = round(charge)
        if charge > 0.0:
            signal = "+"
        elif charge < 0.0:
            signal = "-"
        else:
            # 没有电荷,直接跳过
            return self

        # 修改 `formula` 属性, 加上电荷表示
        if self.fmt == "latex":
            self.formula += "^{%s%s}" % (
                str(abs(charge)) if abs(charge) > 1 else "",
                signal,
            )
        elif self.fmt == "unicode":
            self.formula += "%s%s" % (
                (
                    Formula.SPECIAL_SYMBOL["sup"].get(str(abs(charge)), "␣")
                    if abs(charge) > 1
                    else ""
                ),
                Formula.SPECIAL_SYMBOL["sup"].get(signal, "␣"),
            )
        elif self.fmt == "html":
            self.formula += "<sup>%s%s</sup>" % (
                str(abs(charge)) if abs(charge) > 1 else "",
                signal,
            )
        else:
            self.formula += "%s%s" % (
                str(abs(charge)) if abs(charge) > 1 else "",
                signal,
            )

        return self

    def __str__(self) -> str:
        """在 `print` , `format` 等函数中使用

        Returns:
            str: 化学式字符

        Example:

            >>> print(Formula())

            >>> f'{Formula()}'

            >>> str(Formula())
        """
        return self.formula

    @staticmethod
    def from_dict(dct: dict[str, int], fmt: FORMULA_FORMAT = "") -> form:
        """从 `{'元素符号': 原子数目}` 列表中创建化学式

        Args:
            dct (dict[str, int]): `{'元素符号': 原子数目}` 列表
            fmt (str, optional): 化学式格式. Defaults to "".

        Returns:
            Formula: `Formula` 对象

        Note:
            生成化学式尽量不要有多余的字符,例如 `"_{}"`
        """

        symbol_num = sorted(
            dct.items(),
            key=lambda item: Formula.hill_order(
                item[0], contain_carbon=True if "C" in dct.keys() else False
            ),
        )

        if fmt == "latex":
            formula = "".join(
                [
                    f"{item[0]}{'_{%d}'%item[1] if item[1] > 1 else ''}"
                    for item in symbol_num
                ]
            )
            # return ''.join([f"{item[0]}_{{{item[1] if item[1] > 1 else ''}}}" for item in symbol_num])
        elif fmt == "unicode":
            formula = "".join(
                [
                    f"{item[0]}{Formula.SPECIAL_SYMBOL['sub'].get(str(item[1]),'␣') if item[1] > 1 else ''}"
                    for item in symbol_num
                ]
            )
        elif fmt == "html":
            formula = "".join(
                [
                    f"{item[0]}{'<sub>%d</sub>'%item[1] if item[1] > 1 else ''}"
                    for item in symbol_num
                ]
            )
        else:
            formula = "".join(
                [f"{item[0]}{item[1] if item[1] > 1 else ''}" for item in symbol_num]
            )
        return Formula(formula, fmt)

    @staticmethod
    def from_list(symbols: list[str], fmt: FORMULA_FORMAT = "") -> form:
        """从原子符号列表创建化学式

        先将原子符号列表转换为 `{'元素符号': 原子数目}` 进而调用方法处理

        Args:
            symbols (list[str]): 原子符号列表
            fmt (str, optional): 化学式格式. Defaults to "".

        Returns:
            Formula: `Formula` 对象
        """
        dct = {}
        for sym in symbols:
            if sym in dct.keys():
                dct[sym] += 1
            else:
                dct[sym] = 1
        return Formula.from_dict(dct=dct, fmt=fmt)


class TestFormula:
    def test_hill_order(self) -> None:
        assert Formula.from_list(symbols=["H", "H", "O"]).formula == "H2O"
        assert (
            Formula.from_dict(dct={"Br": 1, "Cl": 1, "H": 2, "Si": 1}).formula
            == "BrClH2Si"
        )
        assert Formula.from_dict(dct={"C": 1, "Cl": 4}).formula == "CCl4"
        assert Formula.from_dict(dct={"C": 1, "H": 3, "I": 1}).formula == "CH3I"
        assert Formula.from_dict(dct={"C": 2, "H": 5, "Br": 1}).formula == "C2H5Br"
        assert Formula.from_dict(dct={"H": 2, "O": 4, "S": 1}).formula == "H2O4S"
        # 测试单双字符的排序
        assert Formula.from_dict(dct={"B": 1, "Be": 1, "Br": 1}).formula == "BBeBr"

    def test_from_dict(self) -> None:
        assert (
            Formula.from_dict(dct={"C": 2, "H": 4, "O": 2, "N": 1}, fmt="latex").formula
            == "C_{2}H_{4}NO_{2}"
        )
        assert (
            Formula.from_dict(
                dct={"C": 2, "H": 4, "O": 2, "N": 1}, fmt="unicode"
            ).formula
            == "C₂H₄NO₂"
        )
        assert (
            Formula.from_dict(dct={"C": 2, "H": 4, "O": 2, "N": 1}, fmt="html").formula
            == "C<sub>2</sub>H<sub>4</sub>NO<sub>2</sub>"
        )

    def test_from_list(self) -> None:
        assert (
            Formula.from_list(
                symbols=["N", "O", "O", "C", "C", "H", "H", "H", "H"], fmt="latex"
            ).formula
            == "C_{2}H_{4}NO_{2}"
        )
        assert (
            Formula.from_list(
                symbols=["N", "O", "O", "C", "C", "H", "H", "H", "H"], fmt="unicode"
            ).formula
            == "C₂H₄NO₂"
        )
        assert (
            Formula.from_list(
                symbols=["N", "O", "O", "C", "C", "H", "H", "H", "H"], fmt="html"
            ).formula
            == "C<sub>2</sub>H<sub>4</sub>NO<sub>2</sub>"
        )

    def test_add_charge(self) -> None:
        assert Formula().add_charge(charge=0.0).formula == "H2O"
        assert Formula().add_charge(charge=1.0).formula == "H2O+"
        assert Formula().add_charge(charge=2.0).formula == "H2O2+"
        assert Formula().add_charge(charge=-1.0).formula == "H2O-"
        assert Formula().add_charge(charge=-2.0).formula == "H2O2-"

        assert Formula(fmt="latex").add_charge(charge=0.0).formula == "H2O"
        assert Formula(fmt="latex").add_charge(charge=1.0).formula == "H2O^{+}"
        assert Formula(fmt="latex").add_charge(charge=2.0).formula == "H2O^{2+}"
        assert Formula(fmt="latex").add_charge(charge=-1.0).formula == "H2O^{-}"
        assert Formula(fmt="latex").add_charge(charge=-2.0).formula == "H2O^{2-}"

        assert Formula(fmt="html").add_charge(charge=0.0).formula == "H2O"
        assert Formula(fmt="html").add_charge(charge=1.0).formula == "H2O<sup>+</sup>"
        assert Formula(fmt="html").add_charge(charge=2.0).formula == "H2O<sup>2+</sup>"
        assert Formula(fmt="html").add_charge(charge=-1.0).formula == "H2O<sup>-</sup>"
        assert Formula(fmt="html").add_charge(charge=-2.0).formula == "H2O<sup>2-</sup>"

        assert Formula(fmt="unicode").add_charge(charge=0.0).formula == "H2O"
        assert Formula(fmt="unicode").add_charge(charge=1.0).formula == "H2O⁺"
        assert Formula(fmt="unicode").add_charge(charge=2.0).formula == "H2O²⁺"
        assert Formula(fmt="unicode").add_charge(charge=-1.0).formula == "H2O⁻"
        assert Formula(fmt="unicode").add_charge(charge=-2.0).formula == "H2O²⁻"


class Mol:
    """单分子类, 表示一个分子或者分子动力学中的一帧。可以输出 `xyz` 格式的文本内容, 也可以进行片段分析

    Example:

        >>> mol = Mol(2, ['H', 'O'], np.array([[1, 2, 3], [4, 5, 6]]), 'molecular OH')
        >>> mol.xyz()
        >>> frag, frag_list = mol.frag()
        >>> mol.frag_info
        >>> mol.gen_chemical_formula(['S', 'C', 'O'])
    """

    def __init__(
        self,
        natom: int,
        atoms: list[str],
        coord: list[list[float]] | np.ndarray,
        title: str = "",
    ) -> None:
        """初始化

        Args:
            natom (int): 原子数
            atoms (list[str]): 原子列表(元素符号)
            coord (list[list[float]] | np.ndarray): 原子坐标, 形状 `(natom, 3)` , 单位 `Angstrom`
            title (str, optional): `xyz` 文件的标题. Defaults to "".
        """
        self.natom = natom
        self.atoms = atoms
        self.coord = np.array(coord)
        self.title = title
        self.frag_info = None
        self.frag_list = None
        # 检查原子数natom与原子坐标的数量是否一致
        assert len(self.atoms) == len(
            self.coord
        ), f"原子数目({len(self.atoms)})与原子坐标数目({len(self.coord)})不匹配"

    def xyz(self, filename: None | str = None) -> str:
        """返回 `xyz` 格式的文本内容

        Args:
            filename (None | str, optional): 输出文件名. Defaults to None. 即为不保存文件。

        Returns:
            str: `xyz` 格式的文本内容
        """
        content = f"{self.natom}\n{self.title.rstrip()}\n"
        for i in range(len(self.atoms)):
            content += "{: <4}{:15.8f}{:15.8f}{:15.8f}\n".format(
                self.atoms[i], *self.coord[i]
            )
        # 保存文件
        if filename is not None:
            with open(filename, "w") as f:
                f.write(content)
            f.close()
        return content

    def frag(
        self,
        max_bond_length: float = 2.5,
        min_atom_number: int = 1,
        get_formula: bool = True,
        fmt: Formula.FORMULA_FORMAT = "",
        index_start: int = 1,
    ) -> tuple[dict[int, dict], list[str]]:
        """根据原子距离使用 `DBSCAN` 算法划分片段

        Args:
            max_bond_length (float, optional): 划分不同片段的键长标准. Defaults to 2.5.
            min_atom_number (int, optional): 每个片段所包含的原子数目的最小值. Defaults to 1.
            get_formula (bool, optional): 是否输出每个片段的分子式. Defaults to True.
            fmt (Formula.FORMULA_FORMAT, optional): 分子式的格式. Defaults to "".
            index_start (int, optional): 片段所包含原子的索引的起始值. Defaults to 1, 与可视化软件一致。

        Returns:
            tuple[dict[int, dict], list[str]]: `frag_info` 没有进行排序的片段划分信息, 片段索引从 `0` 开始。 例如 `{0:{'atoms':['H','H','O'],'index':[1,3,5]}}` 。`frag_list` 是按照字母顺序排序后的片段分子式列表

        Important:
            `frag_info` 中片段索引默认从 `1` 开始。
        """
        clustering = DBSCAN(eps=max_bond_length, min_samples=min_atom_number).fit(
            self.coord
        )
        n_clusters = max(clustering.labels_) + 1
        frag_info = dict()
        frag_list = []
        # 片段数目
        for i in range(n_clusters):
            frag_info[i] = {"atoms": [], "index": []}
        # 按照每个原子的标签将其划分到不同的片段, 索引从1开始与可视化软件一致,使用此索引调用原子坐标时需要减一
        for i, idx_cluster in enumerate(clustering.labels_):
            frag_info[idx_cluster]["atoms"].append(self.atoms[i])
            frag_info[idx_cluster]["index"].append(i + index_start)
        # 添加片段的化学式
        if get_formula:
            # key 即为片段的索引,从0开始
            for idx_frag in frag_info:
                # frag[idx_frag]['formula'] = Mol.gen_chemical_formula(frag[idx_frag]['atoms'])
                frag_info[idx_frag]["formula"] = Formula.from_list(
                    frag_info[idx_frag]["atoms"], fmt
                ).formula
                frag_list.append(frag_info[idx_frag]["formula"])

        frag_list = sorted(frag_list)

        self.frag_info = frag_info
        self.frag_list = frag_list

        return frag_info, frag_list

    def __repr__(self) -> str:
        return f"Mol({self.natom} atoms)"

    def distance_matrix(self) -> np.ndarray:
        """距离矩阵

        Returns:
            np.ndarray: 距离矩阵, 形状 `(natom, natom)` , 单位 `Angstrom` , 沿对角形对称, 对角元素为 0.
        """
        return Mol.cal_distance_matrix(self.coord)

    @staticmethod
    def cal_distance_matrix(coord: np.ndarray) -> np.ndarray:
        r"""计算距离矩阵

        Args:
            coord (np.ndarray): 原子坐标矩阵,形状为 `(natom,3)` , 单位 `Angstroms`

        Returns:
            np.ndarray: 距离矩阵,形状 `(natom, natom)` , 单位 `Angstroms`, 索引为 `(i,j)` 的元素表示索引为 `i` 的原子与索引为 `j` 的原子之间的距离

        Tip:
            计算方法: 将原子坐标按照不同方向重复, 再利用**矢量减法**和求模操作获取距离。

            ![](assets/images/repeat.png){ width="200" }

            $$
            \mathbf{R}^{n_f \times 3} \to \mathbf{R}^{1 \times n_f \times 3} \overset{\text{repeat}}{\longrightarrow} \mathbf{R}^{n_f \times n_f \times 3}
            $$

            Gaussian 的输出文件的 `Distance matrix (angstroms):` 部分也包含此结果
        """
        natom, _ = coord.shape
        coord = coord[np.newaxis, :]
        coord1 = coord.repeat(repeats=natom, axis=0)
        coord2 = coord1.transpose((1, 0, 2))
        distance_matrix = np.linalg.norm(coord1 - coord2, axis=2)
        return distance_matrix

    def threshold_matrix(self, scale: float = 1.15) -> np.ndarray:
        """根据原子的共价半径和计算判断是否成键的阈值矩阵

        Args:
            scale (float, optional): 共价半径的缩放系数. Defaults to 1.15.

        Returns:
            np.ndarray: 用于判断是否成键的阈值矩阵
        """
        return Mol.cal_threshold_matrix(self.atoms, scale)

    @staticmethod
    def cal_threshold_matrix(atoms: list[str], scale: float = 1.15) -> np.ndarray:
        """根据原子的共价半径和计算判断是否成键的阈值矩阵

        可视化软件(如Multiwfn)常以共价半径和的115%判断成键,但需注意电子结构的影响。参考:
            - [原子半径](https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page))
            - [Van der Waals半径](https://en.wikipedia.org/wiki/Van_der_Waals_radius)

        Args:
            atoms (list[str]): 原子列表
            scale (float, optional): 共价半径的缩放系数. Defaults to 1.15.

        Returns:
            np.ndarray: 用于判断是否成键的阈值矩阵, 形状为 `(natom,natom)`, 单位 `Angstroms`, 索引为 `(i,j)` 的元素表示索引为 `i` 的原子与索引为 `j` 的原子的共价键半径之和再乘以 `scale`
        """
        # scale = 1.0  # for debug
        natom = len(atoms)
        radius_list = np.array(
            [elements.symbol(atom_symbol).covalent_radius for atom_symbol in atoms]
        )

        radius_matrix1 = radius_list[:, np.newaxis].repeat(repeats=natom, axis=1)
        radius_matrix2 = radius_list[np.newaxis, :].repeat(repeats=natom, axis=0)
        threshold_matrix = (radius_matrix1 + radius_matrix2) * scale

        return threshold_matrix

    def connected_matrix(self, scale: float = 1.15) -> np.ndarray:
        """连接矩阵, 表示成键关系的矩阵

        Returns:
            np.ndarray: 形状为 `(natom,natom)`, 索引为 `(i,j)` 的元素表示索引为 `i` 的原子与索引为 `j` 的原子是否成键(True or False)

        ![](assets/images/connected_matrix.png){ width="200" }
        """
        return self.distance_matrix() < self.threshold_matrix(scale)

    def frag_idx(self, scale: float = 1.15) -> list[list[int]]:
        """根据距离矩阵与阈值矩阵的关系判断分子是否成键

        Args:
            scale (float, optional): 阈值矩阵的缩放系数. Defaults to 1.15.

        Returns:
            list[list[int]]: 各个片段对应的原子索引列表(**从0开始**)
        """

        def find_frag(connected_matrix: np.ndarray, start_index: int = 0) -> list[int]:
            """根据 `connected_matrix` 搜索包含 `start_index` 的整个片段的索引

            Args:
                connected_matrix (np.ndarray): 表示成键关系的矩阵
                start_index (int): 搜索的初始索引

            Returns:
                list[int]: 包含 `start_index` 的整个片段的索引
            """
            start_nodes = [start_index]
            frag_indices = [start_index]
            connected_nodes = []
            while True:
                for i in start_nodes:
                    # 搜索成键的原子索引
                    for j, is_connected in enumerate(connected_matrix[i]):
                        # 判断是否成键
                        if is_connected and j not in frag_indices:
                            # print(j)
                            connected_nodes.append(j)
                frag_indices += connected_nodes
                start_nodes = connected_nodes
                if len(connected_nodes) == 0:
                    break
                connected_nodes = []

            return frag_indices

        connected_matrix = self.connected_matrix(scale)
        natom = len(connected_matrix)
        atom_indices = list(range(natom))
        remain_indices = atom_indices
        fragidx_list = []
        while len(remain_indices) > 0:
            frag_indices = find_frag(connected_matrix, remain_indices[0])
            # print(frag_indices)
            fragidx_list.append(frag_indices)
            # 使用集合运算去掉已经被找到的原子的索引
            remain_indices = list(set(remain_indices).difference(set(frag_indices)))

        return fragidx_list

    def mass(self) -> np.ndarray:
        """原子质量

        Returns:
            np.ndarray: 原子质量, 形状 `(natom,)` , 单位 `AMU` (`AMU` != `AU`)
        """
        return Mol.cal_atomic_mass(atoms=self.atoms)

    @staticmethod
    def cal_atomic_mass(
        atoms: list[str], mass_number: list[int | None] | None = None
    ) -> np.ndarray:
        """根据元素符号生成对应的原子质量列表, 可以通过指定**质量数**,获取同位素的原子质量

        Args:
            atoms (list[str]): 原子列表

        Returns:
            np.ndarray: 原子质量列表, 形状为 `(natom,)` , 单位为原子质量单位/amu/碳原子质量的十二分之一

        Info:
            [Core table - periodictable](https://periodictable.readthedocs.io/en/latest/api/core.html#periodictable.core.Element.mass_units)
        """
        if mass_number is None:
            mass_number = [None for _ in atoms]
        return np.array(
            list(
                (
                    elements.isotope(f"{mass_num}-{symbol}").mass
                    if mass_num is not None
                    else elements.symbol(f"{symbol}").mass
                )
                for symbol, mass_num in zip(atoms, mass_number)
            )
        )

    @staticmethod
    def cal_center_of_mass(
        coord: np.ndarray, mass: np.ndarray, fragidx_list: list[list[int]]
    ) -> list[np.ndarray]:
        r"""计算质心(center of mass, COM)

        $$
        \mathbf{r_c} = \frac{ \sum_i^n m_i \mathbf{r_i} }{ \sum_i^n m_i }
        $$

        Args:
            coord (np.ndarray): 原子坐标, 形状 `(natom, 3)` , 单位 `Angstrom`
            mass (np.ndarray): 原子质量, 形状 `(natom,)` , 单位 `AMU` , 作用等同于权重
            fragidx_list (list[list[int]]): 各个片段的索引构成的列表

        Returns:
            list[np.ndarray]: 质心位置, 单位 `Angstrom` , 与 `coord` 相同。

        由 `C`, `H` 原子构成的片段的质心用三个以 `A` 开头的金属原子表示

        ![](assets/images/com.jpg){ width="200" }
        """
        return [
            (coord[fragidx] * mass[fragidx][:, np.newaxis]).sum(axis=0)
            / mass[fragidx].sum()
            for fragidx in fragidx_list
        ]

    @staticmethod
    def cal_frag_ekin(
        fragidx_list: list[list[int]],
        mass: np.ndarray,
        veloc: np.ndarray,
        factor: None | str = "mw",
    ) -> list[float]:
        r"""计算片段质心动能
        
        $$
        E_{kin-com} = \frac{ \mathbf{p}_c^2 }{ 2 m_c }, \mathbf{v}_c = \frac{ \sum_i^n m_i \mathbf{v_i} }{ \sum_i^n m_i } 
        $$
        
        使用爱因斯坦求和计算动能
        
        $$
        \mathsf{V}_{i,j,k} \mathbf{m}_j \to \mathsf{P}_{i,j,k} , \mathsf{P}^{n_f \times n_a \times 3} = \mathbf{m}^{n_a} \mathsf{V}^{n_f \times n_a \times 3}
        $$
        
        对于 Gaussian 的质权速度 (`AMU^1/2 Bohr / sec`) , **计算动量需要乘以** $\sqrt{m}$ 而不是 $m$ :
        
        $$
        v_{mw} = \sqrt{m} v \cdots \sqrt{\text{AMU}} \frac{\text{Bohr}}{\text{sec}} \\
        p = \sqrt{m} v_{mw} \cdots \text{AMU} \frac{\text{Bohr}}{\text{sec}}  \\
        E_k = \frac{p^2}{2m} = \frac{1}{2}mv^2 \cdots \text{AMU} (\frac{\text{Bohr}}{\text{sec}})^2
        $$
        
        动能可以用质权速度的动量和表示
        
        $$
        (v_{mw})_i = \sqrt{m_i} v_i \to E_k = \frac{1}{2} \sum_i^n m_i v_i^2 = \frac{1}{2} \sum_i^n (v_{mw})_i^2
        $$
        
        在原子单位制 `AU` 下
        
        - 质量 $m_e$
        
        - **长度** $a_0 = \frac{\hbar}{m_e c \alpha}$
        
        - **时间** $\frac{\hbar}{E_h}$
        
        - 速度 $\frac{a_0 E_h}{\hbar}$
        
        - 能量 $E_h = m_e c^2 \alpha^2$
        
        $$
        p = mv \cdots m_e \frac{a_0 E_h}{\hbar} \\
        E_k = \frac{p^2}{2 m} \cdots m_e (\frac{a_0 E_h}{\hbar})^2 = m_e (\frac{E_h}{m_e c \alpha})^2 = \frac{E_h^2}{m_e c^2 \alpha^2} = E_h
        $$

        Args:
            fragidx_list (list[list[int]]): 各个片段的索引, 从 `0` 开始
            mass (list[float]): 原子质量列表,形状 `(natom,)` , 单位 `AMU`
            veloc (np.ndarray): 速度,形状 `(natom,3)` , 单位 `AMU^1/2 Bohr / sec` (Gaussian)
            factor : 将能量转换为 SI 的因子。对于 **质权速度** , `factor = AMU_IN_KG * ANGSTROM_IN_METER**2 / EV_IN_JOULE` ;
                    使用原子单位下的质量和速度计算的动能只需要转换能量单位即可。

        Returns:
            float: 各个片段动能, 单位 `eV`
        """
        if factor == "mw":
            factor = AMU_IN_KG * BOHR_IN_METER**2 / EV_IN_JOULE
            exponent = 1 / 2
        elif factor == "au":
            # 使用原子单位下的质量和速度计算的动能只需要转换能量单位即可
            factor = HARTREE_IN_EV
            exponent = 1.0
        else:
            factor = 1.0
            exponent = 1.0
        momentum_atom = np.einsum("ij,i->ij", veloc, np.power(mass, exponent))
        # print(momentum_atom)
        ekin_frag = [
            np.sum((momentum_atom[fragidx].sum(axis=0)) ** 2)
            / (2 * mass[fragidx].sum())
            * factor
            for fragidx in fragidx_list
        ]
        return ekin_frag

    # @staticmethod
    # def cal_frag_momentum(
    #     fragidx_list: list[list[int]],
    #     mass: np.ndarray,
    #     veloc: np.ndarray,
    #     factor: None | str = "mw",
    # ) -> list[float]:
    #     """计算片段的动量

    #     Args:
    #         fragidx_list (list[list[int]]): 各个片段的索引, 从 `0` 开始
    #         mass (np.ndarray): 原子质量列表,形状 `(natom,)` , 单位 `AMU`
    #         veloc (np.ndarray): 速度,形状 `(natom,3)` , 单位 `AMU^1/2 Bohr / sec` (Gaussian)
    #         factor (None | str, optional): 将动量转换为 SI 的因子。. Defaults to "mw".

    #     Returns:
    #         list[float]: 各个片段动量

    #     Note:
    #         尚未测试
    #     """
    #     if factor == "mw":
    #         factor = AMU_IN_KG * BOHR_IN_METER**2 / EV_IN_JOULE
    #         exponent = 1 / 2
    #     elif factor == "au":
    #         # 使用原子单位下的质量和速度计算的动能只需要转换能量单位即可
    #         factor = HARTREE_IN_EV
    #         exponent = 1.0
    #     else:
    #         factor = 1.0
    #         exponent = 1.0
    #     momentum_atom = np.einsum("ij,i->ij", veloc, np.power(mass, exponent))
    #     momentum_frag = [
    #         momentum_atom[fragidx].sum(axis=0) * factor for fragidx in fragidx_list
    #     ]
    #     return momentum_frag

    @staticmethod
    def gen_chemical_formula(atoms: list[str] = ["H", "H", "O"]) -> str:
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
            chemical_formula.append(ele[i] + (str(num[i]) if num[i] > 1 else ""))

        return "".join(sorted(chemical_formula))


class TestMol:
    tmol = Mol(
        9,
        atoms=["H", "C", "C", "H", "C", "C", "H", "H", "C"],
        coord=[
            [0.026491000, -1.700245000, -0.003311000],
            [-1.033701000, -1.813950000, -0.001342000],
            [-2.223402000, -1.926584000, 0.000657000],
            [-1.482872000, 0.824532000, 0.000544000],
            [-1.053520000, 1.800542000, 0.000177000],
            [-0.558762000, 2.888348000, -0.000139000],
            [-0.128583000, 3.860262000, -0.000833000],
            [1.457087000, 0.872108000, 0.000116000],
            [2.087408000, 0.012086000, 0.000263000],
        ],
        title="JUST FOR TEST, NOT REAL MOL",
    )

    def test_xyz(self):
        tmol = Mol(2, ["H", "O"], np.array([[1, 2, 3], [4, 5, 6]]), "molecular OH")
        xyz = tmol.xyz()
        assert len(xyz.strip().split("\n")) == 4
        assert "molecular OH" in xyz
        assert "O" in xyz
        assert "H" in xyz

    def test_frag(self):
        frag_info, frag_list = self.tmol.frag(1.2)
        assert frag_info == {
            0: {"atoms": ["H", "C", "C"], "index": [1, 2, 3], "formula": "C2H"},
            1: {
                "atoms": ["H", "C", "C", "H"],
                "index": [4, 5, 6, 7],
                "formula": "C2H2",
            },
            2: {"atoms": ["H", "C"], "index": [8, 9], "formula": "CH"},
        }
        assert frag_list == ["C2H", "C2H2", "CH"]
        # print(mol.xyz())

    def test_distance_matrix(self):
        dm = self.tmol.distance_matrix()

        assert np.isclose(
            dm[0],
            [
                0.0,
                1.06627378,
                2.26125266,
                2.94154558,
                3.66359735,
                4.62576662,
                5.56266952,
                2.94340222,
                2.67945314,
            ],
        ).all()
        assert np.isclose(
            dm[-1],
            [
                2.67945314,
                3.61603809,
                4.72668217,
                3.66155266,
                3.61441608,
                3.9083371,
                4.44061661,
                1.06627502,
                0.0,
            ],
        ).all()

    def test_cal_distance_matrix(self):
        dm = self.tmol.cal_distance_matrix(self.tmol.coord)

        assert np.isclose(
            dm[0],
            [
                0.0,
                1.06627378,
                2.26125266,
                2.94154558,
                3.66359735,
                4.62576662,
                5.56266952,
                2.94340222,
                2.67945314,
            ],
        ).all()
        assert np.isclose(
            dm[-1],
            [
                2.67945314,
                3.61603809,
                4.72668217,
                3.66155266,
                3.61441608,
                3.9083371,
                4.44061661,
                1.06627502,
                0.0,
            ],
        ).all()

    def test_threshold_matrix(self):
        tm = self.tmol.threshold_matrix()

        assert np.isclose(
            [tm[0, 5], tm[2, 4], tm[5, 7]], [1.2305, 1.7479999999999998, 1.2305]
        ).all()

    def test_cal_threshold_matrix(self):
        tm = self.tmol.cal_threshold_matrix(self.tmol.atoms)
        print(tm[0, 5], tm[2, 4], tm[5, 7])

        assert np.isclose(
            [tm[0, 5], tm[2, 4], tm[5, 7]], [1.2305, 1.7479999999999998, 1.2305]
        ).all()

    def test_connected_matrix(self):
        cm = self.tmol.connected_matrix()

        assert np.isclose(
            [cm[0, 5], cm[2, 4], cm[5, 7], cm[3, 4], cm[7, 8]],
            [False, False, False, True, True],
        ).all()

    def test_frag_idx(self):
        assert self.tmol.frag_idx() == [[0, 1, 2], [3, 4, 5, 6], [8, 7]]

    def test_mass(self):
        assert np.isclose(
            self.tmol.mass(),
            [1.008, 12.011, 12.011, 1.008, 12.011, 12.011, 1.008, 1.008, 12.011],
        ).all()

    def test_cal_atomic_mass(self):
        assert np.isclose(
            self.tmol.cal_atomic_mass(self.tmol.atoms),
            [1.008, 12.011, 12.011, 1.008, 12.011, 12.011, 1.008, 1.008, 12.011],
        ).all()

    def test_cal_center_of_mass(self):
        com = self.tmol.cal_center_of_mass(
            self.tmol.coord, self.tmol.mass(), [[0, 1, 2], [3, 4, 5, 6], [8, 7]]
        )

        assert np.isclose(
            com[0], [-1.56190017e00, -1.86341993e00, -4.62046464e-04]
        ).all()
        assert np.isclose(
            com[1], [-8.06108985e-01, 2.34428643e00, 6.34096321e-06]
        ).all()
        assert np.isclose(com[2], [2.03860521e00, 7.86734626e-02, 2.51618481e-04]).all()

    def test_cal_frag_ekin(self):
        pass
        # self.tmol.cal_frag_ekin()

    def test_gen_chemical_formula(self):
        assert Mol.gen_chemical_formula() == "H2O"


# TestMol().test_cal_center_of_mass()


class Traj:
    """轨迹类

    Example:

        >>> traj = Traj('test_files/traj.xyz')
        >>> traj.natom == 84
        >>> traj.bond(1, 2)
        >>> traj.angle(8, 9, 25)
        >>> traj.dihedral(32, 13, 12, 29)
    """

    def __init__(
        self,
        natom: int = 0,
        title: list[str] = [""],
        title_info: list = [],
        atoms: list[str] = [""],
        coord: np.ndarray = np.array([]),
    ) -> None:
        """初始化

        Args:
            natom (int, optional): 原子数目. Defaults to 0.
            title (list[str], optional): 每一帧的标题. Defaults to [""].
            title_info (list, optional): 从每一帧的标题提取出来的信息. Defaults to [].
            atoms (list[str], optional): 原子的元素符号列表. Defaults to [""].
            coord (np.ndarray, optional): 每一帧的原子坐标. Defaults to np.array([]).
        """
        self.natom = natom
        self.title = title
        self.title_info = title_info
        self.atoms = atoms
        self.coord = coord

    @staticmethod
    def from_xyz(trajfile: str) -> traj:
        """从轨迹文件中读取如下信息: 原子数 `natom` , 原子列表 `atoms` : `list[str]` , 标题 `title` : `list(str)`,
        原子坐标 `coord` : `np.ndarray`

        Args:
            trajfile (str): `xyz` 格式的轨迹文件

        Raises:
            ValueError: 轨迹文件读取失败。检查轨迹文件是否满足 `xyz` 格式:
                1. 第一行是原子数
                2. 第二行是标题
                3. 第三行到最后一行是元素符号和原子坐标

        Returns:
            Traj: `Traj` 对象
        """
        assert os.path.exists(trajfile), f"Can NOT find file {trajfile}!"
        # 读取文件
        try:
            with open(trajfile, "r", encoding="utf-8") as f:
                lines = f.readlines()
            f.close()

            # 相同的信息仅读取第一帧
            natom = int(lines[0])
            atoms = [lines[i].split()[0] for i in range(2, 2 + natom)]

            # 读取每一帧的标题与原子坐标
            title = [lines[i] for i in range(1, len(lines), natom + 2)]

            coord = []
            for irow in range(2, len(lines), natom + 2):
                coord.append(Traj.parse_coord(lines[irow : irow + natom]))

            coord = np.array(coord, dtype=np.float64)

            return Traj(natom=natom, title=title, atoms=atoms, coord=coord)
        except IOError as e:
            raise ValueError(f"Failed to read file {trajfile}, please check it!") from e

    def parse_title(self, parse_func: Callable | None = None) -> None:
        """解析 `xyz` 文件的标题

        Args:
            parse_func (object): 解析标题内容的函数, 要求输入为文本, 输出不限

        Example:
            标题为 `i =        5, time =        5.000, E =      -542.5433274126` 解析函数为:

                parse_func = lambda title: [eval(item.split("=")[-1]) for item in title.split(",")]

        Tip:
            `self.title_info` 具体格式与内容由 `parse_func` 决定
        """
        if parse_func is not None:
            self.title_info = [parse_func(t) for t in self.title]

    def xyz(self, filename: None | str = None) -> str:
        """返回 `xyz` 格式的文本内容

        Args:
            filename (None | str, optional): 输出文件名. Defaults to None. 即为不保存文件。

        Returns:
            str: `xyz` 格式的文本内容
        """
        # 调用 `Mol` 的方法生成 `xyz` 格式的文本文件。没有考虑每帧的标题
        content = "".join([Mol(self.natom, self.atoms, c).xyz() for c in self.coord])
        # 保存文件
        if filename is not None:
            with open(filename, "w") as f:
                f.write(content)
            f.close()
        return content

    def bond(self, atom_idx1: int, atom_idx2: int, index_start: int = 1) -> np.ndarray:
        """计算两个原子之间的键长,原子的索引从1开始,与可视化程序一致。

        Args:
            atom_idx1 (int): 第一个原子的索引
            atom_idx2 (int): 第二个原子的索引
            index_start (int, optional): 原子索引的起始值,从1开始与可视化程序一致,从零开始与 `self.coord` 索引一致,索引起始值与索引应该一致. Defaults to 1.

        Returns:
            np.ndarray: 键长数组

        Example:

            >>> Traj.bond(1, 2)

        ![](assets/images/bond.jpg){ width="200" }
        """
        return Traj.cal_bond(
            self.coord, atom_idx1 - index_start, atom_idx2 - index_start
        )

    def angle(
        self,
        side_atom_idx1: int,
        vertex_atom_idx: int,
        side_atom_idx2: int,
        index_start: int = 1,
    ) -> np.ndarray:
        """计算三个原子之间的键角,原子的索引从1开始,与可视化程序一致。

        Args:
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx (int): 顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引
            index_start (int, optional): 原子索引的起始值. Defaults to 1.

        Returns:
            np.ndarray: 键角, 单位为 `degree`

        Example:

            >>> Traj.angle(3, 1, 2)

        ![](assets/images/angle.jpg){ width="200" }
        """
        return Traj.cal_angle(
            self.coord,
            side_atom_idx1 - index_start,
            vertex_atom_idx - index_start,
            side_atom_idx2 - index_start,
        )

    def dihedral(
        self,
        side_atom_idx1: int,
        vertex_atom_idx1: int,
        vertex_atom_idx2: int,
        side_atom_idx2: int,
        index_start: int = 1,
    ) -> np.ndarray:
        """计算二面角(side_atom_idx1)-(vertex_atom_idx1)-(vertex_atom_idx2)-(side_atom_idx2)。二面角公式为两个平面法向量的夹角。

        Args:
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx1 (int): 第一个顶点处原子的索引
            vertex_atom_idx2 (int): 第二顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引
            index_start (int): 原子索引的起始值. Defaults to 1.

        Returns:
            np.ndarray: 二面角,形状为 `(natom,1)` , 单位为 `degree`

        Example:

            >>> Traj.dihedral(4, 1, 2, 3)

        ![](assets/images/dihedral.jpg){ width="200" }
        """
        return Traj.cal_dihedral(
            self.coord,
            side_atom_idx1 - index_start,
            vertex_atom_idx1 - index_start,
            vertex_atom_idx2 - index_start,
            side_atom_idx2 - index_start,
        )

    def distance_matrix(self):
        """计算距离矩阵"""
        return Traj.cal_distance_matrix(self.coord)

    def __len__(self):
        return len(self.coord)

    def __getitem__(self, key: int) -> Mol:
        # 如果不能索引标题就不考虑标题
        if key < len(self.title):
            frame = Mol(self.natom, self.atoms, self.coord[key], self.title[key])
        else:
            frame = Mol(self.natom, self.atoms, self.coord[key])
        return frame

    def __repr__(self) -> str:
        # return f"Trajectory with {len(self.coord)} frames({self.natom} atoms) from file {os.path.abspath(self.trajfile)}."
        return f"Trajectory ({len(self.coord)} frames, {self.natom} atoms)"

    @staticmethod
    def to_idx(atom_idxs: list[int] = [1]) -> list[int] | int:
        """将从 `1` 开始的索引转换为从 `1` 开始的索引

        Args:
            atom_idxs (list[int], optional): 从 `1` 开始的索引. Defaults to [1].

        Returns:
            list[int] | int: 从 `1` 开始的索引
        """
        if isinstance(atom_idxs, int):
            return atom_idxs - 1
        else:
            return list(i for i in atom_idxs)

    @staticmethod
    def parse_coord(coord_str_list: list) -> np.ndarray:
        """将 `xyz` 文件坐标部分的一行(文本)转换为原子坐标

        Args:
            coord_str_list (str): xyz文件坐标部分的一行

        Returns:
            np.ndarray: 原子坐标
        """
        coord_arr = np.array(
            [row.split()[1:] for row in coord_str_list], dtype="float64"
        )
        return coord_arr

    @staticmethod
    def cal_bond(coord: np.ndarray, atom_idx1: int, atom_idx2: int) -> np.ndarray:
        r"""计算键长, 索引从 `0` 开始

        $$
        r_{12} = \sqrt{ (x_1 - x_2)^2 + (y_1 - y_2)^2 + (z_1 - z_2)^2 }
        $$

        Args:
            coord (np.ndarray): 原子坐标 `(nframe,natom,3)`
            atom_idx1 (int): 第一个原子的索引
            atom_idx2 (int): 第二个原子的索引

        Returns:
            np.ndarray: 键长, 形状为 `(natom,1)` ,单位为 `Angstrom`, 与原子坐标的单位相同
        """
        bond = np.sum((coord[:, atom_idx1] - coord[:, atom_idx2]) ** 2, axis=1) ** 0.5
        return bond

    @staticmethod
    def cal_angle(
        coord: np.ndarray,
        side_atom_idx1: int,
        vertex_atom_idx: int,
        side_atom_idx2: int,
    ) -> np.ndarray:
        r"""计算键角(side_atom_idx1)-(vertex_atom_idx)-(side_atom_idx2)。二面角公式为两个键的向量的夹角。原子索引从 `0` 开始。

        $$
        \theta = \arccos ( \frac{ \mathbf{v}_{v s_1} \cdot \mathbf{v}_{v s_2} }{ |\mathbf{v}_{v s_1}| \cdot |\mathbf{v}_{v s_2}| } )
        $$

        Args:
            coord (np.ndarray): 原子坐标, 形状为`(nframe,natom,3)`
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx (int): 顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引

        Returns:
            np.ndarray: 键角,形状为 `(natom,1)` ,单位为 `degree`
        """
        side1 = coord[:, side_atom_idx1] - coord[:, vertex_atom_idx]
        side2 = coord[:, side_atom_idx2] - coord[:, vertex_atom_idx]
        angle = np.arccos(
            np.sum(side1 * side2, axis=1)
            / np.sum(side1**2, axis=1) ** 0.5
            / np.sum(side2**2, axis=1) ** 0.5
        ) * (180 / np.pi)
        return angle

    @staticmethod
    def cal_dihedral(
        coord: np.ndarray,
        side_atom_idx1: int,
        vertex_atom_idx1: int,
        vertex_atom_idx2: int,
        side_atom_idx2: int,
    ) -> np.ndarray:
        r"""计算二面角(side_atom_idx1)-(vertex_atom_idx1)-(vertex_atom_idx2)-(side_atom_idx2)。二面角公式为两个平面法向量的夹角。原子索引从0开始。
        
        $$
        \mathbf{w}_1 = \mathbf{v}_{v_1 v_2} \times \mathbf{v}_{s_1 v_1} , \mathbf{w}_2 = \mathbf{v}_{v_1 v_2} \times \mathbf{v}_{s_2 v_2} \\
        \theta = \arccos ( \frac{ \mathbf{w}_1 \cdot \mathbf{w}_2 }{ |\mathbf{w}_1| \cdot |\mathbf{w}_2| } )
        $$

        Args:
            coord (np.ndarray): 原子坐标 `(nframe,natom,3)`
            side_atom_idx1 (int): 第一条边处原子的索引
            vertex_atom_idx1 (int): 第一个顶点处原子的索引
            vertex_atom_idx2 (int): 第二顶点处原子的索引
            side_atom_idx2 (int): 第二条边处原子的索引

        Returns:
            np.ndarray: 二面角,形状为 `(natom,1)` , 单位为 `degree`
            
        Info:
            Docs ref: [8.12  Internal coordinates definitions](https://sharc-md.org/?page_id=1454)
        """
        side1 = coord[:, side_atom_idx1] - coord[:, vertex_atom_idx1]
        side2 = coord[:, side_atom_idx2] - coord[:, vertex_atom_idx2]
        vertex = coord[:, vertex_atom_idx1] - coord[:, vertex_atom_idx2]
        # 两个平面的法向量
        normal1 = np.cross(vertex, side1)
        normal2 = np.cross(vertex, side2)
        dihedral = np.arccos(
            np.sum(normal1 * normal2, axis=1)
            / np.sum(normal1**2, axis=1) ** 0.5
            / np.sum(normal2**2, axis=1) ** 0.5
        ) * (180 / np.pi)
        return dihedral

    @staticmethod
    def cal_distance_matrix(coord: np.ndarray) -> np.ndarray:
        r"""计算键长矩阵/距离矩阵(distance matrix)

        $$
        M^{\text{bond}}_{ij} = | \vec{R}_i - \vec{R}_j |
        $$

        Args:
            coord (np.ndarray): 原子坐标 `(nframe,natom,3)`

        Returns:
            np.ndarray: 距离矩阵, 形状为 `(nframe,natom,natom)` , 单位为 `Angstrom`
        """
        nframe, natom, _ = coord.shape
        bm = np.zeros((nframe, natom, natom))
        for iframe in range(len(nframe)):
            for irow in range(natom):
                for icol in range(natom):
                    bm[iframe, irow, icol] = np.linalg.norm(coord[irow] - coord[icol])

        return bm

    @staticmethod
    def cal_com(
        atoms: list[str], coord: list[list[float]], idxs_1: list[int], idxs_2: list[int]
    ) -> Any:
        """计算质心距离(center of mass)

        Args:
            atoms (list[str]): 所有原子的元素符号
            coord (list[list[float]]): 所有原子的坐标
            idxs_1 (list[int]): 第一个片段的bool索引
            idxs_2 (list[int]): 第二个片段的bool索引

        Raises:
            ImportError: 未实现

        Returns:
            Any: 未实现
        """

        assert len(atoms) == len(coord)
        # 计算质心 center of mass
        raise ImportError


# For Jupyter
# class Viewer:
# pass


def display_atoms_by_frag(
    Traj, idx_frag: int, idx_display: int, width: int = 500, height: int = 500
):
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

    ![](assets/images/frag_label.png){ width="400" }
    """
    frag_colors = [
        "#FF0000",
        "#FFFF00",
        "#00FF00",
        "#00FFFF",
        "#FF00FF",
        "#FF7F50",
        "#CCCCFF",
    ]

    xyzcontent = Traj[idx_display].xyz()
    coord = Traj[idx_display].coord
    frag, _ = Traj[idx_frag].frag()

    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(xyzcontent, "xyz")

    for idx_frag in frag.keys():
        for idx_atom, atom in zip(frag[idx_frag]["index"], frag[idx_frag]["atoms"]):
            viewer.addLabel(
                f"{atom}{idx_atom}",
                {
                    "position": {
                        "x": coord[idx_atom - 1][0],
                        "y": coord[idx_atom - 1][1],
                        "z": coord[idx_atom - 1][2],
                    },
                    "backgroundColor": frag_colors[idx_frag],
                    "backgroundOpacity": 0.8,
                },
            )

    viewer.setStyle(
        {
            "stick": {"radius": 0.1, "colorscheme": "Jmol"},
            "sphere": {"radius": 0.3, "colorscheme": "Jmol"},
        }
    )
    viewer.setBackgroundColor("black")
    viewer.zoomTo()
    viewer.render()

    return viewer


if __name__ == "__main__":
    pass
