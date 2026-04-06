"""# SHARC 输出文件的数据提取

提取出来的数据将分两类存储在 HDF5 文件中，具体存储方式为:
- 第一类的原始数据是一个字典(普通的键值对)，存储在 `group` 的 `attrs`
- 第二类的原始数据也是一个字典(数据数据)，存储在 `group` 的 `dataset`
- `single_trajdata` 决定数据的具体内容
- 通过 `Dataset` 的不同方法读取不同的数据
"""

import os
import logging
import warnings
import multiprocessing as mp
from glob import glob
from typing import Any
from pprint import pprint as pp

import psutil
import numpy as np
import pandas as pd
import h5py

import xyztraj


def config_logging(log_file: str = "log") -> None:
    """Config logging"""
    logging.basicConfig(
        level=logging.DEBUG,
        filename=log_file,
        filemode="a",  ##模式，有w和a，w就是写模式，每次都会重新写日志，覆盖之前的日志
        # a是追加模式，默认如果不写的话，就是追加模式
        format="%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s",
        # 日志格式
    )


logger = logging.getLogger()


def sharc_molden(
    vibration_info: dict, geom: list, molden_file: str = "freq.molden"
) -> None:
    """生成 SHARC 所需的包含频率计算结果的 MOLDEN 文件。

    Args:
        vibration_info (dict): 频率计算的信息。
        geom (list): 结构信息。
        molden_file (str, optional): 保存的 MOLDEN 文件名. Defaults to "freq.molden".

    Note:
        [The Molden Format](https://www.theochem.ru.nl/molden/molden_format.html)
    """
    molden_content = {}
    molden_content["n_freq"] = len(vibration_info["vibration_number"])

    molden_content["freq"] = "\n".join(
        list(f"{n:>16.4f}" for n in vibration_info["vibration_frequency"])
    )
    #
    molden_content["int"] = "\n".join(
        [f"{n:12.4f}" for n in vibration_info["ir_intensity"]]
    )

    molden_content["natom"] = vibration_info["num_atom"]
    molden_content["fr_coord"] = "\n".join(
        ["{:<2s}{:>12.6f}{:>12.6f}{:>12.6f}".format(*line) for line in geom]
    )

    fr_norm_coord = []
    for i, arr in enumerate(vibration_info["atomic_displacement"]):
        fr_norm_coord.append(f"vibration{i+1:16d}")
        for line in arr:
            fr_norm_coord.append("%12.6f\t%12.6f\t%12.6f" % tuple(line))
    molden_content["fr-norm-coord"] = "\n".join(fr_norm_coord)

    molden_content["rmass"] = "\n".join(
        list(f"{n:>16.12f}" for n in vibration_info["reduced_mass"])
    )

    # print(molden_content)
    with open(molden_file, "w", encoding="utf-8") as f:
        f.write(
            "\n".join(
                ["[Molden Format]"]
                + [
                    f"[{section_name.upper()}]\n{section_content}"
                    for section_name, section_content in molden_content.items()
                ]
            )
        )


def test_SHARC_molden():
    # SHARC_molden(vibration("tests/GAUSSIAN.opt_freq.log"))

    # print(vibration("tests/GAUSSIAN.opt_freq.log"))

    pass


def valid_data():
    r"""验证输出文件 `output.dat` , `output.xyz` 中数据的单位

    Dalitz Plot要对动能归一化,Newton Plot以第一个片段的动量为单位,
    因此,Dalitz Plot和Newton Plot的动能,动量数据都不用考虑单位,只要保证单位一致,
    使得相对大小不变即可。**但是在计算初末动能差值时必需考虑单位转化**。

    * `output.dat` 文件中使用的是原子单位( `a.u.` ), 长度 `Bohr` 和时间 `\hbar/E_h` 都是,具体可参考
    [原子单位制](https://zh.wikipedia.org/wiki/%E5%8E%9F%E5%AD%90%E5%8D%95%E4%BD%8D%E5%88%B6),例如在 `301 step` 的数据

        .. code-block::

            ! Atomic masses
             0.2187466181995E+005
             0.2915694637199E+005
             0.2915694637199E+005

            ! 7 Ekin (a.u.)
             0.4749769065556E+000

            ! 11 Geometry in a.u.
             0.7518592163037E+000 -0.2955786806237E+001 -0.1631350425357E+002
            -0.8152532594626E+000  0.4607814983189E+001  0.2256129130491E+002
             0.5996443792277E-001 -0.3139576909575E+001 -0.1404514387337E+002

            ! 12 veloc_datities in a.u.
            -0.1733564413781E-003 -0.1038086865396E-002 -0.3399721116931E-002
            -0.1428907104805E-003  0.8655511453946E-003  0.4355393676912E-002
             0.2729573069188E-003 -0.8674541507362E-004 -0.1804789933349E-002

    * `output.lis` 文件中用的是常用单位(表头标记)。与上面 `output.dat` 文件 `301 step` 的对应的动能数据为 :code:`Step 301  Energy kin [eV] 12.924779`

    * `output.xyz` 文件中坐标单位为埃, 时间单位为飞秒(如果存储的是速度数据)。最后两帧(300 step, 301 step)的数据如下(时间间隔间隔为0.5 fs)

        .. code-block::

              t=       150.00000    3    2
            C       0.399704831     -1.552796322     -8.595360108
            O      -0.429850617      2.428883893     11.891284945
            O       0.028789909     -1.660432414     -7.412773666
                       3
              t=       150.50000    3    2
            C       0.397866763     -1.564135018     -8.632734683
            O      -0.431413446      2.438350682     11.938921209
            O       0.031731814     -1.661392553     -7.432370063


    两处的数据是一致的,具体验证过程如下,使用了 `scipy.constants` 中的 :code:`Hartree energy in eV = 27.211386245981 eV`
    """
    # 验证动能
    am_dat = np.array(
        [
            0.2187466181995e005,
            0.2915694637199e005,
            0.2915694637199e005,
        ]
    )
    print("原子质量(AU) <- output.dat\n", am_dat)

    veloc_dat = np.array(
        [
            [-0.1733564413781e-003, -0.1038086865396e-002, -0.3399721116931e-002],
            [-0.1428907104805e-003, 0.8655511453946e-003, 0.4355393676912e-002],
            [0.2729573069188e-003, -0.8674541507362e-004, -0.1804789933349e-002],
        ]
    )
    print("原子速度(AU) <- output.dat\n", veloc_dat)

    ekin_total = 0
    for mi, vi in zip(am_dat, veloc_dat):
        ekin_total += np.sum(0.5 * mi * vi**2)
    hartree2ev = 27.211386245981
    print("原子动能和(AU,eV) <- Ek = mv^2/2\n", ekin_total, ekin_total * hartree2ev)

    # 与output.dat和output.lis数据对比
    print(
        "原子动能和(AU <- output.dat, eV <- output.lis)\n",
        0.4749769065556e000,
        12.924779,
    )
    print(ekin_total - 0.4749769065556e000)
    print(ekin_total * hartree2ev - 12.924779)
    # 计算结果吻合

    # 验证坐标与速度
    geom_last1_xyz = np.array(
        [
            [0.397866763, -1.564135018, -8.632734683],
            [-0.431413446, 2.438350682, 11.938921209],
            [0.031731814, -1.661392553, -7.432370063],
        ]
    )

    geom_last2_xyz = np.array(
        [
            [0.399704831, -1.552796322, -8.595360108],
            [-0.429850617, 2.428883893, 11.891284945],
            [0.028789909, -1.660432414, -7.412773666],
        ]
    )

    geom_last1_dat = np.array(
        [
            [0.7518592163037e000, -0.2955786806237e001, -0.1631350425357e002],
            [-0.8152532594626e000, 0.4607814983189e001, 0.2256129130491e002],
            [0.5996443792277e-001, -0.3139576909575e001, -0.1404514387337e002],
        ]
    )

    veloc_last1_dat = np.array(
        [
            [-0.1733564413781e-003, -0.1038086865396e-002, -0.3399721116931e-002],
            [-0.1428907104805e-003, 0.8655511453946e-003, 0.4355393676912e-002],
            [0.2729573069188e-003, -0.8674541507362e-004, -0.1804789933349e-002],
        ]
    )

    # 埃米转玻尔半径,与output.dat数据一致
    print(
        "xyz文件中原子坐标(-> Bohr)及其与dat文件中差值\n",
        geom_last1_xyz / 0.5291772083,
        "\n",
        geom_last1_xyz / 0.5291772083 - geom_last1_dat,
    )

    # xyz的坐标加上模拟步长计算速度
    # 埃米转玻尔半径, fs 转 a.u., 与output.dat数据一致
    print(
        "差分法计算速度(AU: Ang -> Bohr, fs -> a.u.)\n",
        (geom_last1_xyz - geom_last2_xyz) / 0.5291772083 / (0.5 / 2.418884326e-2),
    )
    print("output.dat中的速度(AU)\n", veloc_last1_dat)
    # 计算结果吻合(计算方法本身存在误差)


# valid_data()


def read_output_dat(
    output_data_file: str = "output.dat",
) -> tuple[dict, dict, dict]:
    """读取 `output.dat` 中的所有信息

    Args:
        output_data_file (str, optional): output_data_file (str, optional): `output.dat` 文件路径. Defaults to 'output.dat'.

    Returns:
        tuple[dict, dict, dict]: (SHARC 设置, 头数据, 轨迹数据)

    Example:

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
    with open(output_data_file, "r", encoding="utf-8") as f:
        line = f.readline()

        settings = {}
        while "End of settings" not in line:
            parts = line.split()
            key = parts[0]
            # value = [eval(n) for n in parts[1:]]
            value = [float(n) if "." in n else int(n) for n in parts[1:]]
            if len(value) == 1:
                settings[key] = value[0]
            else:
                settings[key] = value

            line = f.readline()

        # print(settings)

        line = f.readline()

        header_array_data = {}
        while "End of header array data" not in line:
            if line.startswith("!"):
                key = line[1:].strip()
                value = []
                for _ in range(settings["natom"]):
                    line = f.readline()
                    # print(line)
                    if "." in line:
                        value.append(float(line))
                    else:
                        value.append(line.strip())
                header_array_data[key] = value
                line = f.readline()

        # print(header_array_data)

        line = f.readline()
        # print(line)

        traj_data = {}
        while line:
            if "!" in line:
                # 不是所有条目的数据都会存储，所以要去掉编号
                key = line.split(" ", 2)[-1].strip()
                if traj_data.get(key) is None:
                    traj_data[key] = [[]]
                else:
                    traj_data[key].append([])
            else:
                data = [
                    (
                        float("nan")
                        if v.lower() == "nan"
                        else (float(v) if "." in v else int(v))
                    )
                    for v in line.split()
                ]
                traj_data[key][-1].append(data if len(data) > 1 else data[0])

            line = f.readline()

        # print(traj_data)

        # # 转换数据格式
        for key, value in traj_data.items():
            traj_data[key] = np.array(value)

        # # 打印各个量的维度
        # for key, value in traj_data.items():
        #     print(key, value.shape)

    return settings, header_array_data, traj_data


# def test_read_output_dat():
#     sharc_setting, header_data, traj_data = read_output_dat()

#     natom = sharc_setting["natom"]
#     atoms = header_data["Elements"]
#     geometry = traj_data["Geometry in a.u."]
#     velocity = traj_data["Velocities in a.u."]


class TrajData:
    """
    HDF5 File
        Attributes
            max_bond_length
            xyz_path
        Groups
            Doublet_0
                TRAJ_0001
                    Attributes
                        nstep
                        nfrag
                        ...
                    Dataset
                        velocity
                        ...
                TRAJ_000x

    `output_data/expec.out` 中的数据不会用于 Dalitz  Plot 和 Newton Plot
    """

    # 数据来源 output.dat `output_data/expec.out`
    # 格式 dict[k-v,nd]
    def __init__(self, dat_path: str = "output.dat"):
        # 所有路径都是相对路径, 锚定文件 output.dat
        # 一定存在
        self.dat_path = dat_path
        self.mtime = int(os.path.getmtime(dat_path))  # 修改时间
        # 不一定存在，尤其是 `out_path`
        # 根据相对路径关系生成输出文件的路径
        dirname = os.path.dirname(dat_path)
        self.dir_path = dirname
        self.xyz_path = os.path.join(dirname, "output.xyz")
        self.lis_path = os.path.join(dirname, "output.lis")
        self.out_path = os.path.join(dirname, "output_data/expec.out")

    # 读取 SHARC 的输出文件：`output.xyz`, `output.lis`, `output.dat`, `output_data/expec.out`
    # --------------------------------------------------------------------------------------
    def xyz(self) -> xyztraj.Traj:
        """读取 `output.xyz` 文件 (单位为 `Ang` )

        Returns:
            xyztraj.Traj: 轨迹对象。后续可以直接使用 Mol 计算片段划分与初末动能差值
        """
        assert os.path.exists(self.xyz_path)
        return xyztraj.Traj.from_xyz(self.xyz_path)

    def lis(self) -> pd.DataFrame:
        """读取 `output.lis` 文件

        Returns:
            pd.DataFrame: 列表文件 `output.lis` 的内容
        """
        assert os.path.exists(self.lis_path)
        df = pd.read_csv(
            self.lis_path,
            delimiter=r"\s+",
            skiprows=3,
            # `output.lis` 的格式是固定的
            names=[
                "step",
                "time",
                "state_diag",
                "state_mch",
                "energy_kin",
                "energy_pot",
                "energy_tot",
                "angular_momentum",
                "gradient_rms",
                "density_total",
                "expectation_dm",
                "expectation_s",
                "runtime",
            ],
        )

        return df

    def dat(self) -> tuple:
        """读取 `output.dat` 文件(单位为原子单位 `AU` )

        Returns:
            tuple: 必要的数据
        """
        assert os.path.exists(self.dat_path)
        sharc_setting, header_data, traj_data = read_output_dat(self.dat_path)

        natom = sharc_setting["natom"]
        atoms = header_data["Elements"]
        mass = header_data["Atomic masses"]
        geometry = traj_data["Geometry in a.u."]
        velocity = traj_data["Velocities in a.u."]

        return natom, atoms, np.array(mass), geometry, velocity

        # return dict(natom=natom, atoms=atoms, mass=mass, geometry=geometry, velocity=velocity)

    def expec(self) -> pd.DataFrame:
        """先调用 `SHARC` 生成文件 `output_data/expec.out` , 输出(含错误信息)会被屏蔽, 再读取 `output_data/expec.out` 文件

        Raises:
            RuntimeError: 调用 `SHARC` 生成文件 `output_data/expec.out` 失败

        Returns:
            pd.DataFrame: `output_data/expec.out` 内容
        """
        res = os.system(
            f"cd {self.dir_path}; $SHARC/data_extractor.x output.dat > /dev/null 2>&1 "
        )
        if res != 0:
            logger.warning("Failed to generate expec.out in %s", self.dir_path)
            raise RuntimeError
        else:
            with open(self.out_path, "r", encoding="utf-8") as f:
                f.readline()
                title_row = f.readline()
            f.close()
            names = [
                col_header.strip().replace(" ", "")
                for col_header in title_row[1:].split("|")[:-1]
            ]

            df = pd.read_csv(self.out_path, delimiter=r"\s+", skiprows=3, names=names)

            return df


# 为了控制数据源头的数量，便于后续排查，约定只使用两个文件: `output.dat` , `output_data/expec.out`
def single_trajdata(
    analysis_params: dict = dict(dat_path="output.dat", max_bond_length=2.5)
) -> tuple[dict, dict]:
    """从一条轨迹中提取所有所需数据

    提取出来的数据的具体内容将由此函数决定

    Args:
        analysis_params (dict, optional): 提取数据时的参数，包括 `output.xyz` 文件路径以及片段划分时候的最大键长. Defaults to dict(dat_path='output.xyz',max_bond_length=2.5).

    Returns:
        tuple[dict,dict]: `(轨迹整体属性, 原子速度)` , 约定第一个数据以 `attr` 的形式存储，第二个数据以 `dataset` 形式存储
    """
    attrs, datas = {}, {}

    try:

        trajdata = TrajData(analysis_params["dat_path"])
        mtime = trajdata.mtime
        dir_path = trajdata.dir_path
        state = trajdata.xyz_path.split("/")[-3]
        multiplicity = state.split("_")[0]
        attrs.update(
            dict(mtime=mtime, dir_path=dir_path, state=state, multiplicity=multiplicity)
        )

        # 使用 `output.dat` 中的原子质量以获得和SHARC输出文件中更接近的结果
        natom, atoms, mass, geometry, velocity = trajdata.dat()
        nstep = len(geometry)

        mol_first = xyztraj.Mol(natom, atoms, geometry[0] * xyztraj.BOHR_IN_ANGSTROM)
        frag_info, frag_list = mol_first.frag(index_start=0)
        fragidx_first = [frag["index"] for frag in frag_info.values()]
        ekin_first = mol_first.cal_frag_ekin(
            fragidx_first, mass, velocity[0], factor="au"
        )

        mol_last = xyztraj.Mol(natom, atoms, geometry[-1] * xyztraj.BOHR_IN_ANGSTROM)
        frag_info, frag_list = mol_last.frag(index_start=0)
        fragidx_list = [frag["index"] for frag in frag_info.values()]
        nfrag = len(frag_list)

        ekin_last = mol_last.cal_frag_ekin(
            fragidx_list, mass, velocity[-1], factor="au"
        )

        delta_energy = sum(ekin_last) - sum(ekin_first)

        attrs.update(
            dict(
                natom=natom,
                elements=atoms,
                atomic_masses=mass,
                delta_energy=delta_energy,
                frag=str(frag_info),
                frag_list=frag_list,
            ),
            nstep=nstep,
            nfrag=nfrag,
        )
        datas.update(dict(geometry=geometry, velocity=velocity))

        expec = trajdata.expec()
        attrs.update(dict(expec_columns=list(expec.columns)))
        datas.update(dict(expec_data=expec.to_numpy()))

        # pp(attrs)
        # pp(datas)
    except:
        print(f"Failed {dat_path}")
    return attrs, datas


# single_trajdata(dict(dat_path="Singlet_29/TRAJ_00010/output.dat", max_bond_length=2.5))

# single_trajdata(dict(dat_path="State_1/TRAJ_00001/output.dat", max_bond_length=2.5))


def parallel_trajdata(
    dat_paths: list = sorted(glob("*let_*/TRAJ_*/output.dat"))[:10],
    max_bond_length: float = 2.5,
) -> list[tuple[dict, dict]]:
    """并行提取数据

    Args:
        dat_paths (list, optional): 需要提取数据的 `output.dat` 文件路径列表. Defaults to sorted(glob("*let_*/TRAJ_*/output.dat"))[:10]. 即只是用前十个数据进行测试
        max_bond_length (float, optional): 划分片段时的最大键长. Defaults to 2.5.

    Returns:
        list[tuple[dict,dict]]: 提取到的数据, 每一个元素都是 `single_trajdata` 函数输出的格式
    """
    ncpu = max(psutil.cpu_percent(interval=1, percpu=True).count(0), 1)
    ncpu = min(ncpu, len(dat_paths))

    print(f"{ncpu} cores used for extracting data.")

    with mp.Pool(ncpu) as p:
        trajdatas = p.map(
            single_trajdata,
            (dict(dat_path=dp, max_bond_length=max_bond_length) for dp in dat_paths),
        )

    return trajdatas


# parallel_trajdata()

# # HDF 读取测试
# nstep = []

# f = h5py.File("Data.hdf5", "r")
# # 读取分析参数
# analysis_params = dict(f.attrs.items())
# # 遍历每条轨迹并分别读取attribute和dataset
# for multiplicity in f.keys():
#     for trajectory in f[multiplicity].keys():
#         try:
#             # # f[multiplicity][trajectory] 是一个 group
#             # trajattr = dict(f[multiplicity][trajectory].attrs.items())
#             # # trajdata = dict(f[multiplicity][trajectory].items())
#             # # 读取dataset数据需要用[:]获取数据，否则将返回HDF对象
#             # trajdata = dict(
#             #     (key, value[:])
#             #     for key, value in f[multiplicity][trajectory].items()
#             # )
#             # trajdatas.append((trajattr, trajdata))

#             # nstep.append(f[multiplicity][trajectory].attrs["nstep"])

#             if f[multiplicity][trajectory].attrs["nstep"] > 400 and f[multiplicity][trajectory].attrs["nfrag"] == 3:
#                 nstep.append(f[multiplicity][trajectory].attrs["nstep"])

#         except:
#             warnings.warn(
#                 f"Failed to read data of {multiplicity}/{trajectory} from {self.hdf5_file}.",
#                 UserWarning,
#             )
# f.close()


# # plt.hist(nstep,10)
# # plt.xlabel('nstep')
# # plt.ylabel("Count")

# # plt.show()

# print(len(nstep))


class Dataset:
    """`HDF5` 格式的数据文件

    主要功能

        - 读取全部数据或者读取部分数据

        - 更新数据

        - 创建文件

    Note:
        * 单一进程读写数据文件
        * 路径具有主键的作用
        * 保持一定的拓展性，约定数据的格式，但不限制具体有哪些数据：格式固定(相对路径由多重度+轨迹标号构成)，内容自由(属性和dataset内容可自行增加)
    """

    def __init__(self, hdf5_file: str = "Data.hdf5"):
        self.hdf5_file = hdf5_file

    def create(
        self, trajdatas: list[tuple[dict, dict]], analysis_params: dict = dict()
    ) -> None:
        """创建数据文件

        Args:
            trajdatas (list[tuple[dict,dict]]): `parallel_trajdata` 提取到的轨迹数据
            analysis_params (dict, optional): 轨迹分析参数. Defaults to dict().

        Note:
            需要使用 `dir_path` 属性建立group
        """
        # 如果数据文件存在，则删除
        if os.path.exists(self.hdf5_file):
            warnings.warn(
                f"File {self.hdf5_file} already exists and will be deleted.",
                UserWarning,
            )
            os.remove(self.hdf5_file)
        # 创建数据文件并写入数据
        f = h5py.File(self.hdf5_file, "w")
        # 如果分析参数不为空，保存分析参数
        if analysis_params:
            for key, value in analysis_params.items():
                f.attrs[key] = value
        # 保存每条轨迹的数据，不随时间变化的数据存储为attribute，反之则存储为dataset
        for attrs, datas in trajdatas:
            # 尝试逐条轨迹写入数据
            try:
                # print(attrs, data)
                # 创建group对应一条轨迹
                dir_group = f.create_group(attrs["dir_path"])
                # 存储attribute数据
                for key, value in attrs.items():
                    dir_group.attrs[key] = value
                # 存储dataset数据
                for key, value in datas.items():
                    # 将需要存储的数据转换为numpy数组格式
                    dir_group.create_dataset(key, data=np.array(value))
                    # dir_group.create_dataset(key, data=value)
            except:
                warnings.warn(
                    f"Failed to save data from {attrs['dir_path']} to {self.hdf5_file}.",
                    UserWarning,
                )
        f.close()
        print(f"The file {self.hdf5_file} created.")

    def append(self, trajdatas) -> None:
        """直接向 HDF5 文件中追加数据。**注意: 不检查 `analysis_params` 以及是追加数据是否与已有数据重复**

        Args:
            trajdatas (list[tuple]): 需要追加的数据
        """
        f = h5py.File(self.hdf5_file, "r+")
        for attrs, datas in trajdatas:
            # 尝试逐条轨迹写入数据
            try:
                # print(attrs, data)
                # 创建group对应一条轨迹
                dir_group = f.create_group(attrs["dir_path"])
                # 存储attribute数据
                for key, value in attrs.items():
                    dir_group.attrs[key] = value
                # 存储dataset数据
                for key, value in datas.items():
                    # 将需要存储的数据转换为numpy数组格式
                    dir_group.create_dataset(key, data=np.array(value))
                    # dir_group.create_dataset(key, data=value)
            except Exception as e:
                warnings.warn(
                    f"Failed to save data from {attrs['dir_path']} to {self.hdf5_file}. {e}",
                    UserWarning,
                )
        f.close()

    def read(self) -> tuple:
        """读取数据文件中的所有数据

        Returns:
            tuple: 轨迹数据 `trajdatas` 与 分析参数 `analysis_params`
        """
        trajdatas = []
        analysis_params = {}
        if os.path.exists(self.hdf5_file):
            f = h5py.File(self.hdf5_file, "r")
            # 读取分析参数
            analysis_params = dict(f.attrs.items())
            # 遍历每条轨迹并分别读取attribute和dataset
            for multiplicity in f.keys():
                for trajectory in f[multiplicity].keys():
                    try:
                        # f[multiplicity][trajectory] 是一个 group
                        trajattr = dict(f[multiplicity][trajectory].attrs.items())
                        # trajdata = dict(f[multiplicity][trajectory].items())
                        # 读取dataset数据需要用[:]获取数据，否则将返回HDF对象
                        trajdata = dict(
                            (key, value[:])
                            for key, value in f[multiplicity][trajectory].items()
                        )
                        trajdatas.append((trajattr, trajdata))
                    except:
                        warnings.warn(
                            f"Failed to read data of {multiplicity}/{trajectory} from {self.hdf5_file}.",
                            UserWarning,
                        )
            f.close()
        else:
            warnings.warn(f"File {self.hdf5_file} is not found.", UserWarning)
        return trajdatas, analysis_params

    def mtime(self):
        """读取数据文件中所有轨迹的修改时间

        Returns:
            dict: `{xyz_path: mtime}`
        """
        file_mtime = {}
        if os.path.exists(self.hdf5_file):
            f = h5py.File(self.hdf5_file, "r")
            for multiplicity in f.keys():
                for trajectory in f[multiplicity].keys():
                    try:
                        # 读取文件路径和修改时间
                        file_path = f[multiplicity][trajectory].attrs.get("dat_path")
                        mtime = f[multiplicity][trajectory].attrs.get("mtime")
                        # print(type(mtime))
                        # xyz_mtime[file_path] = int(mtime)
                        file_mtime[file_path] = mtime
                    except:
                        warnings.warn(
                            f"Failed to read `mtime` of {multiplicity}/{trajectory} from {self.hdf5_file}.",
                            UserWarning,
                        )
        else:
            warnings.warn(f"File {self.hdf5_file} is not found.", UserWarning)
        # print(xyz_mtime)
        return file_mtime

    def analysis_params(self) -> dict:
        """读取数据文件中的分析参数

        Returns:
            dict: 分析参数
        """
        analysis_params = dict()
        if os.path.exists(self.hdf5_file):
            try:
                f = h5py.File(self.hdf5_file, "r")
                # 读取分析参数
                analysis_params = dict(f.attrs.items())
                f.close()
            except:
                warnings.warn(
                    f"Failed to extract analysis_params from file {self.hdf5_file}.",
                    UserWarning,
                )
        else:
            warnings.warn(f"File {self.hdf5_file} is not found.", UserWarning)
        return analysis_params

    def update(
        self, trajdatas: list[tuple[dict, dict]], analysis_params: dict = dict()
    ) -> None:
        """更新数据文件

        Args:
            trajdatas (list[tuple[dict,dict]]): `parallel_trajdata` 提取到的轨迹数据
            analysis_params (dict, optional): 轨迹分析参数. Defaults to dict().

        Note:
            需要使用 `mtime` 判断是否更新

            不再使用 `mtime` 判断是否需要进行修改。如果存在则删除再重写，如果不存在则直接写入。
        """
        # 数据文件若不存在则创建文件
        if os.path.exists(self.hdf5_file):
            analysis_params_old = self.analysis_params()
            # 如果分析参数改变则重新创建文件
            if analysis_params == analysis_params_old:
                f = h5py.File(self.hdf5_file, "r+")
                for attrs, datas in trajdatas:
                    # 尝试逐条轨迹写入数据
                    try:
                        dir_group_old = f.get(attrs["dir_path"])
                        # if dir_group_old:
                        # if dir_group_old.attrs.get("mtime") != attrs["mtime"]:
                        # del f[attrs["dir_path"]]
                        # logger.info(
                        # f"Delete data from {attrs['dir_path']} in {self.hdf5_file}."
                        # )
                        # else:
                        # pass
                        # # 不存在此轨迹的数据，需要更新

                        # 存在则删除，不存在则直接写入
                        if dir_group_old:
                            del f[attrs["dir_path"]]
                            print(
                                f"Delete data from {attrs['dir_path']} in {self.hdf5_file}."
                            )
                        else:
                            pass
                            # 不存在此轨迹的数据，需要更新

                        # print(attrs, data)
                        # 创建group对应一条轨迹
                        dir_group = f.create_group(attrs["dir_path"])
                        # 存储attribute数据
                        for key, value in attrs.items():
                            dir_group.attrs[key] = value
                        # 存储dataset数据
                        for key, value in datas.items():
                            dir_group.create_dataset(key, data=np.array(value))
                    except:
                        warnings.warn(
                            f"Failed to save data from {attrs['dir_path']} to {self.hdf5_file}.",
                            UserWarning,
                        )
                f.close()
            else:
                warnings.warn(
                    f"Parameters of file {self.hdf5_file} has changed so this file will be deleted.",
                    UserWarning,
                )
                self.create(trajdatas, analysis_params)
        else:
            warnings.warn(
                f"File {self.hdf5_file} not found and will be create.", UserWarning
            )
            self.create(trajdatas, analysis_params)

    def fragdata(self) -> list[dict]:
        """读取数据文件中的所有轨迹的 **属性数据** , 用于 MolFragApp 的表格展示

        Returns:
            list[dict]: 轨迹的数据
        """
        fragdata = []
        try:
            f = h5py.File(self.hdf5_file, "r")
            for multiplicity in f.keys():
                for trajectory in f[multiplicity].keys():
                    try:
                        fragdata_traj = dict(f[multiplicity][trajectory].attrs.items())
                        # 增加一项新的数据。
                        # 后续调用数据时候需要更换 file -> xyz_path
                        fragdata_traj["frag_string"] = " + ".join(
                            fragdata_traj["frag_list"]
                        )
                        #
                        fragdata.append(fragdata_traj)
                    except Exception as e:
                        warnings.warn(
                            f"Failed to extract data from group {multiplicity + '/' + trajectory} in file {self.hdf5_file}. {e}",
                            UserWarning,
                        )
            f.close()
        except FileNotFoundError:
            warnings.warn(
                f"File {self.hdf5_file} is not found. Please create this file before using it.",
                UserWarning,
            )
        return fragdata

    def expec_out(self, path: str) -> pd.DataFrame:
        """根据路径获取 `expec.out` 文件中包含的数据

        Args:
            path (str): `expec.out` 文件对应的路径，例如 `State_1/TRAJ_001`

        Returns:
            pd.DataFrame: `expec.out` 文件中包含的数据
        """
        expec_out = None
        try:
            f = h5py.File(self.hdf5_file, "r")
            expec_out = pd.DataFrame(
                data=f[path]["expec_data"][:], columns=f[path].attrs["expec_columns"]
            )
            f.close()
            return expec_out
        except Exception as e:
            print(e)
            return pd.DataFrame(expec_out)

    # def coord(self, path:str):
    # coord = None
    # try:
    # f = h5py.File(self.hdf5_file, "r")
    # coord = f[path]["geometry"][:]*xyztraj.BOHR_IN_ANGSTROM
    # f.close()
    # except Exception as e:
    # print(e)
    # return coord

    def traj4show(self, path: str) -> tuple[xyztraj.Traj, pd.DataFrame, Any] | None:
        """获取 `MolFragApp` 中展示轨迹动画与能量曲线所需的数据

        Args:
            path (str): 对应的路径，例如 `State_1/TRAJ_001`

        Returns:
            tuple[xyztraj.Traj, pd.DataFrame, Any] | None: 示轨迹动画与能量曲线所需的数据
        """
        traj = None
        try:
            f = h5py.File(self.hdf5_file, "r")
            traj = (
                xyztraj.Traj(
                    natom=f[path].attrs["natom"],
                    atoms=f[path].attrs["elements"],
                    coord=f[path]["geometry"][:] * xyztraj.BOHR_IN_ANGSTROM,
                    # 原子单位的坐标转换为埃米
                ),
                pd.DataFrame(
                    data=f[path]["expec_data"][:],
                    columns=f[path].attrs["expec_columns"],
                ),
                eval(f[path].attrs["frag"]),
            )
            f.close()
        except Exception as e:
            print(e)
        return traj

    # def get(self, path:str):
    # data = None
    # try:
    # f = h5py.File(self.hdf5_file, "r")
    # data = dict(f[path].attrs.items())
    # data.update(dict((key, value[:]) for key, value in f[path].items()))
    # except Exception as e:
    # print(e)
    # return data

    # def frag_ek_p(self, paths:list):
    # # 使用最后一帧的片段划分
    # ekin_frag, momentum_frag = [], []

    # for path in paths:
    # try:
    # f = h5py.File(self.hdf5_file, "r")
    # frag_info, mass, veloc = eval(f[path].attrs["frag"]), f[path].attrs["atomic_masses"], f[path]["velocity"][-1]
    # momentum_atom = np.einsum("ij,i->ij", veloc, mass)
    # fragidx_list = [fv["index"] for fv in frag_info.values()]
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
    # print(e)
    # return np.array(ekin_frag), np.array(momentum_frag)

    def frag_ek_p(
        self, paths: list[str], idx: int, frag_orders: list[str]
    ) -> tuple[np.ndarray, np.ndarray]:
        """按照最后一帧的片段划分计算指定路径的轨迹的在特定帧的片段动能与动量

        Args:
            paths (list[str]): 指定需要计算动量和动能的轨迹的路径
            frag_orders (list[str]): 片段的顺序，每个元素由化学式和原子索引构成，例如 `H2O | [2, 1, 0]`
            idx (int, optional): 指定计算哪一帧的动量和动能. Defaults to -1.

        Returns:
            tuple[np.ndarray, np.ndarray]: 片段动能 `(num_selected, 3)` 与动量 `(num_selected, 3, 3)`
        """
        # 使用最后一帧的片段划分
        ekin_frag, momentum_frag = [], []

        for path in paths:
            try:
                f = h5py.File(self.hdf5_file, "r")
                frag_info, mass, veloc = (
                    eval(f[path].attrs["frag"]),
                    f[path].attrs["atomic_masses"],
                    f[path]["velocity"][idx],
                )
                frag_info = dict(
                    (f"{fv['formula']} | {sorted(fv['index'])}", fv["index"])
                    for fv in frag_info.values()
                )
                fragidx_list = [frag_info.get(formula) for formula in frag_orders]
                # fragidx_list = [[0],[1],[2]]
                momentum_atom = np.einsum("ij,i->ij", veloc, mass)
                momentum_frag.append(
                    [momentum_atom[fragidx].sum(axis=0) for fragidx in fragidx_list]
                )
                ekin_frag.append(
                    [
                        np.sum((momentum_atom[fragidx].sum(axis=0)) ** 2)
                        / (2 * mass[fragidx].sum())
                        for fragidx in fragidx_list
                    ]
                )
                f.close()
            except Exception as e:
                print(f"{path} {e}")
        return np.array(ekin_frag), np.array(momentum_frag)

    def union(self, other) -> None:
        """对两个HDF5文件中数据取并集, 相同的新数据不会覆盖(现实警告信息)

        Args:
            other (Dateset): 另一个 `Dataset` 对象
        """
        assert (
            self.analysis_params()["max_bond_length"]
            == other.analysis_params()["max_bond_length"]
        ), "Different `max_bond_length`"
        union_file = "%s_%s_union.hdf5" % (
            self.hdf5_file.rsplit(".", 1)[0],
            other.hdf5_file.rsplit(".", 1)[0],
        )
        Dataset(union_file).create(
            self.read()[0] + other.read()[0], self.analysis_params()
        )

    def add(self, other) -> None:
        """对两个HDF5文件中数据进行合并, 重复 `key` 添加随机序列后缀

        Args:
            other (Dateset): 另一个 `Dataset` 对象
        """
        # 内容合并，
        assert (
            self.analysis_params()["max_bond_length"]
            == other.analysis_params()["max_bond_length"]
        ), "Different `max_bond_length`"
        add_file = "%s_%s_add.hdf5" % (
            self.hdf5_file.rsplit(".", 1)[0],
            other.hdf5_file.rsplit(".", 1)[0],
        )
        trajdatas, analysis_params = self.read()
        trajdatas_add = other.read()[0]

        for attrs_add, datas_add in trajdatas_add:
            # print(attrs_add["dir_path"])
            if attrs_add["dir_path"] in (attrs["dir_path"] for attrs, _ in trajdatas):
                while True:
                    key_add_new = attrs_add["dir_path"] + f"_{randseq()}"
                    if key_add_new not in (attrs["dir_path"] for attrs, _ in trajdatas):
                        attrs_add["dir_path"] = key_add_new
                        break
            trajdatas.append((attrs_add, datas_add))

        Dataset(add_file).create(trajdatas, analysis_params)


class TestDataset:

    # def test_read_create(self):
    # # 读取原有数据并创建用于测试的文件
    # ds = Dataset('Data.hdf5')
    # pp(ds.read()[0][-1])
    # trajdatas, analysis_params = ds.read()
    # d2 = Dataset('Data2.hdf5')
    # d2.create(trajdatas[:3]+trajdatas[-3:], analysis_params)

    def test_read(self):
        ds = Dataset("Data.hdf5")
        trajdatas, analysis_params = ds.read()
        pp(trajdatas)
        pp(analysis_params)

    def test_analysis_params(self):
        ds = Dataset("Data.hdf5")
        analysis_params = ds.analysis_params()
        pp(analysis_params)

    def test_fragdata(self):
        ds = Dataset("Data.hdf5")
        fragdata = ds.fragdata()
        pp(fragdata)

    def test_create(self):
        ds = Dataset("Data.hdf5")
        trajdatas, analysis_params = ds.read()

        dsc = Dataset("Data_created.hdf5")
        dsc.create(trajdatas, analysis_params)

    def test_update(self):
        ds = Dataset("Data.hdf5")
        trajdatas, analysis_params = ds.read()

        dsc = Dataset("Data_update.hdf5")
        dsc.create(trajdatas[:2] + trajdatas[-2:], analysis_params)

        trajdatas[2][0]["mtime"] = 1024
        dsc.update(trajdatas, analysis_params)

    # def test_read_output_dat(self):
    #     sc, hd, td = TrajData.read_output_dat()
    #     pp(sc)
    #     pp(hd)
    #     pp(td)

    def test_full(self):
        analysis_params = dict(
            max_bond_length=2.5, max_nstep=2000, dat_pattern="*let_*/TRAJ_*/output.dat"
        )
        Dataset().create(
            parallel_trajdata(
                sorted(glob(analysis_params["dat_pattern"]))[:10],
                analysis_params["max_bond_length"],
            ),
            analysis_params,
        )

    def test_batch_extract(self):
        analysis_params = dict(
            max_bond_length=2.5, max_nstep=2000, dat_pattern="State_*/TRAJ_*/output.dat"
        )
        # 批量提取数据并写入文件
        HDF_FILE = "Data.hdf5"
        ds = Dataset(HDF_FILE)
        BATCH_SIZE = 40
        dat_files = sorted(glob(analysis_params["dat_pattern"]))
        for i in range(0, len(dat_files), BATCH_SIZE):
            print(f"PROCESS {i*BATCH_SIZE} - {(i+1)*BATCH_SIZE}")
            if not os.path.exists(HDF_FILE):
                ds.create(
                    parallel_trajdata(
                        dat_files[i : i + BATCH_SIZE],
                        analysis_params["max_bond_length"],
                    ),
                    analysis_params,
                )
            else:
                ds.append(
                    parallel_trajdata(
                        dat_files[i : i + BATCH_SIZE],
                        analysis_params["max_bond_length"],
                    )
                )


if __name__ == "__main__":
    import time

    time_start = time.time()

    analysis_params = dict(
        max_bond_length=2.5, max_nstep=2000, dat_pattern="State_*/TRAJ_*/output.dat"
    )
    # 批量提取数据并写入文件
    HDF_FILE = "Data.hdf5"
    ds = Dataset(HDF_FILE)
    BATCH_SIZE = 100
    dat_files = sorted(glob(analysis_params["dat_pattern"]))
    for i in range(0, len(dat_files), BATCH_SIZE):
        print(f"PROCESS {i} - {i+BATCH_SIZE}")
        if not os.path.exists(HDF_FILE):
            ds.create(
                parallel_trajdata(
                    dat_files[i : i + BATCH_SIZE],
                    analysis_params["max_bond_length"],
                ),
                analysis_params,
            )
        else:
            ds.append(
                parallel_trajdata(
                    dat_files[i : i + BATCH_SIZE],
                    analysis_params["max_bond_length"],
                )
            )

    time_stop = time.time()
    print(f"{len(dat_files)} files / {time_stop - time_start:.2f} sec")
