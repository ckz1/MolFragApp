# MolFragApp

## Introduction

MolFragApp is a Python-based project that uses molecular structure files in 'xyz' format for fragmentation analysis and trajectory visualization.

## Demonstration

### Interface

![](images/overview.png)

### Setting Parameters

![](images/parameter.png)

> - The configuration of [bond length](https://baike.baidu.hk/item/%E9%8D%B5%E9%95%B7/2442392) and ~~number of atoms~~ should be considered carefully based on the specific system and problem.
> - The configuration for the number of atoms was not suitable and has been removed.

### Analysis Results

![](images/result.png)

#### Fragmentation

![](images/frag_split.gif)

#### Trajectory Statistics

![](images/multitraj_stat.gif)

#### Viewing the Animation

![](images/singletraj_geom_ene.gif)

## Usage

### Install

```shell
# Create environment
conda create --name molfrag python=3.11
conda activate molfrag

git clone https://github.com/ckz1/MolFragApp.git

# pip freeze > requirements.txt
# Install dependencies
pip install -r requirements.txt
```

### Run

```shell
# Stop firewall (Linux)
systemctl stop firewalld

# Run directly
streamlit run MolFragApp.py
# Run in background
nohup streamlit run MolFragApp.py > MolFragApp.log 2>&1 &
```

### Example

1. [Install](#install)
2. Unzip the demo.zip file, which contains several trajectory files: [Doublet_0/TRAJ_00004](Doublet_0/TRAJ_00004) and [Doublet_2/TRAJ_00012](Doublet_2/TRAJ_00012)
3. [Run](#run)
4. Modify parameters:
   - Change the trajectory file (`xyz` format) path to: `Singlet_*/TRAJ_*/output.xyz`
   - ~~Change the `template` file path to: `MOLCAS.template`~~
