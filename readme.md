# MolFragApp

## 介绍

MolFragApp 是一个基于Python语言的项目，使用 `xyz`格式的分子结构文件进行碎片分析和轨迹可视化。

MolFragApp is a project based on Python language that uses molecular structure files in 'xyz' format for fragmentation analysis and visualization of trajectories

## 使用效果

### 界面

![](images/overview.png)

### 设置参数

![](images/parameter.png)

> - [键长](https://baike.baidu.hk/item/%E9%8D%B5%E9%95%B7/2442392)和~~原子数~~的设置需要考虑具体的体系和问题。
> - 原子数的设置不合适，已删除。

### 分析结果

![](images/result.png)

#### 片段划分

![](images/frag_split.gif)

#### 轨迹统计

![](images/multitraj_stat.gif)

#### 查看动画

![](images/singletraj_geom_ene.gif)

## 使用方法

### Install

```shell
# 创建环境
conda create --name molfrag python=3.11
conda activate molfrag

git clone 

# pip freeze > requirements.txt
# 安装依赖
pip install -r requirements.txt
```

### Run

```shell
# 关闭防火墙(Linux)
systemctl stop firewalld

# 直接运行
streamlit run MolFragApp.py
# 后台运行
nohup streamlit run MolFragApp.py > MolFragApp.log 2>&1 &
```

### demo

1. [安装](#install)
2. 解压 `demo.zip` 文件，其中包含若干轨迹文件
3. [运行](#run)
4. 修改参数
   - 修改轨迹文件(`xyz`格式)路径为: `Singlet_*/TRAJ_*/output.xyz`
   - 修改 `template` 文件路径为: `MOLCAS.template`
