import os
from glob import glob
import numpy as np
import itertools
import pandas as pd 
import streamlit as st
import matplotlib.pyplot as plt 



def reorder(arr,max_delta=0.5):
    for i in range(1,len(arr)):
        # 任何一个差值大于阈值就进行排序, 为了处理连续异常，mask必需在每次循环前重算
        mask = np.abs(arr[i] - arr[i-1]) > max_delta
        if mask.sum() >= 1:
            # 计算当前步所有排列与上一步的差值，选择差值最小的排列作为排序后的结果
            all_comb = itertools.permutations(arr[i],len(arr[i]))
            delta_comb = [[np.abs(np.array(comb)-arr[i-1]).sum(),np.array(comb)] for comb in all_comb]
            delta_comb.sort(key=lambda item:item[0])
            arr[i] = delta_comb[0][1]
    
    return arr

def get_all_traj_state(lis_pattern="*let_*/TRAJ_*/output.lis", save_data = True):
    traj_state = {}
    
    for lis_file in sorted(glob(lis_pattern)):
        #print(lis_file,os.path.dirname(lis_file))
        data = np.loadtxt(lis_file)
        #print(data[:,2]) # State diag
        traj_state[os.path.dirname(lis_file)] = dict(time=data[:,1],state=data[:,2].astype(int))
    
    if save_data:
        np.save("traj_state.npy", traj_state)
        
    return traj_state


def get_traj_state(lis_file, save_data = True):
    # >>> get_traj_state("Singlet_0/TRAJ_00001/output.lis")
    data = np.loadtxt(lis_file)
    traj_state = dict(time=data[:,1],state=data[:,2].astype(int))
    print(traj_state)
    if save_data:
        np.save("traj_state.npy", traj_state)
    return traj_state

def load_traj_state(npy_file='traj_state.npy'):
    traj_state = np.load(npy_file, allow_pickle=True)
    return traj_state.item()

# ===============================================

def get_all_state_chg(csv_pattern="*let_*/TRAJ_*/*.csv",save_data=True):
    state_chg = {}
    for csv_file in sorted(glob(csv_pattern)):
        # print(os.path.dirname(csv_file),pd.read_csv(csv_file).columns)
        data = pd.read_csv(csv_file)
        state_chg[os.path.dirname(csv_file)] = data.to_dict('list')
        
    if save_data:
        np.save("state_chg.npy", state_chg)
        
    return state_chg
        
def load_state_chg(npy_file="state_chg.npy"):
    state_chg = np.load(npy_file, allow_pickle=True)
    return state_chg.item()

# ===============================================

def get_traj_chg(traj_state=None,state_chg=None,save_data=True):
    if traj_state == None:
        traj_state = load_traj_state()
    if state_chg == None:
        state_chg = load_state_chg()

    def row_op(row):
        # return row['time'], row[f's{int(row.state)}-chg1'], row[f's{int(row.state)}-chg2'], row[f's{int(row.state)}-chg3']
        return dict(time=row['time'], chg1=row[f's{int(row.state)}-chg1'], chg2=row[f's{int(row.state)}-chg2'], chg3=row[f's{int(row.state)}-chg3'], Epot=row[f'state{int(row.state)}'])
    
    traj_chg = {}
    for traj, time_chg in state_chg.items():
        time_chg = pd.DataFrame(time_chg)
        time_state = pd.DataFrame(traj_state[traj])
        time_all = pd.merge(time_chg,time_state,how='left',on='time')
        # print(time_all)
        # traj_chg = time_all.apply(row_op,axis=1,result_type='expand')
        # print(traj_chg)
        traj_chg[traj] = time_all.apply(row_op,axis=1,result_type='expand').to_dict('list')
        
    if save_data:
        np.save("traj_chg.npy", traj_chg)
        
    return traj_chg
    
def load_traj_chg(npy_file="traj_chg.npy"):
    traj_chg = np.load(npy_file, allow_pickle=True)
    return traj_chg.item()
    
    
# get_traj_chg()
# print(load_traj_chg())

# ===============================================


def display_traj_chg(plot_all_ene=False, plot_all_chg=False, frag_name={'1':'#1','2':'#2','3':'#3'}):
    traj_chg = load_traj_chg()
    
    
    
    traj = st.radio(
        "选择轨迹",
        traj_chg.keys(),
        horizontal = True   
    )
    
    
    fig, axs = plt.subplots(2,1,sharex=True)
    
    traj_plot = pd.DataFrame(traj_chg[traj])
    
    if plot_all_ene:
        state_chg = pd.DataFrame(load_state_chg()[traj])
        
        for col in state_chg.columns:
            if 'state' in col:
                axs[0].plot(state_chg.time, state_chg[col],'k-',color='gray',alpha=0.4)
                
    
    axs[0].plot(traj_plot['time'], traj_plot['Epot'],'.')
    axs[0].set_ylabel('Epot')
    
    if plot_all_chg:
        state_chg = pd.DataFrame(load_state_chg()[traj])
        
        for col in state_chg.columns:
            if 'chg' in col:
                axs[1].plot(state_chg.time, state_chg[col],'k-',color='gray',alpha=0.4)
        
    chg_plot = ['chg1','chg2','chg3']
    # st.write(traj_plot[chg_plot])
    # lines = axs[1].plot(traj_plot['time'],traj_plot[chg_plot],'.')
    lines = axs[1].plot(traj_plot['time'],reorder(traj_plot[chg_plot].to_numpy()),'.')
    # traj_plot[chg_plot].to_csv('special_data.csv')
    # axs[1].legend(lines,chg_plot)
    axs[1].legend(lines,frag_name.values())
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Charge')
    
    st.pyplot(fig)
    
    # 一条轨迹的电荷dalitz 图
    fig, ax = plt.subplots()
    charge_total = 2
    # d1, d2, d3 = traj_plot['chg1']/charge_total, traj_plot['chg2']/charge_total, traj_plot['chg3']/charge_total
    d1_d2_d3 = reorder(traj_plot[chg_plot].to_numpy())/charge_total
    d1, d2, d3 = d1_d2_d3[:,0],d1_d2_d3[:,1],d1_d2_d3[:,2]
    x = (d3 - d2) / 3**0.5
    y = d1 
    
    ax.plot(x,y,color='gray',alpha=0.4)
    # ax.plot(x.iloc[0],y.iloc[0],'r*')
    # ax.plot(x.iloc[-1],y.iloc[-1],'bx')
    ax.plot(x[0],y[0],'r*')
    ax.plot(x[-1],y[-1],'bx')
    
    
    ax.plot([1/3**0.5,0,-1/3**0.5,1/3**0.5],[0,1,0,0],'k-')
    # 三边中点三角形
    ax.plot([0,1/3**0.5/2,-1/3**0.5/2,0],[0,1/2,1/2,0],'k-')
    
    ax.text(0,0,'Charge %s'%frag_name['1'], ha='center',va='top')
    ax.text(1/3**0.5/2,1/2,'Charge %s'%frag_name['2'])
    ax.text(-1/3**0.5/2,1/2,'Charge %s'%frag_name['3'],ha='right')

    ax.axis('equal')
    # ax.legend()
    ax.set_xticks([])
    ax.set_yticks([])
    
    st.pyplot(fig)

display_traj_chg(plot_all_ene=True,plot_all_chg=True,frag_name={'1':'$H^1$','2':'$H^2$','3':'$O^3$'})

# ===============================================


def dalitz_chg(frag_name={'1':'#1','2':'#2','3':'#3'}):
    traj_chg = load_traj_chg()
    charge_total = 2
    
    fig, ax = plt.subplots()
    
    for traj, traj_plot in traj_chg.items():
        traj_plot = pd.DataFrame(traj_plot)
        # d1, d2, d3 = traj_plot['chg1']/charge_total, traj_plot['chg2']/charge_total, traj_plot['chg3']/charge_total
        d1_d2_d3 = reorder(traj_plot[['chg1','chg2','chg3']].to_numpy())/charge_total
        d1, d2, d3 = d1_d2_d3[:,0],d1_d2_d3[:,1],d1_d2_d3[:,2]
        x = (d3 - d2) / 3**0.5
        y = d1 
        
        # st.write(traj.split('/')[0].split('_')[-1])
        traj_color = traj.split('/')[0].split('_')[-1]
        
        # ax.plot(x,y,color='gray',alpha=0.4)
        ax.plot(x,y,color=f'C{traj_color}',alpha=0.8)
        
        # ax.plot(x.iloc[0],y.iloc[0],'r*')
        # ax.plot(x.iloc[-1],y.iloc[-1],'bx')
        
        
        ax.plot(x[0],y[0],'r*',color=f'C{traj_color}')
        ax.plot(x[-1],y[-1],'bx',color=f'C{traj_color}')
        
        # break
        
    # ax.plot([1/3**0.5,0,-1/3**0.5,1/3**0.5],[0,1,0,0],'k-')
    ax.plot([0,1/3**0.5/2,-1/3**0.5/2,0],[0,1/2,1/2,0],'k-')
    
    ax.text(0,0,'Charge %s'%frag_name['1'], ha='center',va='top')
    ax.text(1/3**0.5/2,1/2,'Charge %s'%frag_name['2'])
    ax.text(-1/3**0.5/2,1/2,'Charge %s'%frag_name['3'],ha='right')
    
    ax.axis('equal')
    # ax.legend()
    ax.set_xticks([])
    ax.set_yticks([])
    
    plt.savefig('charge_dalitz.png')
    
    st.pyplot(fig)

dalitz_chg(frag_name={'1':'$H^1$','2':'$H^2$','3':'$O^3$'})


