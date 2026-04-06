import multiprocessing as mp
import time
import re
import os

# def sleep_fun(n):
    # print(f'{n} start')
    # time.sleep(n)
    # print(f'{n} stop')

    
    
def parse_all_run_traj(filename:str='all_run_traj.sh'):
    with open(filename,'r') as f:
        content = f.read()
    f.close()
    
    cwd = re.search("CWD=(.*?)\n",content).group(1)
    
    # print(cwd)
    
    # print(re.findall("cd \$CWD.*?\nbash run.sh",content))
    run_script = re.findall("cd \$CWD.*?DONE",content,re.S)
    run_rows = []
    
    for row in run_script:
        run_rows.append(row.replace('\n',';').replace('$CWD',cwd))
        
    return run_rows
    
def run_one_row(row):
    res = os.system(row)
    
    # print(f"START{row}")
    # time.sleep(3)
    # print(f"STOP {row}")
    # res = 1
    return res 
        
def run_mp(filename:str='all_run_traj.sh',num_process:int=0):

        
    inp = parse_all_run_traj(filename)
    # print(inp)

    
    if num_process == 0:
        num_process = mp.cpu_count()
        
    pool = mp.Pool(num_process)
    
    out = pool.map(run_one_row, inp)
    
    return out
    
# def test_run_mp():
    # run_mp(filename='all_run_traj_test.sh',num_process=3)
    
if __name__ == '__main__':
    
    run_mp(num_process=32)

