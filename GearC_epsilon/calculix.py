import os
import multiprocessing
import pandas as pd
def run(directory, rA1, rB1, rl, rD1, eqP, eqW, eqPI, eqAP, eqAW):
    os.chdir(directory)
    
    dfs = pd.read_csv('mm.dfl', quotechar=',')
    dfm = pd.read_csv('allinone.inp', quotechar=',', error_bad_lines=False)
    
    print(dfs)
    print(dfm)
    # with open('mm.dfl', 'r+') as surf:
    #     lines = zip(*[line.split() for line in surf])[0]
    # with open('allinone.inp', 'r') as mesh:
    #     linem = zip(*[line.split() for line in mesh])[0]
    # print(lines)
    # if first_col in first_colm:
        
        # if rcon < rA1:
        #     heat = 0
        # elif rA1 <= rcon < rB1:
        #     heat = np.polyval(eqP[0], rcon)
        # elif rB1 <= rcon < rl[0,0]:
        #     heat = np.polyval(eqP[1], rcon)
        # elif rl[0,0] <= rcon < rD1:
        #     heat = np.polyval(eqP[2], rcon)
        # else:
        #     heat = np.polyval(eqP[3], rcon) 
        # currentline[2] = heat
    # os.environ['OMP_NUM_THREADS'] = str(multiprocessing.cpu_count())
    # os.system("cgx -b pre.fbd")
    # os.system("ccx _main_temp")
    # os.system("monitor.py Hertz")
    # os.system("cgx -b post.fbd")
    # os.system("cgx -v _main_temp.frd")
