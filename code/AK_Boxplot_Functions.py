#%%
import numpy as np
from scipy import stats
from Models import *
from Epistasis_calc_functions import *
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
from statistics import mean

data_out = pd.read_excel('../data/OutputDF.xlsx')
data_sen = pd.read_excel('../data/SensorDF.xlsx')
data_reg = pd.read_excel('../data/RegulatorDF.xlsx')
#%%
def get_boxplot(data): 
    df = data
    eps = list(df['Epistasis'])
    gen = list(df['Genotype'])
    if gen[0].find('O') != -1:
        x = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'O10']
    elif gen[0].find('S') != -1:
        x = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10']
    elif gen[0].find('R') != -1:
        x = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10']
    x_ticks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    # list instead of df
    gen1 = eps[0:360]
    gen2 = eps[360:720]
    gen3 = eps[720:1080]
    gen4 = eps[1080:1440]
    gen5 = eps[1440:1800]
    gen6 = eps[1800:2160]
    gen7 = eps[2160:2520]
    gen8 = eps[2520:2880]
    gen9 = eps[2880:3240]
    gen10 = eps[3240:3600]
    new_df = [gen1,gen2,gen3,gen4,gen5,gen6,gen7,gen8,gen9,gen10]
    a = plt.boxplot(new_df, showfliers=False, patch_artist=False)
    l, r = plt.xlim()
    plt.hlines(y= 0, xmin= l, xmax= r, linestyles='dashed', linewidth = 1, colors='k')
    for i in range(10):
        y1 = np.array(new_df[i][0:60])
        y2 = np.array(new_df[i][60:360])
        x1 = np.array(np.random.normal(1+i, 0.04, size=len(y1)))
        x2 = np.array(np.random.normal(1+i, 0.04, size=len(y2)))
        plt.scatter(x1, y1, color = 'red', marker='o', alpha=0.2, linewidth = 1.5, edgecolors='black')
        plt.scatter(x2, y2, color = 'blue', marker='o', alpha=0.2, linewidth = 1.5,edgecolors='black')
    plt.xticks(x_ticks, x)

#%%        
def get_allBoxplots(data1, data2, data3):
   fig = plt.figure()
   plt.subplot(3,1,1)
   get_boxplot(data_out)
   plt.subplot(3,1,2)
   get_boxplot(data_sen)
   plt.ylabel('Epistasis', fontsize=10)
   plt.subplot(3,1,3)
   get_boxplot(data_reg)
   plt.xlabel('Genotypes', fontsize=10)


#%%
get_allBoxplots(data_out, data_sen, data_reg)