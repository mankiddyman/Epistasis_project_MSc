#%%
import numpy as np
from scipy import stats
from Models import *
from Epistasis_calc_functions import *
import matplotlib.pyplot as plt
import statistics
from statistics import mean
# COMMENT
#%%
def get_eData(model):
    # observed epistasis
    if model == 'observed':
        eps_exp_df = get_Eps('observed')
        eps = eps_exp_df['Ep']
    # model hill epistasis    
    elif model == 'model_hill':
        eps_hill_df = get_Eps(model_hill.model)
        eps = eps_hill_df['Ep']
    # elif model == 'other model'
    elif model == 'model_thermodynamic':
        eps_td_df = get_Eps(model_thermodynamic.model)
        eps = eps_td_df['Ep']
    return eps 
#%%
eps_lad = pd.read_excel('../results/Eps_observed.xlsx')['Ep']
eps_hm = pd.read_excel('../results/Eps_model_hill_all.xlsx')['Ep']
eps_td = pd.read_excel('../results/Eps_model_thermodynamic_all.xlsx')['Ep']
eps_hm_sens = pd.read_excel('../results/Eps_model_hill_sensor.xlsx')['Ep']
eps_td_sens = pd.read_excel('../results/Eps_model_thermodynamic_sensor.xlsx')['Ep']
eps_hm_stripe = pd.read_excel('../results/Eps_model_hill_stripe.xlsx')['Ep']
eps_td_stripe = pd.read_excel('../results/Eps_model_thermodynamic_stripe.xlsx')['Ep']
#%%
def get_histograms(eps1, eps2, eps3, eps4, eps5, n, b, strat): # x, y
    option = n
    if option == 1:
        label_name = 'Hill Model'
    elif option == 2:
        label_name = 'Hill Shaky Model'
        # pvalue = float(pvalue2)
    elif option == 3:
        label_name = 'Thermodynamic Model'
        # pvalue = float(pvalue3)
    elif option == 4:
        label_name = 'Competition Model'
    option2 = b
    if option2 == 1:
        label__name = 'Hill Model'
    elif option2 == 2:
        label__name = 'Hill Shaky Model'
        # pvalue = float(pvalue2)
    elif option2 == 3:
        label__name = 'Thermodynamic Model'
        # pvalue = float(pvalue3)
    elif option2 == 4:
        label__name = 'Competition Model'
        # pvalue = float(pvalue4)
    binwidth = 0.1

    d, u = plt.ylim()
    mean_lad = statistics.mean(eps1)
    mean_hm = statistics.mean(eps2)
    mean_td = statistics.mean(eps3)
    #
    mean_td_sens = statistics.mean(eps4)
    mean_hm_sens = statistics.mean(eps5)
    mean_td_sens1 = round(mean_td_sens,3)
    mean_hm_sens1 = round(mean_hm_sens,3)
    mean_lad1 = round(mean_lad,3)
    mean_hm1 = round(mean_hm,3)
    mean_td1 = round(mean_td,3)
    #
    new_mean_lad = 'Mean: ' + str(mean_lad1)
    new_mean_hm = 'Mean: ' + str(mean_hm1)
    new_mean_td = 'Mean: ' + str(mean_td1)
    new_mean_td_sens = 'Mean: ' + str(mean_td_sens1)
    new_mean_hm_sens = 'Mean: ' + str(mean_hm_sens1)

    fig = plt.figure()
    ax = fig.subplots(1,3)
    #plt.show()
    #plt.subplot(1,3,1) 
    ax[0].hist(eps1, bins = np.arange(min(eps1), max(eps1) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5,  label=f'Experimental ,{new_mean_lad}')
    #plt.legend('', frameon=False)
    ax[0].hist(eps2, bins = np.arange(min(eps2), max(eps2) + binwidth, binwidth), alpha=0.5,fill = False, histtype = 'step', edgecolor = 'g', linewidth = 1.5,  label=label_name + f", {new_mean_hm}")
    ax[0].hist(eps5, bins = np.arange(min(eps5), max(eps5) + binwidth, binwidth), alpha=0.5,fill = False, histtype = 'step', edgecolor = 'mediumaquamarine', linewidth = 1.5,  label='Hill Sensor'+ f", {new_mean_hm_sens}")
    ax[0].vlines(x= mean_lad, ymin= d, ymax= 1000, linestyles='dashed', linewidth = 0.6, colors='r')
    ax[0].vlines(x= mean_hm, ymin= d, ymax= 1000, linestyles='dashed', linewidth = 0.6, colors='g')
    ax[0].vlines(x= mean_hm_sens, ymin= d, ymax= 1000, linestyles='dashed', linewidth = 0.6, colors='lime')
    #plt.legend(loc='upper left', bbox_to_anchor=(3.55, 1), frameon=False)
    #plt.ylabel('Mutants', fontsize=10)
    #plt.xlim(-1,1)
    
    #
    #plt.subplot(1,3,2)
    ax[1].hist(eps1, bins = np.arange(min(eps1), max(eps1) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5)
    ax[1].hist(eps3, bins = np.arange(min(eps3), max(eps3) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'b', linewidth = 1.5, label=label__name+ f", {new_mean_td}")
    ax[1].hist(eps4, bins = np.arange(min(eps4), max(eps4) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'navy', linewidth = 1.5,  label='Thermodynamic Sensor'+ f", {new_mean_td_sens}")
    ax[1].vlines(x= mean_lad, ymin= d, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='r')
    ax[1].vlines(x= mean_td, ymin= d, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='b')
    ax[1].vlines(x= mean_td_sens, ymin= d, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='cornflowerblue')
    ax[1].legend(loc='upper left', bbox_to_anchor=(2.29, 0.83), frameon=False)
    ax[1].xlabel('Epistasis', fontsize=10)
    plt.xlim(-1,1)
    plt.title('Distribution of Epistasis of \n'+strat + ' data', fontweight='bold')
    #
    plt.subplot(1,3,3)
    plt.hist(eps1, bins = np.arange(min(eps1), max(eps1) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5,  label=f'Experimental ,{new_mean_lad}')
    plt.hist(eps2, bins = np.arange(min(eps2), max(eps2) + binwidth, binwidth), alpha=0.5,fill = False, histtype = 'step', edgecolor = 'g', linewidth = 1.5,  label=label_name + f", {new_mean_hm}")
    plt.hist(eps4, bins = np.arange(min(eps4), max(eps4) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'navy', linewidth = 1.5,  label='Thermodynamic Sensor'+ f", {new_mean_td_sens}")
    plt.vlines(x= mean_lad, ymin= d, ymax= 900, linestyles='dashed', linewidth = 0.6, colors='r')
    plt.vlines(x= mean_hm, ymin= d, ymax= 900, linestyles='dashed', linewidth = 0.6, colors='g')
    plt.vlines(x= mean_td_sens, ymin= d, ymax= 900, linestyles='dashed', linewidth = 0.6, colors='cornflowerblue')
    plt.xlim(-1,1)
    plt.subplots_adjust(wspace=0.26)
    #plt.
    """
    plt.hist(eps1, bins = np.arange(min(eps1), max(eps1) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5,  label=f'Experimental ,{new_mean_lad}')
    plt.hist(eps2, bins = np.arange(min(eps2), max(eps2) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'g', linewidth = 1.5,  label=label_name + f", {new_mean_hm}")
    plt.hist(eps3, bins = np.arange(min(eps3), max(eps3) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'b', linewidth = 1.5,  label=label__name+ f", {new_mean_td}")
    plt.hist(eps4, bins = np.arange(min(eps4), max(eps4) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'navy', linewidth = 1.5,  label='Thermodynamic Sensor'+ f", {new_mean_td_sens}")
    plt.hist(eps5, bins = np.arange(min(eps5), max(eps5) + binwidth, binwidth), alpha=0.5,fill = False, histtype = 'step', edgecolor = 'mediumaquamarine', linewidth = 1.5,  label='Hill Sensor'+ f", {new_mean_hm_sens}")
"""
"""
    plt.vlines(x= mean_lad, ymin= d, ymax= 2500, linestyles='dashed', linewidth = 0.6, colors='r')
    plt.vlines(x= mean_hm, ymin= d, ymax= 2500, linestyles='dashed', linewidth = 0.6, colors='g')
    plt.vlines(x= mean_td, ymin= d, ymax= 2500, linestyles='dashed', linewidth = 0.6, colors='b')
    plt.vlines(x= mean_td_sens, ymin= d, ymax= 2500, linestyles='dashed', linewidth = 0.6, colors='cornflowerblue')
    plt.vlines(x= mean_hm_sens, ymin= d, ymax= 2500, linestyles='dashed', linewidth = 0.6, colors='lime')
  """  
"""
    plt.legend(loc='upper right')
    plt.title('Distribution of Epistasis of \n'+strat + 'data', fontweight='bold')
    plt.xlabel('Epistasis', fontsize=10)
    plt.ylabel('Mutants', fontsize=10)
    plt.xlim(-1,1)
    """
    #%%
    # RSS
    # pvalue, statistic, sm, boxplot
    #%%
get_histograms(eps_lad, eps_hm, eps_td, eps_td_sens, eps_hm_sens, 1, 3, strat = ' sensor')
    # eps4 = td_stripe, eps5 = td_sens, eps6 = hm_stripe, eps7 = hm_sens
    #%%
    # RSS 
    
#%%
# testing lad epistasis values 
# def 
#df1 = pd.read_excel('../data/Source_Data.xlsx', 'Figure 2', usecols="J")
#df2 = df1.drop([0])
#df3 = df2.dropna()
#eps_rub = df3['Unnamed: 9']
def compare_df(df1, df2):
    new_df = pd.DataFrame(list(zip(df1, df2)),columns =['old', 'new'])
    new_df['different'] = np.where(abs((new_df['old']/new_df['new'])-1)>0.2, 'yes', '-')
    return new_df
def get_epsInfo(eps, df):
    # find row 
    # in that row access genotype category, genotype, inducer level
    # return dictionary
    i = np.isclose(df['Ep'], eps).argmax()
    output = []
    gen_cat = df.iloc[i]['genotype category']
    gen = df.iloc[i]['genotype']
    ind_lvl = df.iloc[i]['inducer level']
    output.append(gen_cat)
    output.append(gen)
    output.append(ind_lvl)
    return output
#%%
# to be used in Figure 3
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
       
def get_allBoxplots(data1, data2, data3):
   fig = plt.figure()
   plt.subplot(3,1,1)
   get_boxplot(data1)
   plt.subplot(3,1,2)
   get_boxplot(data2)
   plt.ylabel('Epistasis', fontsize=10)
   plt.subplot(3,1,3)
   get_boxplot(data3)
   plt.xlabel('Genotypes', fontsize=10)

#%%
fig, axs = plt.subplots(1,3, figsize=(15,5.5))
#ax = fig.subplots(1,3) 
#ax[0].figure(figsize=(8,8))
axs[0].hist(eps_lad, bins = np.arange(min(eps_lad), max(eps_lad) + 0.1, 0.1), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5,  label='Experimental, Mean:-0.198')
axs[0].hist(eps_hm, bins = np.arange(min(eps_hm), max(eps_hm) + 0.1, 0.1), alpha=0.5,fill = False, histtype = 'step', edgecolor = 'g', linewidth = 1.5,  label='Hill Model, Mean: -0.172')
axs[0].hist(eps_hm_sens, bins = np.arange(min(eps_hm_sens), max(eps_hm_sens) + 0.1, 0.1), alpha=0.5,fill = False, histtype = 'step', edgecolor = 'mediumaquamarine', linewidth = 1.5,  label='Hill Sensor, Mean:-0.04')
axs[0].vlines(x= -0.198, ymin= 0, ymax= 1000, linestyles='dashed', linewidth = 0.6, colors='r')
axs[0].vlines(x= -0.172, ymin= 0, ymax= 1000, linestyles='dashed', linewidth = 0.6, colors='g')
axs[0].vlines(x= -0.04, ymin= 0, ymax= 1000, linestyles='dashed', linewidth = 0.6, colors='lime')
axs[0].axis(xmin=-1,xmax=1)
axs[0].tick_params(labelsize=20)
axs[0].set_ylabel('Frequency', fontsize=24)
#fig.legend()
#
axs[1].hist(eps_lad, bins = np.arange(min(eps_lad), max(eps_lad) + 0.1, 0.1), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5)
axs[1].hist(eps_td, bins = np.arange(min(eps_td), max(eps_td) + 0.1, 0.1), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'b', linewidth = 1.5, label='Thermodynamic Model, Mean: 0.249')
axs[1].hist(eps_td_sens, bins = np.arange(min(eps_td_sens), max(eps_td_sens) + 0.1, 0.1), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'navy', linewidth = 1.5,  label='Thermodynamic Sensor, Mean:-0.093')
axs[1].vlines(x= -0.198, ymin= 0, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='r')
axs[1].vlines(x= 0.249, ymin= 0, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='b')
axs[1].vlines(x= -0.093, ymin= 0, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='cornflowerblue')
axs[1].set_xlabel('Epistasis', fontsize=10)
#axs[1].axis(xmin=-1,xmax=1)
axs[1].set_title('Distribution of Epistasis', fontweight='bold')
axs[1].title.set_fontsize(24)
axs[1].xaxis.label.set_fontsize(24)
axs[1].tick_params(labelsize=20)
fig.legend(bbox_to_anchor=(1.29,0.9),fontsize=18)
#
axs[2].hist(eps_lad, bins = np.arange(min(eps_lad), max(eps_lad) + 0.1, 0.1), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5,  label='Experimental, Mean:-0.198')
axs[2].hist(eps_hm, bins = np.arange(min(eps_hm), max(eps_hm) + 0.1, 0.1), alpha=0.5,fill = False, histtype = 'step', edgecolor = 'g', linewidth = 1.5,  label='Hill Model, Mean: -0.172')
axs[2].hist(eps_td_sens, bins = np.arange(min(eps_td_sens), max(eps_td_sens) + 0.1, 0.1), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'navy', linewidth = 1.5,  label='Thermodynamic Sensor, Mean:-0.093')
axs[2].vlines(x= -0.198, ymin= 0, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='r')
axs[2].vlines(x= -0.172, ymin= 0, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='b')
axs[2].vlines(x= -0.093, ymin= 0, ymax= 700, linestyles='dashed', linewidth = 0.6, colors='cornflowerblue')
axs[2].axis(xmin=-1,xmax=1)
axs[2].tick_params(labelsize=20)