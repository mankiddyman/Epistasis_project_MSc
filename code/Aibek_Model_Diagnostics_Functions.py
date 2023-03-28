#%%
import numpy as np
from scipy import stats
from Models import *
from Epistasis_calc_functions import *
import matplotlib.pyplot as plt
import statistics
from statistics import mean
import seaborn as sns
"""
Aibek Kappassov
March 21, 2023
"""
#%%
# 1
# retrieving epistasis data
def get_eData(model):
    if model == 'observed':
        eps_exp_df = pd.read_excel('../results/Eps_observed.xlsx')
        eps = eps_exp_df['Ep']   
    elif model == 'model_hill':
        eps_hill_df = pd.read_excel('../results/Eps_model_hill.xlsx')
        eps = eps_hill_df['Ep']
    elif model == 'model_thermodynamic':
        eps_td_df = pd.read_excel('../results/Eps_model_thermodynamic.xlsx')
        eps = eps_td_df['Ep']
    return eps 
#%%
def getEps(model_desc:str = "observed"):
    df =  pd.read_excel(f"../results/Eps_{model_desc}.xlsx").Ep
    return df
eps_lad = getEps("observed")
eps_hm = getEps('model_hill_all')
eps_td = getEps('model_thermodynamic_all')
eps_hm_sens = getEps('model_hill_sensor')
eps_td_sens = getEps('model_thermodynamic_sensor')

#%%
# 2
# epistasis disribution 
def get_histograms(eps_lad,eps_hm,eps_hm_sens,eps_td,eps_td_sens):
    #define paramaeters for the plot
    exp_col = 'dimgray'
    hill_all_col = 'springgreen'
    hill_sens_col = 'mediumseagreen'
    therm_all_col = 'turquoise'
    therm_sens_col = 'slateblue'
    mean_lad = round(eps_lad.mean(),3)
    mean_hm= round(eps_hm.mean(),3)
    mean_td= round(eps_td.mean(),3)
    mean_hm_sens= round(eps_hm_sens.mean(),3)
    mean_td_sens= round(eps_td_sens.mean(),3)
    hist_width = 2 #mean of hist outlines
    mean_width = 1.5 #width of mean lines
    sub_size = 17 #subtitle size

    fig, axs = plt.subplots(1,3, figsize=(15,5.5))
    y_obs, x, _ = axs[0].hist(eps_lad, bins = np.arange(min(eps_lad), max(eps_lad) + 0.1, 0.1), fill = False, histtype = 'step', edgecolor = exp_col, linewidth = hist_width,  label=f'Experimental, Mean: {mean_lad}')
    #histograms, hill
    y_hill_all, x , _ = axs[0].hist(eps_hm, bins = np.arange(min(eps_hm), max(eps_hm) + 0.1, 0.1), fill = False, histtype = 'step', edgecolor = hill_all_col, linewidth = hist_width,  label=f'Hill all. Mean: {mean_hm}')
    y_hill_sens, x, _ = axs[0].hist(eps_hm_sens, bins = np.arange(min(eps_hm_sens), max(eps_hm_sens) + 0.1, 0.1),fill = False, histtype = 'step', edgecolor = hill_sens_col, linewidth = hist_width,  label=f'Hill Sensor, Mean: {mean_hm_sens}')
    #mean lines, hill
    axs[0].axvline(x= mean_lad, ls='dashed', linewidth = mean_width, c=exp_col)
    axs[0].axvline(x= mean_hm,  ls='dashed', linewidth = mean_width, c=hill_all_col)
    axs[0].axvline(x= mean_hm_sens, ls='dashed', linewidth = mean_width, c=hill_sens_col)
    #axes, hill
    axs[0].axis(xmin=-1,xmax=1)
    axs[0].tick_params(labelsize=20)
    axs[0].set_ylabel('Frequency', fontsize=24)
    axs[0].set_title(f'Hill functions', fontsize = sub_size)
    #hist, therm
    axs[1].hist(eps_lad, bins = np.arange(min(eps_lad), max(eps_lad) + 0.1, 0.1), fill = False, histtype = 'step', edgecolor = exp_col, linewidth = hist_width)
    y_therm_all, x, _ = axs[1].hist(eps_td, bins = np.arange(min(eps_td), max(eps_td) + 0.1, 0.1), fill = False, histtype = 'step', edgecolor = therm_all_col, linewidth = hist_width, label=f'Thermodynamic All, Mean: {mean_td}')
    y_therm_sens, x, _ = axs[1].hist(eps_td_sens, bins = np.arange(min(eps_td_sens), max(eps_td_sens) + 0.1, 0.1), fill = False, histtype = 'step', edgecolor = therm_sens_col, linewidth = hist_width,  label='Thermodynamic Sensor, Mean:-0.093')
    #means, therm
    axs[1].axvline(x= mean_lad, ymin= 0, ymax= 700, ls='dashed', linewidth = mean_width, c=exp_col)
    axs[1].axvline(x= mean_td, ls='dashed', linewidth = mean_width, c=therm_all_col)
    axs[1].axvline(x= mean_td_sens, ls='dashed', linewidth = mean_width, c=therm_sens_col)
    #axes, therm
    axs[1].set_xlabel('Epistasis', fontsize=10)
    axs[1].xaxis.label.set_fontsize(24)
    axs[1].tick_params(labelsize=20)
    axs[1].set_title(f'Thermodynamic functions', fontsize = sub_size)
    #hist, top
    axs[2].hist(eps_lad, bins = np.arange(min(eps_lad), max(eps_lad) + 0.1, 0.1), fill = False, histtype = 'step', edgecolor = exp_col, linewidth = hist_width,  label=f'Experimental')
    axs[2].hist(eps_hm, bins = np.arange(min(eps_hm), max(eps_hm) + 0.1, 0.1),fill = False, histtype = 'step', edgecolor = hill_all_col, linewidth =hist_width,  label=f'Hill Model')
    axs[2].hist(eps_td_sens, bins = np.arange(min(eps_td_sens), max(eps_td_sens) + 0.1, 0.1), fill = False, histtype = 'step', edgecolor = therm_sens_col, linewidth = hist_width,  label=f'Thermodynamic Sensor')
    #means, top
    axs[2].axvline(x= mean_lad, ls='dashed', linewidth = mean_width, c=exp_col)
    axs[2].axvline(x= mean_hm,ls='dashed', linewidth = mean_width, c='b')
    axs[2].axvline(x= mean_td_sens,ls='dashed', linewidth = mean_width, c=therm_sens_col)
    #axes, top
    #axs[2].axis(xmin=-1,xmax=1)
    axs[2].tick_params(labelsize=20)
    axs[2].set_title(f'Top functions', fontsize = sub_size)
    #ylims
    ylim = max(y_obs.max(), y_hill_all.max(), y_hill_sens.max(), y_therm_all.max(), y_therm_sens.max())*1.1
    #axs[0].set_ylim(0, ylim)
    #axs[1].set_ylim(0, ylim)
    #axs[2].set_ylim(0, ylim)
    #legend
    axs[0].legend(loc = 'lower left',fontsize=18, bbox_to_anchor = (-0.5,-0.5), frameon=False)
    axs[1].legend(loc = 'lower left',fontsize=18, bbox_to_anchor = (-0.4,-0.4), frameon=False)
    #axs[2].legend(loc = 'lower left',fontsize=18, bbox_to_anchor = (0.05,-0.5))
    #title
    fig.suptitle('Distribution of Epistasis ', fontweight='bold', fontsize = 24, y = 1.05)
    #save to jpg
    fig.savefig('../results/Ep_compare_hist.jpg', bbox_inches='tight')
#%%
get_histograms(eps_lad, eps_hm, eps_hm_sens, eps_td, eps_td_sens)
    
#%%
# 3
# RSS 
def get_rss(df1, df2): # hm, td
    rss1 = df1['rss'].tolist()
    rss2 = df2['rss'].tolist()
    rss1.extend(rss2)
    mod = ['hill']*30 + ['thermodynamic']*30
    node = ['Output']*10 +['Regulator']*10 + ['Sensor']*10+['Output']*10 +['Regulator']*10 + ['Sensor']*10
    df = pd.DataFrame(list(zip(node, mod, rss1)),columns =['Node', 'Model', 'RSS'])
    plt.figure()
    vp = sns.violinplot(data=df, x="Node", y="RSS", hue="Model", split=True)
    plt.xlabel('Node',fontsize=11)
    plt.ylabel('RSS',fontsize=11)
    plt.title("Violin Plots for RSS", fontweight='bold',fontsize=12)
    
# data_hm = pd.read_excel('../results/model_hill.modelSM_params_all.xlsx')
# data_td = pd.read_excel('../results/model_thermodynamic.modelSM_params_sensor.xlsx')
# get_rss(data_hm,data_td) 
#%%
# 4
# two-sample t test
def get_pvalue(eps1, eps2, n):
    eps1_list = list(eps1)
    eps2_list = list(eps2)
    # observed epistasis
    eps1_test = np.array(eps1_list)
    # model epistasis
    eps2_test = np.array(eps2_list)
    # check equal variance assumption
    var_eps1 = np.var(eps1)
    var_eps2 = np.var(eps2)
    ratio_var = np.divide(var_eps1, var_eps2)
    # conduct two sample t test 
    alpha = 0.05
    option = n
    # change model name
    if ratio_var < 4:
        test = stats.ttest_ind(a=eps1_test, b=eps2_test, equal_var=True)
        p_value = getattr(test, 'pvalue')
        if option == 1:
            model_name = 'Hill model'
        elif option == 2:
            model_name = 'Hill Shaky model'
        elif option == 3:
            model_name = 'Thermodynamic model'
        elif option == 4:
            model_name == 'Competition Model'
        if p_value<alpha:
            message = "p value = " +str(p_value)+  " < 0.05"
            message = message + " .We reject the null hypothesis. We have sufficient evidence to conclude that the experimental and " + model_name + " epistasis means are not equal."
        else:
            message = "p value = " +str(p_value)+  " > 0.05"
            message = message + " .We fail to reject the null hypothesis. We do not have sufficient evidence to conclude that the experimental and " + model_name + " epistasis means are not equal."
    return message
# 5
def get_statistic(eps1, eps2):
    eps1_list = list(eps1)
    eps2_list = list(eps2)
    # observed epistasis
    eps1_test = np.array(eps1_list)
    # model epistasis
    eps2_test = np.array(eps2_list)
    test = stats.ttest_ind(a=eps1_test, b=eps2_test, equal_var=True)
    p_value = getattr(test, 'pvalue')
    p_value_str = str(p_value)
    word_list = p_value_str.split('e', 1)
    stat = float(word_list[0])
    power = int(word_list[1])
    power_str = str(power)
    stat1 = round(stat,2)
    pvalue = str(stat1) + 'e' + power_str
    line = 'p<' + pvalue
    return line  

#%%
# functions used/incorporated in chi_figure3_func.py and other plots/epistasis calculations
# 6, 7
# testing lad epistasis values 
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
# functions used/incorporated in chi_figure3_func.py and other plots/epistasis calculations
# 8, 9
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
# test functions above
# eps_exp_df = get_Eps('observed')
# eps = eps_exp_df['Ep']
# df1 = pd.read_excel('../data/Source_Data.xlsx', 'Figure 2', usecols="J")
# df2 = df1.drop([0])
# df3 = df2.dropna()
# eps_rub = df3['Unnamed: 9']
# pizda = compare_df(eps_rub,eps)

#%%
# variables used in model diagnostics
# eps_lad = get_eData('observed')
# eps_hm = get_eData('model_hill')
# eps_td = get_eData('model_thermodynamic')