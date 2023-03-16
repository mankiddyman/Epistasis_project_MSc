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
eps_lad = get_eData('observed')
eps_hm = get_eData('model_hill')
eps_td = get_eData('model_thermodynamic')

#%%
def get_histograms(eps1, eps2, eps3, n, b): # x, y
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
    plt.hist(eps1, bins = np.arange(min(eps1), max(eps1) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'r', linewidth = 1.5,  label='Log-Additive Model')
    plt.hist(eps2, bins = np.arange(min(eps2), max(eps2) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'g', linewidth = 1.5,  label=label_name)
    plt.hist(eps3, bins = np.arange(min(eps3), max(eps3) + binwidth, binwidth), alpha=0.5, fill = False, histtype = 'step', edgecolor = 'b', linewidth = 1.5,  label=label__name)
    
    d, u = plt.ylim()
    mean_lad = statistics.mean(eps_lad)
    mean_hm = statistics.mean(eps_hm)
    mean_td = statistics.mean(eps_td)
    mean_lad1 = round(mean_lad,3)
    mean_hm1 = round(mean_hm,3)
    mean_td1 = round(mean_td,3)
    new_mean_lad = 'Mean: ' + str(mean_lad1)
    new_mean_hm = 'Mean: ' + str(mean_hm1)
    new_mean_td = 'Mean: ' + str(mean_td1)
    
    plt.vlines(x= mean_lad, ymin= d, ymax= 610, linestyles='dashed', linewidth = 0.6, colors='r')
    plt.vlines(x= mean_hm, ymin= d, ymax= 880, linestyles='dashed', linewidth = 0.6, colors='g')
    plt.vlines(x= mean_td, ymin= d, ymax= 830, linestyles='dashed', linewidth = 0.6, colors='b')
    
    plt.text(-0.43, 615, new_mean_lad, fontsize = 7)
    plt.text(-0.16, 875, new_mean_hm, fontsize = 7)
    plt.text(-0.12, 820, new_mean_td, fontsize = 7)
    
    plt.legend(loc='upper left')
    plt.title('Distribution of Epistasis', fontweight='bold')
    plt.xlabel('Epistasis', fontsize=10)
    plt.ylabel('Mutants', fontsize=10)
    #%%
    # RSS
    # pvalue, statistic, sm, boxplot
    #%%
    get_histograms(eps_lad, eps_hm, eps_td, 1, 3)
    
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