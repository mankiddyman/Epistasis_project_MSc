import numpy as np
from scipy import stats
from Models import *
from Epistasis_calc_functions import *
from Model_Diagnostics_Functions import *
import matplotlib.pyplot as plt

#part of model_diagnostic_functions
# eps_obs = get_Data('observed')
# eps_hill = get_Data('model_hill')

# Observed_Eps_df = get_Eps('observed') #no need to run if eps already generated
# Hill_Eps_df = get_Eps('model_hill')

#%%
def Sort_mutants(model): #observed, model_hill, and thermodynamic
    if model == str('observed'):
        Eps_df = pd.read_excel('../results/Eps_observed.xlsx') #takes from results.
    elif model == str('model_hill'):
        Eps_df = get_Eps('model_hill') #saves time if already generated, just modift the df.
    elif model == str('thermodynamic'):
        Eps_df = pd.read_excel('../results/Eps_model_thermodynamic.xlsx')
    else:
        print('the inputted model is invalid, please selec from observed, model_hill or thermodynamic')

    #output
        #filter out O10
    o10 = Eps_df[Eps_df['genotype'].str.contains(str('O10'))] 
        #filter out 01 by removing duplicates from o10, merge two dataframes and assign boolean column for duplicated then filter out.
    o1 = Eps_df[Eps_df['genotype'].str.contains(str('O1'))] 
    merge = pd.concat([o1,o10], axis=0)
    merge['dups'] = merge.duplicated(keep=False)
    o1 = merge[(merge['dups']==False)]

    number = [2,3,4,5,6,7,8,9] #Filter out rest of mutant genotypes
    Output_df = pd.DataFrame()
    Output_df['O1'] = o1['Ep']
    Output_df.reset_index(drop=True, inplace=True) #gets rid of index from original table
    for num in number:
        genotype = 'O'+ str(num)
        temp = Eps_df[Eps_df['genotype'].str.contains(genotype)] #could simply with .query for genotypes
        temp.reset_index(drop=True, inplace=True)
        Output_df[genotype] = temp['Ep']
    o10.reset_index(drop=True, inplace=True)
    Output_df['O10'] = o10['Ep']

    #Sensor
    s10 = Eps_df[Eps_df['genotype'].str.contains(str('S10'))] 
    s1 = Eps_df[Eps_df['genotype'].str.contains(str('S1'))] 
    mergeS = pd.concat([s1,s10], axis=0) #merge dataframes
    mergeS['dups'] = mergeS.duplicated(keep=False)
    s1 = mergeS[(mergeS['dups']==False)]

    Sensor_df = pd.DataFrame()
    Sensor_df['S1'] = s1['Ep']
    Sensor_df.reset_index(drop=True, inplace=True)
    for num in number:
        genotype = 'S'+ str(num)
        temp = Eps_df[Eps_df['genotype'].str.contains(genotype)]
        temp.reset_index(drop=True, inplace=True)
        Sensor_df[genotype] = temp['Ep']
    s10.reset_index(drop=True, inplace=True)
    Sensor_df['S10'] = s10['Ep']

    #Regulator
    r10 = Eps_df[Eps_df['genotype'].str.contains(str('R10'))] 
    r1 = Eps_df[Eps_df['genotype'].str.contains(str('R1'))] 
    mergeR = pd.concat([r1,r10], axis=0)
    mergeR['dups'] = mergeR.duplicated(keep=False)
    r1 = mergeR[(mergeR['dups']==False)]

    #Filter out rest of mutant genotypes
    Regulator_df = pd.DataFrame()
    Regulator_df['R1'] = r1['Ep']
    Regulator_df.reset_index(drop=True, inplace=True)
    for num in number:
        genotype = 'R'+ str(num)
        temp = Eps_df[Eps_df['genotype'].str.contains(genotype)]
        temp.reset_index(drop=True, inplace=True)
        Regulator_df[genotype] = temp['Ep']
    r10.reset_index(drop=True, inplace=True)
    Regulator_df['R10'] = r10['Ep']

    return Output_df, Regulator_df, Sensor_df 


# %%

# def Epistasis_hist(Output_epi_df,Regulator_epi_df,Sensor_ep_df):
#histograms
