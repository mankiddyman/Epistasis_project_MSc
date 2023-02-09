import math
import pandas as pd
import numpy as np
from data_wrangling import meta_dict
from Models import *
import re
import itertools as iter

#Calculate epistasis for all double and triple mutants from observed and expected flourcsence at low, medium and high inducer concs

df_DM = meta_dict['DM']
df_TM = meta_dict['TM']
df_M = (pd.concat([df_DM, df_TM])).drop_duplicates()

df_fits = pd.read_excel('../data/SM_params.xlsx').rename(columns={'Unnamed: 0': 'mutant'})

I_conc = {'low':0.0, 'medium': 0.000195, 'high':0.2}
WT_params = df_fits.iloc[0].values.flatten().tolist()[1:-1]

#expected output flourecence for a single mutant predicted form a chosen model at concentrations in I_conc dictionary
def g_hat(mutant:str, model):
    I_conc_np = np.array(list(I_conc.values()))
    params = df_fits.loc[df_fits['mutant'] == mutant].values.flatten().tolist()[1:-1]
    g_hat = np.divide(model(I_conc_np, *params)[-1],model(I_conc_np, *WT_params)[-1])
    return g_hat

#log additive expected GFP for double or triple mutants for a given model
def G_hat(model, mutants:list):
    G_hat = g_hat(mutants[0], model)
    for mut in mutants[1:]:
        G_hat = np.multiply(g_hat(mut, model),G_hat)
    G_hat = np.log10(G_hat)
    return G_hat

#get list of possible mutant id names - e.g S1_R1 vs R1_S1
def get_mut_id(mutants:list):
    muts = []
    for i, mut in enumerate(mutants):
        mut_index = re.search("[0-9]",mut).start()
        mutant1_id = f"{mut[:1]}{mut[mut_index:]}"
        muts.extend([mutant1_id])
    #permutations of mutants to search against in df_TM
    mut_perms = (list(iter.permutations(muts)))
    for i,ids in enumerate(mut_perms):
        string = ''
        for j in range(len(ids)):
            string += ids[j] +'_'
        
        mut_perms[i] = string[:-1]
    return mut_perms

#observed GFP for double or triple mutants
def G_obs(mutants:list):
    mut_perms = get_mut_id(mutants)
    df_mutant = df_M.loc[df_M['genotype'] == mut_perms[0]]
    for mut in mut_perms[1:]:
        dfa = df_M.loc[df_M['genotype'] == mut]
        df_mutant = pd.concat([df_mutant, dfa])
    
    #get mean and std for low, med, high inducer concs
    #assumes low inducer concs come above medium, which come before high in mutant data
    MT_mean = np.array(df_mutant['obs_fluo_mean'])
    MT_sd = np.array(df_mutant['obs_SD'])

    #WT flourecence mean and sd
    WT_mean = np.array(df_M['obs_fluo_mean'].loc[df_M["genotype"]== "WT"])
    WT_sd = np.array(df_M['obs_SD'].loc[df_M["genotype"]== "WT"])

    G_obs = np.divide(MT_mean, WT_mean)
    return G_obs

def Epistasis(G_obs, G_hat):
    return G_obs - G_hat

