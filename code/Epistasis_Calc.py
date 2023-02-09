import math
import pandas as pd
import numpy as np
from data_wrangling import meta_dict
from Models import *
import re
from itertools import permutations

#Calculate epistasis for all double and triple mutants from observed and expected flourcsence at low, medium and high inducer concs

df_DM = meta_dict['DM']
df_TM = meta_dict['TM']
df_M = (pd.concat([df_DM, df_TM], ignore_index = True)).drop_duplicates()

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
def G_hat_logadd(model, mutants:list):
    G_hat = g_hat(mutants[0], model)
    for mut in mutants[1:]:
        G_hat = np.multiply(g_hat(mut, model),G_hat)
    G_hat = np.log10(G_hat)
    return G_hat

#get list of possible mutant id names - e.g ['Sensor1','Regulator1'] --> ['S1_R1' , 'R1_S1']
def get_mut_ids(mutants:list):
    muts = []
    for mut in mutants:
        mut_index = re.search("[0-9]",mut).start()
        mutant1_id = f"{mut[:1]}{mut[mut_index:]}"
        muts.extend([mutant1_id])
    #permutations of mutants to search against in df_M
    mut_perms = list(permutations(muts))
    for i,ids in enumerate(mut_perms):
        string = ''
        for j in range(len(ids)):
            string += ids[j] +'_'
        mut_perms[i] = string[:-1]
    return mut_perms

#the opposite of get_mut_ids that turns id into list of names e.g. S1_O1 --> ['Sensor1','Regulator1']
def get_mut_names(mut_id:str):
    mutant_names = []
    ids = mut_id.split('_')
    for mut in ids:
        if mut.startswith('S'):
            mut = f"Sensor{mut[1:]}"
            mutant_names.extend([mut])
        if mut.startswith('R'):
            mut = f"Regulator{mut[1:]}"
            mutant_names.extend([mut])
        if mut.startswith('O'):
            mut = f"Output{mut[1:]}"
            mutant_names.extend([mut])
    return mutant_names

#observed GFP for double or triple mutants
#takes list of possible combinations of mutation ids as parameter e.g. ['O1_R1', 'R1_O1']
def G_obs(mutants:list):
    mut_perms = get_mut_ids(mutants)
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

    G_obs = np.log10(np.divide(MT_mean, WT_mean))
    return G_obs

def Epistasis(model, mutants:list):
    Epsilon = G_obs(mutants) - G_hat_logadd(model, mutants)
    return Epsilon

#calculate epistasis for all double mutants
def get_Eps(model=model_4_pred):
    Ep_low = []
    Ep_medium = []
    Ep_high = []
    for mut_id in df_M['genotype'][(df_M['inducer level'] == 'low') | (df_M['inducer level'] == 'low')][1:]:
        mut_names = get_mut_names(mut_id)
        Ep_low += [Epistasis(model_4_pred, mut_names)[0]]
        Ep_medium += [Epistasis(model_4_pred, mut_names)[1]]
        Ep_high += [Epistasis(model_4_pred, mut_names)[2]]
    Eps = [math.nan]*3 +Ep_low + Ep_medium + Ep_high
    return Eps

df_M["mean epistasis"] = get_Eps()