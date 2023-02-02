# the point of this script is to collect all our data in all its forms and output an object to enable easy interrogation.
#the structure will be a dictionary with 4 keys Wt, SM,DM, TM the values for these keys will be dataframes containing the appropriate data
import os
import math
import pandas as pd
import numpy as np
import re


meta_dict={"WT":None,"SM":None,"DM":None,"TM":None}


meta_dict["WT"]=pd.read_csv('../data/WT_single.csv')
#now single mutants (SM)
path_main = "../data/mutants_separated"

#and the SM stripe data
Stripe_data = pd.ExcelFile('../data/Source_Data.xlsx')
Stripe_data.sheet_names
SM_stripes = pd.read_excel(Stripe_data, 'Figure 1')

#to get list of means for 'sensor 1':
#SM_stripes['Sensor 1'][1:17]
#get column number of 'Sensor 1' column
#list(SM_stripes.columns).index('Sensor 1')

#generating the names of mutants
all_files=os.listdir(os.path.join(path_main,"mutants_separated"))
dat_files=[x for x in all_files if '.dat' in x]

SM_names=[x[:-4] for x in dat_files]

#%%
mean=[]
stdev=[]
Output_mean=[]
Output_stdev=[]
Regulator_mean=[]
Regulator_stdev=[]
Sensor_mean=[]
Sensor_stdev=[]
inducer = []
genotype = []
Stripe_mean = []
Stripe_stdev = []

for mutant in SM_names:
    mean=[]
    stdev=[]

    a = re.search("[0-9]",mutant).start() #index of mutant identifier
    mutant_stripe = f"{mutant[:a]} {mutant[a:]}"
    
    mutant_table=pd.read_table(os.path.join(path_main,"mutants_separated",(mutant+".dat")),header=None)
    indexes=int(len(mutant_table)/3)
    ind = mutant_table.iloc[0:indexes,0] #concentrations
    genotype += [mutant]*indexes
    inducer.extend(ind)
    for i in range(0,indexes):       
        mean.append(mutant_table.iloc[[i,i+indexes,i+(indexes*2)],1].mean())
        stdev.append(mutant_table.iloc[[i,i+indexes,i+indexes*2],1].std())

    Stripe_mean += list(SM_stripes[mutant_stripe][1:indexes+1])
    std_col = list(SM_stripes.columns).index(mutant_stripe)+1
    Stripe_stdev += list(SM_stripes.iloc[1:indexes+1, std_col])

    if mutant.startswith("Output"):

        Output_mean.extend(mean)        
        Output_stdev.extend(stdev)

        Regulator_mean.extend([math.nan]*len(mean))
        Regulator_stdev.extend([math.nan]*len(mean))
        Sensor_mean.extend([math.nan]*len(mean))
        Sensor_stdev.extend([math.nan]*len(mean))

    if mutant.startswith("Regulator"):

        Regulator_mean.extend(mean)
        Regulator_stdev.extend(stdev)

        Output_mean.extend([math.nan]*len(mean))
        Output_stdev.extend([math.nan]*len(mean))
        Sensor_mean.extend([math.nan]*len(mean))
        Sensor_stdev.extend([math.nan]*len(mean))
        
    if mutant.startswith("Sensor"):

        Sensor_mean.extend(mean)
        Sensor_stdev.extend(stdev)

        Output_mean.extend([math.nan]*len(mean))
        Output_stdev.extend([math.nan]*len(mean))
        Regulator_mean.extend([math.nan]*len(mean))
        Regulator_stdev.extend([math.nan]*len(mean))
#%%
#There was a missing value for sensor 7 at inducer conc 0.00020, this was entered manually as nan in Sensor7.dat


sm_df = pd.DataFrame({"Inducer" :  inducer,"Mutant_ID": genotype,'Output_mean': Output_mean, 'Output_stdev': Output_stdev, 'Regulator_mean': Regulator_mean,"Regulator_stdev": Regulator_stdev,"Sensor_mean": Sensor_mean,"Sensor_stdev": Sensor_stdev,"Stripe_mean": Stripe_mean,"Stripe_stdev": Stripe_stdev})

meta_dict["SM"] = sm_df

#%%