# the point of this script is to collect all our data in all its forms and output an object to enable easy interrogation.
#the structure will be a dictionary with 4 keys Wt, SM,DM, TM the values for these keys will be dataframes containing the appropriate data
import os
import math



meta_dict={"WT":None,"SM":None,"DM":None,"TM":None}


meta_dict["WT"]=pd.read_csv('../data/WT_single.csv')
#now single mutants
path_main = "../data/mutants_separated"

#generating the names of mutants
all_files=os.listdir(os.path.join(path_main,"mutants_separated"))
dat_files=[x for x in all_files if '.dat' in x]

SM_names=[x[:-4] for x in dat_files]


sm_df = pd.DataFrame(columns=["Inducer_conc","Mutant_ID",'Output_mean', 'Output_stdev', 'Regulator_mean',"Regulator_stdev","Sensor_mean","Sensor_stdev","Stripe_mean","Stripe_stdev"])
mean=[]
stdev=[]
Output_mean=[]
Output_stdev=[]
Output_genotype=[]
Regulator_mean=[]
Regulator_stdev=[]
Regulator_genotype=[]
Sensor_mean=[]
Sensor_stdev=[]
Sensor_genotype=[]

for mutant in SM_names:
    mean=[]
    stdev=[]
    mutant_table=pd.read_table(os.path.join(path_main,"mutants_separated",(mutant+".dat")),header=None)
    indexes=int(len(mutant_table)/3)
    for i in range(0,indexes):
        
        mean.append(mutant_table.iloc[[i,i+indexes,i+(indexes*2)],1].mean())
        stdev.append(mutant_table.iloc[[i,i+indexes,i+indexes*2],1].std())
    Genotype=(SM_names*3)
    Genotype.sort()
    if mutant.startswith("Output"):
        print(mutant)
        Output_mean+=mean
        
        Output_stdev.extend(stdev)
        Output_genotype.extend([mutant]*indexes)
    if mutant.startswith("Regulator"):
        Regulator_mean.extend(mean)
        Regulator_stdev.extend(stdev)
        Regulator_genotype.extend([mutant]*indexes)
        
    if mutant.startswith("Sensor"):
        Sensor_mean.extend(mean)
        Sensor_stdev.extend(stdev)
        Sensor_genotype.extend([mutant]*indexes)
        
#for output
print(len(Output_mean))
Output_mean+=[math.nan]*int(2*(len(SM_names)/3)*indexes)    
Output_stdev+=[math.nan]*int(2*(len(SM_names)/3)*indexes)
Output_genotype+=[math.nan]*int(2*(len(SM_names)/3)*indexes)
