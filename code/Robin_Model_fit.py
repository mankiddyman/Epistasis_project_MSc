#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from Models import *
from data_wrangling import *


 #change dictionary to list of coeffs
def coef_dict_to_list(coef_dict):
    return list(coef_dict.values())
#the single mutant to be studied
   # function to input paramaters as a list
def coeff_input(param_list):
    if len(param_list) != len(coeffs):
        print('input paramaters list is not the same length as coefficents: COEFFIECIENTS UNCHANGED')
    else:
        for i,key in enumerate(list(coeffs.keys())):
            coeffs[str(key)] = param_list[i]
    return coeffs
def get_data_SM(mutation:str):
        df_MT = df
        data = meta_dict["SM"]
        data = data.loc[data['Mutant_ID'] == mutation]
        #WT data missing measurements for inducer = 0.2 so drop last columnn
        data = data[:-1]

        a = re.search("[0-9]",mutation).start()
        mut_col = f"{mutation[:a]}"
        mutation_means = f"{mutation[:a]}_mean"
        df_MT[[mut_col, "Stripe"]] = data[[mutation_means, "Stripe_mean"]].values
        df_MT[["Mutant_ID"]]=mutation
        if mutation.startswith("Output"):
            df_MT["Sensor"]=meta_dict["WT"].Sensor

        return df_MT
#data=pd.read_csv('../data/WT_single.csv')
def min_fun(coefficients,data):
        log_sen = np.log10(data.Sensor)
        log_reg = np.log10(data.Regulator)
        log_out = np.log10(data.Output)
        log_stripe = np.log10(data.Stripe)   
        ind = data.S
        Sensor_est, Regulator_est, Output_est, Stripe_est = model_4_pred(ind,*coefficients)
        log_sen_est = np.log10(Sensor_est)
        log_reg_est = np.log10(Regulator_est)
        log_out_est = np.log10(Output_est)
        log_stripe_est = np.log10(Stripe_est)

        #need to know what variables exist for each given mutant
        if "Mutant_ID" in data:
            mutant_id=data.Mutant_ID[0]
            if mutant_id.startswith("Sensor"):
                #need to ignore reg and output in fitting
                log_reg,log_reg_est,log_out,log_out_est=0,0,0,0
            # #if mutant_id.startswith("Regulator"):
            #     #need to ignore reg and output in fitting
            #     log_sen,log_sen_est,log_out,log_out_est=0,0,0,0
            # #if mutant_id.startswith("Output"):
            #     #need to ignore reg and sensor in fitting
            #     log_reg,log_reg_est,log_sen,log_sen_est=0,0,0,0
                
        

           

        result = np.power((log_sen - log_sen_est), 2)
        result += np.power((log_reg - log_reg_est), 2)
        result += np.power((log_out - log_out_est), 2)
        result += np.power((log_stripe - log_stripe_est), 2)
        return np.sum(result)
    #now gonna do the scipy.optimize.minimize
def WT_fit_plot(ax, y,params,label:str):
        return ax.plot(Signal, y, **params,label=f"{label}")


    #define scatter plotting function with log scales
def WT_Plotter(ax,y, params):
    out = ax.scatter(Signal, data_[y], **params, marker = 'o')
    xScale = ax.set_xscale('log')
    yScale = ax.set_yscale('log')

    return out, xScale, yScale 

all_files=os.listdir(os.path.join(path_main,"mutants_separated"))
dat_files=[x for x in all_files if '.dat' in x]

SM_names=[x[:-4] for x in dat_files]
#SM_names=SM_names[0:10] #remove later

#paramater estimates from WT:
WT_params = [6.08397103e+02, 1.52504577e+04, 1.66805905e+03, 1.19893355e+00, 6.87964693e+02, 2.34976114e+04, 6.23671728e-02, 3.91730917e-01, 5.90606548e+02, 3.52871257e+04, 5.29890033e-04, 8.29829791e-01, 4.28817019e+00, 3.13322189e+00, 1.80901848e+00]   
coeffs = {'A_s': 1,'B_s': 1,'C_s':1,'N_s':0, 'A_r':1,'B_r':1,'C_r':1,'N_r':0,'A_h':1, 'B_h':1, 'C_h':1, 'A_o':1,'B_o':1,'C_o': 1,'N_o': 0}
init_coeffs = coeff_input(WT_params)

#df to compare SM parameter estimates with WT
df_param_compare = pd.DataFrame(init_coeffs, index = ['WT'])
rss_initial=min_fun(data=meta_dict['WT'],coefficients=WT_params)
df_param_compare['rss']=[rss_initial]

#%%
for i in SM_names[0:4]:
    SM_mutant_of_interest=i
    print("Fitting Mutant:",SM_mutant_of_interest)
    #use model chosen from Models.py to obtain fit of model to data

    df = pd.read_csv('../data/WT_single.csv')

    #substitute mutation data into the WT dataframe for a particular set of mutations
    #for single mutations to start:
    

    #example of above function to get SM dat for Sensor1 mutation:
    S1_df = get_data_SM(SM_mutant_of_interest)

    # minimizing sum of sqrd deviations of log model and log data
    data_ = S1_df

    init_coeffs = coeff_input(WT_params)

    print("starting paramaters:", WT_params)
    bnds = ((0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None))
    min_result=minimize(min_fun,args=(data_),x0= WT_params ,method='Nelder-Mead',bounds=bnds,options={"maxiter":1e5,"disp":True})
    print("finished fitting")

    #plotting the predictions now
    #generating estimates
    Sensor_est_array,Regulator_est_array,Output_est_array, Stripe_est_array = model_4_pred(data_.S,*min_result.x)

    Sensor_est_array_initial,Regulator_est_array_initial,Output_est_array_initial, Stripe_est_array_initial = model_4_pred(data_.S,*WT_params)

    Signal=data_.S

    #now need to calculate residual sum of logsquares
    #defined as value of minimisation function for converged parameter set
    rss_converged=min_fun(data=data_,coefficients=min_result.x)
    rss_initial=min_fun(data=data_,coefficients=WT_params)
    rss_relative=rss_initial/(rss_converged+rss_initial) #% of rss explained by convergence
    WT_params = [6.08397103e+02, 1.52504577e+04, 1.66805905e+03, 1.19893355e+00, 6.87964693e+02, 2.34976114e+04, 6.23671728e-02, 3.91730917e-01, 5.90606548e+02, 3.52871257e+04, 5.29890033e-04, 8.29829791e-01, 4.28817019e+00, 3.13322189e+00, 1.80901848e+00]

    df_parameters=pd.DataFrame({"epsilon":min_result.x-WT_params,"Initial_guesses":WT_params,"Converged":min_result.x})
    df_parameters.index=coeffs.keys()

    
    #add mutant parameter values to df_param_compare to compare
    a = min_result.x
    a = np.append(a, rss_converged)
    dfa = pd.DataFrame(a, index =list(df_param_compare)).T
    dfa = dfa.set_index(pd.Series([SM_mutant_of_interest]))
    df_param_compare = pd.concat([df_param_compare, dfa])


    fig, ((Sensor, Regulator), (Output, Stripe)) = plt.subplots(2,2, constrained_layout=True)
    
    WT_Plotter(Sensor, "Sensor", {'color':'red'})
    WT_fit_plot(Sensor, Sensor_est_array_initial, {'color':'black'},label="Initial Guess")
    WT_fit_plot(Sensor, Sensor_est_array,{'color':'red'},label="Converged")
    Sensor.set_title(r'inducer -> sensor (GFP output)')
    Sensor.set_xscale('log')
    Sensor.set_yscale('log')
    
    if SM_mutant_of_interest.startswith("Sensor")!=True:
        WT_Plotter(Regulator, "Regulator", {'color': 'blue'})
    WT_fit_plot(Regulator, Regulator_est_array, {'color':'blue'},label="Converged")
    WT_fit_plot(Regulator, Regulator_est_array_initial, {'color':'black'},label="")
    Regulator.set_title(r'inducer ->S -|R (GFP output)')
    Regulator.set_xscale('log')
    Regulator.set_yscale('log')
    

    if SM_mutant_of_interest.startswith("Sensor")!=True:
         WT_Plotter(Output, "Output", {'color': 'purple'})
    WT_fit_plot(Output, Output_est_array, {'color':'purple'},label="Converged")
    WT_fit_plot(Output, Output_est_array_initial, {'color':'black'},label="")
    Output.set_title(r'inducer -> S -| Output (GFP)')
    Output.set_xscale('log')
    Output.set_yscale('log')
    

    WT_Plotter(Stripe,"Stripe", {'color': 'green'})
    WT_fit_plot(Stripe, Stripe_est_array, {'color':'green'},label="Converged")
    WT_fit_plot(Stripe, Stripe_est_array_initial, {'color':'black'},label="")
    Stripe.set_title(r'Full circuit with stripe')
    fig.legend(bbox_to_anchor=(1.3, 1))
    title = ["SM data type data plots for mutation", SM_mutant_of_interest]
    txt=f''' Across all four plots: \n
    RSS (converged)={round(rss_converged,3)} \n
    RSS (initial)={round(rss_initial,3)}\n
    RSS (% reduction)={round(rss_relative,3)}\n
    '''
    txt+=str(df_parameters)
    fig.text(.11,-.81,txt,wrap=True)
    txt=f'''{min_result}
    '''
    fig.text(0.7,-.81,txt,wrap=True)


    plt.suptitle(title,fontweight="bold")

    plt.show()
    print("final parameter estimates:", min_result.x)

    fig.savefig(os.path.join("..","results",f"{SM_mutant_of_interest}"+".pdf"), bbox_inches='tight')

df_param_compare.to_excel('../data/SM_params.xlsx')
    # %%
from PyPDF2 import PdfMerger

pdfs= [s+".pdf" for s in SM_names]
merger = PdfMerger()

for pdf in pdfs:
    merger.append(pdf)

merger.write(os.path.join("..","results","SM_models.pdf"))
merger.close()