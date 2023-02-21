#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from Models import *
from data_wrangling import *
from itertools import chain
from PyPDF2 import PdfMerger
import inspect

 #change dictionary to list of coeffs
def coef_dict_to_list(coef_dict):
    return list(coef_dict.values())
#the single mutant to be studied
   # function to input paramaters as a list

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

    #now gonna do the scipy.optimize.minimize
def WT_fit_plot(ax, y,params,label:str):
        return ax.plot(Signal, y, **params,label=f"{label}")


    #define scatter plotting function with log scales
def WT_Plotter(ax,y, params):
    out = ax.scatter(Signal, data_[y], **params, marker = 'o')
    xScale = ax.set_xscale('log')
    yScale = ax.set_yscale('log')

    return out, xScale, yScale 
def dict_to_list(params_dict,return_keys=False):
    if return_keys==True:
       a=[list(i.keys()) for i in list(params_dict.values())]
    elif return_keys==False:
        a=[list(i.values()) for i in list(params_dict.values())]

    return list(chain.from_iterable(a))
        

all_files=os.listdir(os.path.join(path_main,"mutants_separated"))
dat_files=[x for x in all_files if '.dat' in x]

SM_names=[x[:-4] for x in dat_files]
#SM_names=SM_names[0:10] #remove later

#paramater estimates from WT:
#WT_params = [6.59963463e+02, 1.63471394e+04, 1.25925588e+03, 1.16043953e+00, 1.99831036e+03, 2.04000859e+11, 2.77180774e+06, 8.37522575e-01, 5.47787795e-06, 6.71081447e+04, 1.41294256e-03, 5.41433758e+07,2.12643877e+00, 2.72060461e+00, 1.25044310e+00]

#WT_params = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]   

model_type_in=model_hill
params_dict={"sen_params":{"As":618.047086,"Bs":16278.856600,"Cs":1300.653790,"Ns":1.096541},"reg_params":{"Ar":1916.175610,"Br":18874.240800,"Cr":0.009030,"Nr":0.820433},"out_h_params":{"Ah":683.835638,"Bh":32464.380200,"Ch":0.000473},"out_params":{"Fo":2.821352,"Ao":0.632148,"Bo":0.972768 ,"Co":2.640174,"No":1.919339}}
SM_mutant_of_interest="Regulator1"
#%%
def model_fitting_SM(model_type=model_type_in,n_iter:float=1e5,mutant_range:slice=slice(0,len(SM_names))):
        #n_iter is how many iterations u want to evaluate
        #mutant_range is if u only want to do a specific few mutants at a time eg for testing
        
    df_param_compare = pd.DataFrame(columns=dict_to_list(params_dict,True), index = ['WT'])
    df_param_compare.iloc[0]=dict_to_list(params_dict)
    rss_initial=min_fun(params_list=dict_to_list(params_dict),data=meta_dict['WT'],model_type=model_type)
    df_param_compare['rss']=[rss_initial]

    for i in SM_names[mutant_range]: 
        SM_mutant_of_interest=i
        print("Fitting Mutant:",SM_mutant_of_interest)
        #use model chosen from Models.py to obtain fit of model to data

        df = pd.read_csv('../data/WT_single.csv')

        #substitute mutation data into the WT dataframe for a particular set of mutations
        #for single mutations to start:
        

        #example of above function to get SM dat for Sensor1 mutation:
        S1_df = get_data_SM(SM_mutant_of_interest)

      
        data_ = S1_df

        
        WT_params=dict_to_list(params_dict)
        print("starting paramaters:\n", params_dict)
        
        if SM_mutant_of_interest.startswith("Sensor"):
            bnds = ((0, None), (0, None), (0, None), (0, None), (WT_params[4], WT_params[4]), (WT_params[5],WT_params[5]), (WT_params[6], WT_params[6]),(WT_params[7], WT_params[7]), (WT_params[8],WT_params[8]), (WT_params[9],WT_params[9]), (WT_params[10],WT_params[10]), (0,None), (WT_params[12],WT_params[12]), (WT_params[13],WT_params[13]), (WT_params[14],WT_params[14]),(WT_params[15],WT_params[15]))
        elif SM_mutant_of_interest.startswith("Regulator"):
            bnds = ((WT_params[0],WT_params[0]), (WT_params[1],WT_params[1]), (WT_params[2],WT_params[2]) ,(WT_params[3],WT_params[3]),(0, None) ,(0, None),(0, None),(0, None), (WT_params[8],WT_params[8]), (WT_params[9],WT_params[9]), (WT_params[10],WT_params[10]), (0,None), (WT_params[12],WT_params[12]), (WT_params[13],WT_params[13]), (WT_params[14],WT_params[14]),(WT_params[15],WT_params[15]))
        elif SM_mutant_of_interest.startswith("Output"):
            bnds = ((WT_params[0],WT_params[0]), (WT_params[1],WT_params[1]), (WT_params[2],WT_params[2]) ,(WT_params[3],WT_params[3]),(WT_params[4], WT_params[4]), (WT_params[5],WT_params[5]), (WT_params[6], WT_params[6]),(WT_params[7], WT_params[7]), (0, None),(0, None),(0, None),(0, None),(0, None),(0, None),(0, None),(0, None))

        min_result=minimize(min_fun,args=(data_,model_type),x0= WT_params ,method='Nelder-Mead',bounds=bnds,options={"maxiter":n_iter,"disp":True})
        print("finished fitting")

        #plotting the predictions now
        #generating estimates
        Sensor_est_array,Regulator_est_array,Output_est_array, Stripe_est_array = model_type(params_list=min_result.x,I_conc=data_.S)

        Sensor_est_array_initial,Regulator_est_array_initial,Output_est_array_initial, Stripe_est_array_initial = model_type(I_conc=data_.S,params_list=WT_params)

        
        #now need to calculate residual sum of logsquares
        #defined as value of minimisation function for converged parameter set
        rss_converged=min_fun(data=data_,params_list=min_result.x,model_type=model_type)
        rss_initial=min_fun(data=data_,params_list=WT_params,model_type=model_type)
        rss_relative=rss_initial/(rss_converged+rss_initial) #% of rss explained by convergence
        
        df_parameters=pd.DataFrame({"epsilon":min_result.x-WT_params,"Initial_guesses":WT_params,"Converged":min_result.x})
        df_parameters.index=dict_to_list(params_dict,return_keys=True)

        
        #add mutant parameter values to df_param_compare to compare
        # a = min_result.x
        # a = np.append(a, rss_converged)
        # dfa = pd.DataFrame(a, index =list(df_param_compare)).T
        # dfa = dfa.set_index(pd.Series([SM_mutant_of_interest]))
        # df_param_compare = pd.concat([df_param_compare, dfa])

        Signal=data_.S
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
        title = ["SM data type data plots for mutation", SM_mutant_of_interest,"using model:",str(model_type.__name__)]
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

        fig.savefig(os.path.join("..","results",f"{SM_mutant_of_interest}"+"_"+str(model_type.__name__)+".pdf"), bbox_inches='tight')
        print("done fitting, now exporting")

        
        a=list(df_parameters.transpose().iloc[2]) #collecting converged parameters for export
        a.append(rss_converged)
        df_param_compare.loc[f'{SM_mutant_of_interest}']=a
        print(i,df_param_compare)
        # %%
    df_param_compare.to_excel(os.path.join("..","data",f"{model_type.__name__}"+"SM_params.xlsx"))
       

    pdfs= [os.path.join("..","results",(s+"_"+str(model_type.__name__)+".pdf")) for s in SM_names]
    merger = PdfMerger()

    for pdf in pdfs:
        merger.append(pdf)

    merger.write(os.path.join("..","results","SM_models_"+str(model_type.__name__+".pdf")))
    merger.close()
    print("done with all")
    return 1
#thermodynamic
params_dict={"sen_params":{"a_s":1,"Kp_s":1,"P":1,"Ki_s":1,"C_pi_s":1},"reg_params":{"Ar":1,"Br":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1},"out_params":{"Fo":1,"Ao":1,"Bo":1,"Co":1,"No":1}}


WT_params=dict_to_list(params_dict)
WT_params=list(np.array(WT_params)+100)
WT_params=[1]*13
WT_params=list(np.array(WT_params)+100)
bnds = ((0, None), (0, None), (0, None), (0, None), (0, None), (0, None), (0, None),(0, None), (0, None), (0, None),(1, None), (0,1), (0,1))
def get_WT_params(model_type=model_type_in,n_iter:float=1e5):
    #this function will estimate the wild type parameters for a given model.
    #now loading wt dataframe
    data_=meta_dict['WT']
    min_result=minimize(min_fun,args=(data_,model_type),x0=WT_params,method='Nelder-Mead',bounds=bnds,options={"maxiter":n_iter,"disp":True})

    Sensor_est_array,Regulator_est_array,Output_est_array, Stripe_est_array = model_type(params_list=min_result.x,I_conc=data_.S)
    
    Sensor_est_array_initial,Regulator_est_array_initial,Output_est_array_initial, Stripe_est_array_initial = model_type(I_conc=data_.S,params_list=WT_params)

    rss_converged=min_fun(data=data_,params_list=min_result.x,model_type=model_type)
    rss_initial=min_fun(data=data_,params_list=WT_params,model_type=model_type)
    rss_relative=rss_initial/(rss_converged+rss_initial) #% of rss explained by 
  
    fig, ((Sensor, Regulator), (Output, Stripe)) = plt.subplots(2,2, constrained_layout=True)
    
    WT_Plotter(Sensor, "Sensor", {'color':'red'})
    WT_fit_plot(Sensor, Sensor_est_array_initial, {'color':'black'},label="Initial Guess")
    WT_fit_plot(Sensor, Sensor_est_array,{'color':'red'},label="Converged")
    Sensor.set_title(r'inducer -> sensor (GFP output)')
    Sensor.set_xscale('log')
    Sensor.set_yscale('log')
    
    
    WT_Plotter(Regulator, "Regulator", {'color': 'blue'})
    WT_fit_plot(Regulator, Regulator_est_array, {'color':'blue'},label="Converged")
    WT_fit_plot(Regulator, Regulator_est_array_initial, {'color':'black'},label="")
    Regulator.set_title(r'inducer ->S -|R (GFP output)')
    Regulator.set_xscale('log')
    Regulator.set_yscale('log')
    

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
    title = ["WT data plots ","using model:",str(model_type.__name__)]
    df_parameters=pd.DataFrame({"epsilon":min_result.x-WT_params,"Initial_guesses":WT_params,"Converged":min_result.x})
    #df_parameters.index=dict_to_list(params_dict,return_keys=True)

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

    print("finished fitting")

#dictionary structure for input to model
#to add, feature that will generate initial guess from wildtype
#which u then feed to the actual running of the function

#fitting

#thermodynamic model
#rewrite this
params_dict={"sen_params":{"a_s":1,"Kp_s":1,"P":1,"Ki_s":1,"C_pi_s":1},"reg_params":{"a_r":1,"Kp_r":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1}}


params_dict={"sen_params":{"As":1,"Bs":1,"Cs":1,"Ns":1},"reg_params":{"Ar":1,"Br":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1},"out_params":{"Fo":1,"Ao":1,"Bo":1,"Co":1,"No":1}}

params_dict["sen_params"]


#now generating initial guesses

model_type_in=thermodynamic_model



#functionality to add ,
#  automatically find out number of arguements, automatically set bounds on paramaeters