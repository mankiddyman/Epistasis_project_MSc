import numpy as np
from scipy.integrate import odeint
import pandas as pd

# plots make sense
#should investigate the relationships between inducer and sensor (calculated using regression of inducer )
#function to do subsequent non linear regression is :
#%%
#we will now create model functions for each of the steady state values of Sensor, Regulator and Output.
#Henceforth referred to as S, R and O

#hill function model
def init_model(I_conc,A_s,B_s,C_s,N_s, A_r,C_r,N_r,A_o,C_o,N_o):
    Sensor= (A_s+B_s*(C_s*I_conc)**N_s)/(1+np.power(C_s*I_conc,N_s))
    Regulator = (A_r)/(1+ np.power(C_r*Sensor,N_r))
    Sens_Reg = Sensor + Regulator
    Output = (A_o)/(1+np.power(C_o*Sensor,N_o))
    Stripe = (A_o)/(1+np.power(C_o*Sens_Reg,N_o))
    return Sensor, Regulator, Output, Stripe

def model_leaky(I_conc,A_s,B_s,C_s,N_s,L_r,A_r,C_r,N_r,L_o,A_o,C_o,N_o):
    Sensor= (A_s+B_s*(C_s*I_conc)**N_s)/(1+np.power(C_s*I_conc,N_s))
    Regulator = L_r+(A_r)/(1+ np.power(C_r*Sensor,N_r))
    Sens_Reg = Sensor + Regulator
    Output = (A_o)/(1+np.power(C_o*Sensor,N_o))
    Stripe = L_o+(A_o)/(1+np.power(C_o*Sens_Reg,N_o))
    return Sensor, Regulator, Output, Stripe

#example input to model_hill
#assumes degredation of lacI, TetR and GFP are constant
#params_dict={"sen_params":{"As":1,"Bs":1,"Cs":1,"Ns":1},"reg_params":{"Ar":1,"Br":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1},"out_params":{"Fo":1,"Ao":1,"Bo":1,"Co":1,"No":1}}
#%%
class model_hill:
    def __init__(self,params_list:list,I_conc):
        self.params_list=params_list
        self.I_conc=I_conc
        self.example_dict={"sen_params":{"A_s":1,"B_s":1,"C_s":1,"N_s":1},"reg_params":{"A_r":1,"B_r":1,"C_r":1,"N_r":1},"out_h_params":{"A_h":1,"B_h":1,"C_h":1},"out_params":{"A_o":1,"B_o":1,"C_o":1,"N_o":1},"free_params":{"F_o":1}}
        self.n_parameters=16
    @staticmethod    
    def model(params_list,I_conc):
        correct_length=16
        #S is subscript for parameters corresponding to Sensor
        #R is subscript for parameters corresponding to Regulator
        #H is subscript for parameters corresponding to the half network I->S -| O
        #O is subscript for parameters corresponding to Output
        # As=params_dict['sen_params']['As']
        # Bs=params_dict['sen_params']['Bs']
        # Cs=params_dict['sen_params']['Cs']
        # Ns=params_dict['sen_params']['Ns']
        # Ar=params_dict['reg_params']['Ar']
        # Br=params_dict['reg_params']['Br']
        # Cr=params_dict['reg_params']['Cr']
        # Nr=params_dict['reg_params']['Nr']
        # Ah=params_dict['out_h_params']['Ah']
        # Bh=params_dict['out_h_params']['Bh']
        # Ch=params_dict['out_h_params']['Ch']
        # Fo=params_dict['out_params']['Fo']
        # Ao=params_dict['out_params']['Ao']
        # Bo=params_dict['out_params']['Bo']
        # Co=params_dict['out_params']['Co']
        # No=params_dict['out_params']['No']
        
        

        if len(params_list)!=correct_length:
            print("params_list of incorrect length should be of length ",correct_length)
            return 0
        #sensor
        A_s=params_list[0]
        B_s=params_list[1]
        C_s=params_list[2]
        N_s=params_list[3]
        #regulator
        A_r=params_list[4]
        B_r=params_list[5]
        C_r=params_list[6]
        N_r=params_list[7]
        #out_half
        A_h=params_list[8]
        B_h=params_list[9]
        C_h=params_list[10]
        #output
        A_o=params_list[11]
        B_o=params_list[12]
        C_o=params_list[13]
        N_o=params_list[14]
        #free
        F_o=params_list[15]
        
        Sensor = A_s+B_s*np.power(C_s*I_conc,N_s)
        Sensor /= 1+np.power(C_s*I_conc,N_s)

        Regulator = B_r/(1+np.power(C_r*Sensor,N_r))
        Regulator += A_r

        Output_half = B_h/(1+np.power(C_h*Sensor,N_o))
        Output_half += A_h

        Output = A_o*A_h + B_o*B_h/(1+np.power(C_h*(Sensor+C_o*Regulator),N_o))
        Output*=F_o
        #I wonder why we describe different repression strengths for repression by LacI_regulator and LacI_sensor?
        return Sensor,Regulator,Output_half, Output

    

#%%
    

class model_thermodynamic:
    def __init__(self,params_list:list,I_conc):
        self.params_list=params_list
        self.I_conc=I_conc
        self.example_dict={"sen_params":{"a_s":1,"Kp_s":1,"P":1,"Ki_s":1,"C_pi_s":1},"reg_params":{"a_r":1,"Kp_r":1,"Cr":1,"Nr":1},"out_h_params":{"Ah":1,"Bh":1,"Ch":1}}
    
    # thermodynamics model
    @staticmethod
    def model(params_list,I_conc):
        correct_length=13
        if len(params_list)!=correct_length:
            print("params_list of incorrect length should be of length ",correct_length)
            return 0
        # a_s, a_r, a_o represent the production rates divided by the degradation rates
        # I represents arabinose
        # P represents concentration of polymerase
        # Ki, Kii, Kps, Kpr, Ks, Kpo, K_lacI represent the binding affinity of the interaction
        # C_pi, C_ps, C_po_lacI represent the level of cooperative binding between activator or repressor with polymerase
        a_s=params_list[0]
        a_r=params_list[1]
        a_o=params_list[2]
        Ki_s=params_list[3]
        Kp_s=params_list[4]
        Kp_r=params_list[5]
        Ks_r=params_list[6]
        Kp_o=params_list[7]
        K_lacI_o=params_list[8]
        P=params_list[9]
        C_pi_s=params_list[10]
        C_ps_r=params_list[11]
        C_po_lacI_o=params_list[12]
        I=I_conc

        Sensor = a_s*(Kp_s*P+2*Ki_s*I*Kp_s*P+C_pi_s*Kp_s*Ki_s*P*np.power(I,2))
        Sensor /= (1+2*Ki_s*I+Kp_s*P+2*Ki_s*I*Kp_s*P+np.power(Ki_s*I,2)+C_pi_s*Kp_s*Ki_s*P*np.power(I,2))

        Regulator = a_r*(Kp_r*P+Kp_r*P*Ks_r*Sensor*C_ps_r)
        Regulator /= (1+Kp_r*P+Kp_r*P*Ks_r*Sensor*C_ps_r+Ks_r*Sensor)

        Output_half = a_o*(Kp_o*P+K_lacI_o*Sensor*Kp_o*P*C_po_lacI_o)
        Output_half /= (1+Kp_o*P+K_lacI_o*Sensor*Kp_o*P*C_po_lacI_o+K_lacI_o*Sensor)

        Output = a_o*(Kp_o*P+K_lacI_o*(Sensor+Regulator)*Kp_o*P*C_po_lacI_o)
        Output /= (1+Kp_o*P+K_lacI_o*(Sensor+Regulator)*Kp_o*P*C_po_lacI_o+K_lacI_o*(Sensor+Regulator))

        return Sensor, Regulator, Output_half, Output
#%%

class model_hill_shaky:
    def __init__(self,params_list:list,I_conc):
        self.params_list=params_list
        self.I_conc=I_conc
        self.example_dict={"sen_params":{"A_s":1,"B_s":1,"C_s":1,"N_s":1},"reg_params":{"A_r":1,"B_r":1,"C_r":1,"N_r":1},"out_h_params":{"A_h":1,"B_h":1,"C_h":1},"out_params":{"A_o":1,"B_o":1,"C_o":1,"N_o":1},"free_params":{"F_o":1}}
        self.correct_length=16
    @staticmethod
    def model(params_list:list,I_conc):
        #S is subscript for parameters corresponding to Sensor
        #R is subscript for parameters corresponding to Regulator
        #H is subscript for parameters corresponding to the half network I->S -| O
        #O is subscript for parameters corresponding to Output
        
        #creates variables described in params_dict 

        #sensor
        A_s=params_list[0]
        B_s=params_list[1]
        C_s=params_list[2]
        N_s=params_list[3]
        #regulator
        A_r=params_list[4]
        B_r=params_list[5]
        C_r=params_list[6]
        N_r=params_list[7]
        #out_half
        A_h=params_list[8]
        B_h=params_list[9]
        C_h=params_list[10]
        #output
        A_o=params_list[11]
        B_o=params_list[12]
        C_o=params_list[13]
        N_o=params_list[14]
        #free
        F_o=params_list[15]
        
        Sensor = np.array([])
        Regulator = np.array([])
        Output_half = np.array([])
        Output = np.array([])
        #initial conditions assumed as steady state with no inducer present
        S0 = A_r
        R0 = B_r/(1+np.power(C_r*S0,N_r))+ A_r
        H0 = B_h/(1+np.power(C_h*S0,N_o))+A_h
        O0 = A_o*A_h + B_o*B_h/(1+np.power(C_h*(S0+C_o*R0),N_o))*F_o
        SRHO0 = np.array([S0,R0,H0,O0])
        #arbitrary time point to integrate ODE up to
        t = np.array([1])
        #define system of ODEs to be solved by odeint, for a each inducer concentration
        for conc in I_conc:
            def ODEs(SRHO):
                #for set in params_dict.values():
                #    locals().update(set) 
                S = SRHO[0]
                R = SRHO[1]
                H = SRHO[2]
                O = SRHO[3]
                #S for sensor concentration at time t, prod for production
                S_prod = A_s+B_s*np.power(C_s*conc,N_s)
                S_prod /= 1+np.power(C_s*conc,N_s)
                #change in S concentration w.r.t. time, deg for degredation rate
                dSdt = S_prod - S

                R_prod = B_r/(1+np.power(C_r*S,N_r))
                R_prod += A_r

                dRdt = R_prod - R

                O_half_prod = B_h/(1+np.power(C_h*S,N_o))
                O_half_prod += A_h

                dHdt = O_half_prod - H

                O_prod = A_o*A_h + B_o*B_h/(1+np.power(C_h*(S+C_o*R),N_o))
                O_prod*=F_o #should Fo scale degredation as well?

                dOdt = O_prod - O
                return np.array([dSdt, dRdt, dHdt, dOdt])

            SRHO = odeint(ODEs, SRHO0, t)
            Sensor = np.append(Sensor, SRHO[0,0])
            Regulator = np.append(Regulator, SRHO[0,1])
            Output_half = np.append(Output_half, SRHO[0,2])
            Output = np.append(Output, SRHO[0,3])
        return np.array(Sensor),np.array(Regulator) ,np.array(Output_half), np.array(Output)
#%%

#function to minimize while fitting, expects a dictionary of parameters corresponding to the model of interest, 
def min_fun(params_list:list,data,model_type):
        log_sen = np.log10(data.Sensor)
        log_reg = np.log10(data.Regulator)
        log_out = np.log10(data.Output)
        log_stripe = np.log10(data.Stripe)   
        ind = data.S
        Sensor_est, Regulator_est, Output_est, Stripe_est = model_type(params_list,I_conc=ind)
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

#model_hill_shakey: 
#uses hill functions as in model_hill, but does not assume that a steady state has been reached in the system
#example input to model_hill_shakey:

#%% redefine model_hill_shaky
#I wonder why we describe different repression strengths for repression by LacI_regulator and LacI_sensor?
#%%
# Competition Model

# F(S,R(S),O(S,R(S)))=0
# F should be a function that takes one variable (S), and it's 0 when that's true 

# f(I), in our case I, represents production of S
# g(S) represents protduction of R
# K, d - constants 

def competition_model(params_list:list,I_conc):
    
    g=params_list[0]
    K=params_list[1]
    d=params_list[2]
    S=params_list[3]
    I=I_conc
    
    F = d*S*I-d*np.power(S,2)*K*g+d*np.power(S,2)*K*I-np.power(I,2)+I*K*g*S-2*K*np.power(I,2)*S+np.power(K,2)*I*g*np.power(S,2)-np.power(K,2)*np.power(I,2)*np.power(S,2)
    F /= K*I*(d*np.power(K,2)*I-d*S-I+K*g*S-K*I*S)
    
    return F