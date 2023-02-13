import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

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

def model_4_pred(I_conc,As,Bs,Cs,Ns,Ar,Br,Cr,Nr,Ah,Bh,Ch,Fo,Ao,Bo,Co,No):
    #S is subscript for parameters corresponding to Sensor
    #R is subscript for parameters corresponding to Regulator
    #H is subscript for parameters corresponding to the half network I->S -| O
    #O is subscript for parameters corresponding to Output

    Sensor = As+Bs*np.power(Cs*I_conc,Ns)
    Sensor /= 1+np.power(Cs*I_conc,Ns)

    Regulator = Br/(1+np.power(Cr*Sensor,Nr))
    Regulator += Ar

    Output_half = Bh/(1+np.power(Ch*Sensor,No))
    Output_half += Ah

    Output = Ao*Ah + Bo*Bh/(1+np.power(Ch*(Sensor+Co*Regulator),No))
    Output*=Fo
    #I wonder why we describe different repression strengths for repression by LacI_regulator and LacI_sensor?
    return Sensor,Regulator,Output_half, Output



# thermodynamics model
def thermodynamic_model(a_s,a_r,a_h,a_o,Ki_s,Kii_s,Kp_s,Kp_r,Ks_r,Kp_o,K_lacI_o,I,P,C_ip_s,C_iip_s,C_ps_r,C_po_lacI_o):

    # a_s, a_r, a_oh, a_o represent the production rates divided by the degradation rates
    # I represents arabinose
    # P represents concentration of polymerase
    # Ki, Kii, Kps, Kpr, Ks, Kpo, K_lacI represent the binding affinity of the interaction
    # C_ip, C_iip, C_ps, C_po_lacI represent the level of cooperative binding between activator or repressor with polymerase

    Sensor = a_s*(Ki_s*I*Kp_s*P*C_ip_s+Kp_s*P+Kii_s*np.power(I,2)*Kp_s*P*C_iip_s)
    Sensor /= (1+Ki_s*I+Ki_s*I*Kp_s*P*C_ip_s+Kp_s*P+Kii_s*np.power(I,2)*Kp_s*P*C_iip_s)

    Regulator = a_r*(Kp_r*P+Kp_r*P*Ks_r*Sensor*C_ps_r)
    Regulator /= (1+Kp_r*P+Kp_r*P*Ks_r*Sensor*C_ps_r+Ks_r*Sensor)

    Output_half = a_h*(Kp_o*P+K_lacI_o*Sensor*Kp_o*P*C_po_lacI_o)
    Output_half /= (1+Kp_o*P+K_lacI_o*Sensor*Kp_o*P*C_po_lacI_o+K_lacI_o*Sensor)

    Output = a_o*(Kp_o*P+K_lacI_o*(Sensor+Regulator)*Kp_o*P*C_po_lacI_o)
    Output /= (1+Kp_o*P+K_lacI_o*(Sensor+Regulator)*Kp_o*P*C_po_lacI_o+K_lacI_o*(Sensor+Regulator))

    return Sensor, Regulator, Output_half, Output

    