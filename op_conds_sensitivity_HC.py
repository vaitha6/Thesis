import numpy as np
from math import *
from Modified_VLM_Parameters import *
from Coupled_model import coupled_model
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.patheffects as patheffects


Pow_arr=np.linspace(1,9,50)
P_init_arr=np.linspace(1.5e+05, 6e+05, 50)
Pow_len=np.shape(Pow_arr)[0]
P_init_len=np.shape(P_init_arr)[0]
T_0_arr=np.linspace(300,500,20)
obj_arr=np.zeros([P_init_len, Pow_len])
T_sat=np.zeros(P_init_len)

F_arr=np.zeros([P_init_len, 20])
Isp_arr=np.zeros([P_init_len, 20])
m_dot_arr=np.zeros([P_init_len, 20])
x_v_arr=np.zeros([P_init_len, 20])
T_v_arr=np.zeros([P_init_len, 20])
P_v_arr=np.zeros([P_init_len, 20])

m_dot_flat=np.zeros([P_init_len*20, Pow_len])
T_v_flat=np.zeros([P_init_len*20, Pow_len])
P_v_flat=np.zeros([P_init_len*20, Pow_len])
F_flat=np.zeros([P_init_len*20, Pow_len])
Isp_flat=np.zeros([P_init_len*20, Pow_len])

opcond_iterations()

def opcond_iterations():
    
    n=0

    for h in range(Pow_len):
        
        for i in range(P_init_len):
        
            for j in range(np.shape(T_0_arr)[0]):
                
                T_sat[i]=((B/(A-log10(P_init_arr[i]*760/101325)))-C)+273.15
                m_dot_arr[i,j], F_arr[i,j], Isp_arr[i,j], T_v_arr[i,j], P_v_arr[i,j], x_v_arr[i,j] = coupled_model(Pow_arr[h], T_0_arr[j]+100, P_init_arr[i], T_sat[i], T_0_arr[j], LS_geom_arr, nozzle_geom_arr, fluid_prop_arr)
                
                n+=1
        print('Iteration: ', n)
        
        m_dot_flat[:,h+10]=m_dot_arr.flatten()
        T_v_flat[:,h+10]=T_v_arr.flatten()
        P_v_flat[:,h+10]=P_v_arr.flatten()
        F_flat[:,h+10]=F_arr.flatten()
        Isp_flat[:,h+10]=Isp_arr.flatten()
    

    np.savetxt('F_flat_HC_opcond_5(mod).dat', F_flat)
    np.savetxt('Isp_flat_HC_opcond_5(mod).dat', Isp_flat)
    np.savetxt('m_dot_flat_HC_opcond_5(mod).dat', m_dot_flat)
    np.savetxt('T_v_flat_HC_opcond_5(mod).dat', T_v_flat)
    np.savetxt('P_v_flat_HC_opcond_5(mod).dat', P_v_flat)
