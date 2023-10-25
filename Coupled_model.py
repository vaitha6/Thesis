import numpy as np
from math import *
from Modified_VLM_Parameters import *
from Modified_VLM_model_ver3 import heating_chamber
from Modified_VLM_nozzle_model import nozzle_model
import matplotlib.pyplot as plt

'''
ASSUMPTIONS:
1) Flow through nozzle is isentropic. This implies that the adaibatic index is constant for all temperatures.
2) Critical conditions are present at the nozzle throat. (M_e at throat=1)
'''

'''
gamma=1.2176
R_H2O=375.3
AR=49.2
T_v=3267
P_v=5.93e+06
D_e=1.6
A_e=pi*(D_e**2)/4
'''

T_v_count=10

x_v=np.zeros(T_v_count)
T_v=np.zeros(T_v_count,dtype=float)
P_v=np.zeros(T_v_count,dtype=float)
m_dot_arr=np.zeros(T_v_count,dtype=float)
F=np.zeros(T_v_count)
I_sp=np.zeros(T_v_count)

T_exit=np.zeros([19,100])
P_exit=np.zeros([19,100])

F_out=np.zeros([19,100])
Isp_out=np.zeros([19,100])
m_out=np.zeros([19,100])
x=np.zeros([19,100])

Tw_init=np.arange(374,574,2)
power_inp=np.hstack([np.linspace(0.1,1,10),np.linspace(2,50,10)])

def coupled_model(power_input, Tw_initial, P_init, T_crit, T_0, LS_geom_arr, nozzle_geom_arr, fluid_prop_arr):
    
    #no_it=0
    
    for i in range(T_v_count-1):
                    
                    if(i==0):
                        m_dot_arr[i]=5e-06
                        
                    T, T_w, P, x_tp, x_vp, x, h_b, s_plot = heating_chamber(power_input, Tw_initial, P_init, T_crit, T_0, m_dot_arr[i], LS_geom_arr, fluid_prop_arr, Nu_arr)
                    T_v[i]=T[count-1]
                    P_v[i]=P[count-1]
                    x_v[i]=x[count-1]
                    
                    m_dot_arr[i+1], F[i+1], I_sp[i+1],Re_t,C_d,Isp_eff,theta,delta,M_e,Re_e,S_mean,M_e,P_e,P_t,T_t,T_e,v_ex = nozzle_model(gamma_H2O, R_H2O, T_v[i], P_v[i], nozzle_geom_arr, fluid_prop_arr)
                    
                    if(abs(m_dot_arr[i+1]-m_dot_arr[i])<1e-10):
                        
                        '''
                        print(' Iteration', no_it+1, '\n', 'mass flow rate : ', m_dot[i], ' Exit temp. : ', T_v[i], ' Exit pressure: ', P_v[i])
                        print('F : ',F[i+1], ' I_sp : ', I_sp[i+1], ' Exit velocity: ', v_ex)
                        '''
                        
                        #print('Solution converged!!!')
                        break
                    
                    '''
                    print(' Iteration', no_it+1, '\n', 'mass flow rate : ', m_dot[i], ' Exit temp. : ', T_v[i], ' Exit pressure: ', P_v[i])
                    print('F : ',F[i+1], ' I_sp : ', I_sp[i+1])
                    print('Solution NOT converged!!!')
                    '''
                    
                    #no_it+=1
                        
                        
    return m_dot_arr[i+1], F[i+1], I_sp[i+1], T_v[i], P_v[i], x_v[i] #, v_ex, P_e, x_tp, x_vp, s_plot
