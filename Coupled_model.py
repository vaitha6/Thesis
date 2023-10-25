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

#u_v=20 #ARBITRARY VELOCITY AT EXIT [m/s]


def iterations(P_init,T_crit):
    
    #----------------------Iterations for nozzle -----------------------#

    plot_switch=0
    no_it=0
    
    
    for p in range(np.shape(power_inp)[0]):
        
        for q in range(np.shape(Tw_init)[0]):
            
            
            m_out[p,q], F_out[p,q], Isp_out[p,q], T_exit[p,q], P_exit[p,q], x[p,q] = coupled_model(power_inp[p],Tw_init[q],P_init,T_crit)
            
            print(' Iteration', no_it+1, '\n', 'mass flow rate : ', m_out[p,q], ' Exit temp. : ', T_exit[p,q], ' Exit pressure: ', P_exit[p,q])
            print('F : ',F_out[p,q], ' I_sp : ', Isp_out[p,q])
                        
                        #plotting(T,T_w,P,x_tp,x_vp,s_plot)          
                        #plotting(T,T_w,P,x_tp,x_vp,s_plot)
            no_it+=1
    
    np.savetxt('P_exit.dat',P_exit)
    np.savetxt('T_exit.dat',T_exit)
    np.savetxt('F_out.dat',F_out)
    np.savetxt('Isp_out.dat',Isp_out)
    np.savetxt('m_out.dat',m_out)

def final_plot():
    
    T_switch=0
    
    T_exit = np.genfromtxt('T_exit.dat',dtype=float,delimiter=' ')
    P_exit = np.genfromtxt('P_exit.dat',dtype=float,delimiter=' ')
    F_out = np.genfromtxt('F_out.dat',dtype=float,delimiter=' ')
    Isp_out = np.genfromtxt('Isp_out.dat',dtype=float,delimiter=' ')
    m_out = np.genfromtxt('m_out.dat',dtype=float,delimiter=' ')
    
    for i in range(np.shape(T_exit)[0]):
        for j in range(np.shape(T_exit)[1]):
            if(T_exit[i,j]>=647.15 and T_switch==0):
                T_index=i
                T_switch=1
                break

    #pow_lim=power_inp[T_index]


    fig1, ax1 = plt.subplots()
    CS=ax1.contour(Tw_init,power_inp,T_exit,levels=[350, 400, 425, 450, 470, 490, 520, 550, 575, 600, 647] ,colors='k')
    ax1.clabel(CS, inline=True, fontsize=10)
    ax1.set_xlabel('Initial wall temperature (in K)')
    ax1.set_ylabel('Input power (in W)')
    ax1.set_ylim([0, 10]) #3.6])
    #ax1.hlines(pow_lim,xmin=374,xmax=473,colors='b',linestyles='dashdot')
    

    fig2, ax2 = plt.subplots()
    CS=ax2.contour(Tw_init,power_inp,P_exit,levels=[499980.55, 499980.6, 499980.7, 499980.75, 499980.8, 499990, 499999.9, 499999.99, 499999.9999],colors='k')
    ax2.clabel(CS, inline=True, fontsize=10)
    ax2.set_xlabel('Initial wall temperature (in K)')
    ax2.set_ylabel('Input power (in W)')
    ax2.set_ylim([0, 10]) #3.6])
    
    fig3, ax3 = plt.subplots()
    CS=ax3.contour(Tw_init,power_inp,F_out*1000,levels=[1.425, 1.427, 1.429, 1.431, 1.433, 1.435, 1.437, 1.439, 1.441, 1.443, 1.445, 1.446, 1.447875, 1.45, 1.455, 1.46],colors='k')
    ax3.clabel(CS, inline=True, fontsize=10)
    ax3.set_xlabel('Initial wall temperature (in K)')
    ax3.set_ylabel('Input power (in W)')
    ax3.set_ylim([0, 10]) #3.6])
    
    fig4, ax4 = plt.subplots()
    CS=ax4.contour(Tw_init,power_inp,Isp_out,levels=[94, 96, 97.25, 101, 103, 105, 108, 111, 114, 117, 120, 123],colors='k')
    ax4.clabel(CS, inline=True, fontsize=10)
    ax4.set_xlabel('Initial wall temperature (in K)')
    ax4.set_ylabel('Input power (in W)')
    ax4.set_ylim([0, 10]) #3.6])

    
    fig5, ax5 = plt.subplots()
    CS=ax5.contour(Tw_init,power_inp,m_out,levels=[1.08e-06, 1.1e-06, 1.25e-06, 1.3e-06, 1.35e-06, 1.4e-06, 1.45e-06, 1.5e-06, 1.51765e-06, 1.6e-06, 1.67e-06],colors='k')
    ax5.clabel(CS, inline=True, fontsize=10)
    ax5.set_xlabel('Initial wall temperature (in K)')
    ax5.set_ylabel('Input power (in W)')
    ax5.set_ylim([0, 10]) #3.6])

    
    return T_exit, P_exit, F_out, Isp_out, m_out

