from math import *
import numpy as np
import matplotlib.pyplot as plt
from Modified_VLM_Parameters import P_atm,g_0,mu_v,W_t,L,t_noz,W_nd,A_t, AR, nozzle_geom

'''
W_t=1.5e-04
D_h_t=151.38801524e-06
r_c=W_t/2
t_noz=120e-06
W_nd=1760e-06
L_slant=((0.5*W_nd)-W_t)/sin(15*pi/180)
'''


def nozzle_model(gamma, R, T_v, P_v, nozzle_geom_arr, fluid_prop_arr):
    
    D_h_t=2*sqrt(nozzle_geom_arr[5]*nozzle_geom_arr[6]/pi)
    D_h_nd=2*sqrt(nozzle_geom_arr[3]*nozzle_geom_arr[6]/pi)

    A_e, r_c, L_slant = nozzle_geom(nozzle_geom_arr)
    k_l, k_v, mu_l, mu_v, gamma, R = fluid_prop_arr
    #AR=1.5
    #W_t=2.2e-04
    
    A_t=A_e/AR
    rho_v=P_v/(R*T_v)
    
    #Constants for calculating mach number
    A_temp=gamma+1
    B_temp=gamma-1
    C_temp=A_temp/(2*B_temp)
    
    vkchv=sqrt(gamma*(((2/(1+gamma))**((1+gamma)/(gamma-1)))))
    
    k_AR=((vkchv/AR)**2)*((gamma-1)/(2*gamma))
            
    P_e=1000
    PR=0.8
    temp=0
    
    while(abs(temp-P_e)>1e-6):
        PR=(k_AR/(1-(PR**((gamma-1)/gamma))))**(gamma/2)
        temp=P_e
        P_e=PR*P_v
        #temp=P_v*((1-(P_e/P_v))**((1-gamma)/gamma))*(k_AR**(gamma/2))
        #P_e=temp
    
    P_t=P_v*(2/(gamma+1))**(gamma/(gamma-1))
    
    T_e=T_v*(P_e/P_v)**((gamma-1)/gamma)
    T_t=T_v*(P_t/P_v)**((gamma-1)/gamma)
    
    M_e=sqrt((2*T_v/(T_e*(gamma-1)))*(1-((P_e/P_v)**((gamma-1)/gamma))))

    
    u_t=sqrt(gamma*R*T_t)
    rho_t=rho_v*((2/(gamma+1))**(1/(gamma-1)))
    Re_t=rho_t*u_t*D_h_t/mu_v
    
    m_dot_IRT=rho_t*A_t*u_t #vkchv*P_v*A_t/sqrt(R*T_v)
    
    u_e=sqrt((2*gamma/(gamma-1))*R*T_v*(1-(P_e/P_v)**((gamma-1)/gamma))) 
    rho_e=m_dot_IRT/(u_e*W_nd*t_noz)
    Re_e=rho_e*u_e*D_h_nd/mu_v

    v_ex_IRT=M_e*sqrt(gamma*R*T_e)#*Isp_eff #sqrt((2*gamma/(gamma-1))*R*T_v*((1-(P_e/P_v))**((gamma-1)/gamma)))

    F_mom=(m_dot_IRT*v_ex_IRT)
    F_pres=((P_e-P_atm)*A_e)   

    F_IRT=(F_mom+F_pres)
    
    C_d_const_arr=[0.05, 0.5, 0.75, 0.05, 0.019, 0.1, 0.5, 0.5, 0.21, 0.97, 0.86]
    
    C_d_1, C_d_2, C_d_3, C_d_4, C_d_5, C_d_6, C_d_7, C_d_8, C_d_9, C_d_10, C_d_11 = C_d_const_arr
    
    C_d=(((r_c+(C_d_1*C_d_2*W_t))/(r_c+(C_d_3*C_d_4*W_t)))**C_d_5)*(1-(((r_c+(C_d_6*C_d_7*W_t))/(C_d_8*W_t))**C_d_9)*(1/sqrt(Re_t))*(C_d_10+(C_d_11*gamma)))
    
    C_f_vis=17.6*exp(0.0032*A_e/A_t)/sqrt(0.773*Re_t)  
    C_f_IRT=F_IRT/(P_v*A_t)
    
    flow_qual=(C_f_IRT - C_f_vis)/C_f_IRT
    
    m_dot=m_dot_IRT*C_d
    
    Isp_eff_const_arr = [0.048, 0.037, 0.92, 1.49, 0.9, 7.09, 0.113, 0.29, 2]
    
    
    delta_mean=L_slant*Isp_eff_const_arr[0]/(Re_e**0.2) #Incompressible displacement thickness
    theta_mean=L_slant*Isp_eff_const_arr[1]/(Re_e**0.2) #Incompressible momentum thickness
    S_mean=delta_mean/theta_mean
    theta=theta_mean*(1-((Isp_eff_const_arr[2]*(tanh(Isp_eff_const_arr[3]*(S_mean-Isp_eff_const_arr[4])))*(M_e**2))/(Isp_eff_const_arr[5]+(M_e**2)))) #Compressible momentum thickness
    delta=((S_mean*(1+(Isp_eff_const_arr[6]*(M_e**2))))+(Isp_eff_const_arr[7]*(M_e**2)))*theta #Compressible displacement thickness
    Isp_eff=(1-((Isp_eff_const_arr[8]*delta/W_nd)*(1+(theta/delta)))) #Isp efficiency at nozzle exit
    
    v_ex=v_ex_IRT*Isp_eff
    
    F=(m_dot*v_ex)+F_pres #*flow_qual
    
    v_ex_eff=F/m_dot
 
    I_sp=v_ex_eff/g_0
    
    #print(' v_ex : ', v_ex, 'P_e : ', P_e, 'T_e: ', T_e)
    
    return m_dot,F,I_sp,Re_t,C_d,Isp_eff,theta,delta,M_e,Re_e,S_mean,M_e,P_e,P_t,T_t,T_e,v_ex #,F_pres,F_mom,PR, AR, Re_t, flow_qual



def backtrack_model(gamma, R, F, I_sp, m_dot, nozzle_geom_arr):
    
    A_e, r_c, L_slant = nozzle_geom(nozzle_geom_arr)
    
    #AR=49.2 
    
    vkchv=sqrt(gamma*(((2/(1+gamma))**((1+gamma)/(gamma-1)))))
    print('Vanderkenchov : ', vkchv)
    
    k_AR=((vkchv/AR)**2)*((gamma-1)/(2*gamma))

    PR=0.1
    temp=0
    
    while(abs(temp-PR)>1e-7):
        
        temp=PR
        PR=(k_AR/(1-(PR**((gamma-1)/gamma))))**(gamma/2)
     
    #P_v=4.8E+05
    #P_e=PR*P_v

    P_v=np.arange(1e+05, 1.01e+07, 1e+04)
    T_v=np.zeros(1000)
    P_e=np.zeros(1000)
    u_t=np.zeros(1000)
    trans_coeff=np.zeros(1000)
    v_ex=np.zeros(1000)
    
    print(np.shape(P_v)[0])
    
    for i in range(np.shape(P_v)[0]):
        
        P_e[i]=PR*P_v[i]
    
        v_ex_eff=I_sp*g_0
        v_ex[i]=v_ex_eff - (P_e[i]*A_e/m_dot) 
            
        T_v[i]=(gamma-1)*(v_ex[i]**2)/(2*gamma*R*(1-(PR**((gamma-1)/gamma))))
        
        P_t=P_v[i]*(2/(gamma+1))**(gamma/(gamma-1))
    
        T_t=T_v[i]*(P_t/P_v[i])**((gamma-1)/gamma)
        
        u_t[i]=sqrt(abs(gamma*R*T_t))
        
        trans_coeff[i]=F/(P_v[i]*A_t)
        
    plt.plot(P_v, T_v)
    plt.plot(P_v_2, T_v_2)
    plt.plot(4.8e+05, 436.03, marker="o", markersize=5, markeredgecolor="red", markerfacecolor="green")
    plt.plot(4.8e+05, 423.03, marker="o", markersize=5, markeredgecolor="blue", markerfacecolor="yellow")
    plt.legend(['Model', 'Literature', 'Model value for fixed chamber pressure', 'Literature value for fixed chamber pressure'])
    plt.xlabel('Chamber pressure [in bar]', fontsize=16)
    plt.ylabel('Chamber temperature [in K]', fontsize=16)
    '''
    v_ex_eff=I_sp*g_0
    v_ex=v_ex_eff - (P_e*A_e/m_dot) 
        
    T_v=(gamma-1)*(v_ex**2)/(2*gamma*R*(1-(PR**((gamma-1)/gamma))))
    
    P_t=P_v*(2/(gamma+1))**(gamma/(gamma-1))
    T_t=T_v*(P_t/P_v)**((gamma-1)/gamma)
    u_t=sqrt(gamma*R*T_t)
    '''
    
    return T_v, P_v, P_e, u_t # trans_coeff, v_ex

37.308922*tan(12*pi/180)

