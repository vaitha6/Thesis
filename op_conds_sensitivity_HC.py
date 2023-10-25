import numpy as np
from math import *
from Modified_VLM_Parameters import *
from Coupled_model import coupled_model
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.patheffects as patheffects


Pow_arr=np.linspace(1,9,9)
P_init_arr=np.linspace(1e+05, 6e+05, 20)
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

'''
m_dot_flat=np.zeros([P_init_len*20, Pow_len])
T_v_flat=np.zeros([P_init_len*20, Pow_len])
P_v_flat=np.zeros([P_init_len*20, Pow_len])
F_flat=np.zeros([P_init_len*20, Pow_len])
Isp_flat=np.zeros([P_init_len*20, Pow_len])
'''

def opcond_iterations():
    
    n=0

    for h in range(Pow_len-10):
        
        for i in range(P_init_len):
        
            for j in range(np.shape(T_0_arr)[0]):
                
                T_sat[i]=((B/(A-log10(P_init_arr[i]*760/101325)))-C)+273.15
                m_dot_arr[i,j], F_arr[i,j], Isp_arr[i,j], T_v_arr[i,j], P_v_arr[i,j], x_v_arr[i,j] = coupled_model(Pow_arr[h+10], T_0_arr[j]+100, P_init_arr[i], T_sat[i], T_0_arr[j], LS_geom_arr, nozzle_geom_arr, fluid_prop_arr)
                
                n+=1
        print(n)
        
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
    
    
def objective(weight_arr, F_to_Pow_val, Isp_val, mdot_val, P_init_val):
    
    #obj = 100*((weight_arr[0]*F_val) + (weight_arr[1]*Isp_val) + (weight_arr[2]*(mdot_val)) + (weight_arr[3]*Tv_val) + (weight_arr[4]*(1-P_init_val)) + (weight_arr[5]*(1-Pow_val)))/sum(weight_arr)
    obj = 100*((weight_arr[0]*F_to_Pow_val) + (weight_arr[1]*Isp_val) + (weight_arr[2]*(mdot_val)) + (weight_arr[3]*(1-P_init_val)))/sum(weight_arr)
    
    
    return obj
    
def obj_iterations():
    
    n=0
    increment=9
    switch=0
    n_pow=0
    
    F_arr_1=np.zeros([P_init_len, Pow_len])
    Pow_arr_1=np.zeros([P_init_len, Pow_len])
    F_lim=[]
    F_lim_x=[]
    F_lim_y=[]
    T_v_lim_x=[]
    T_v_lim_y=[]
    T_v_arr_1=np.zeros([P_init_len, Pow_len])
    F_to_Pow_arr=np.zeros([P_init_len, Pow_len])
    max_val=0
    
    for i in range(P_init_len):
        
        switch=0
        #switch_0=0
        
        for j in range(Pow_len):
            
            if(T_v_flat[increment+20*i,j]>582.435):
                T_v_arr_1[i,j]=nan
                Pow_arr_1[i,j]=nan
                
                if(switch==0):
                    T_v_lim_y.append(P_init_arr[i]*1e-05)
                    T_v_lim_x.append(Pow_arr[j]-1)
                    switch=1
                
            else:
                T_v_arr_1[i,j]=T_v_flat[increment+20*i,j]
                Pow_arr_1[i,j]=Pow_arr[j]
            '''
            if(Pow_arr[j]>4):
                Pow_arr_1=nan
            
            else:
                Pow_arr_1=Pow_arr[j]
            '''
            
            if(F_flat[increment+20*i,j]>1.2e-03):
                F_arr_1[i,j]=nan
                
                if(n_pow<=9):
                    n_pow+=1
                    F_lim_x.append(n_pow)
                    F_lim_y.append(P_init_arr[i]*1e-05 - 0.26315)      
            
            else:
                F_arr_1[i,j]=F_flat[increment+20*i,j]
            
            F_to_Pow_arr[i,j]=F_arr_1[i,j]/Pow_arr_1[i,j]
                
            #obj_arr[i,j]=objective([1, 5, 1, 1, 1, 1], F_arr_1[i,j]/np.amax(F_flat), Isp_flat[increment+20*i,j]/np.amax(Isp_flat), mdot_flat[increment+20*i,j]/np.amax(mdot_flat), T_v_arr_1[i,j]/np.amax(T_v_flat), P_init_arr[n]/np.amax(P_init_arr), Pow_arr_1/np.amax(Pow_arr))
            obj_arr[i,j]=objective([5, 5, 0, 0], F_to_Pow_arr[i,j]/np.amax(F_to_pow_flat), Isp_flat[increment+20*i,j]/np.amax(Isp_flat), mdot_flat[increment+20*i,j]/np.amax(mdot_flat), P_init_arr[n]/np.amax(P_init_arr))
            
            if(obj_arr[i,j]<100 and obj_arr[i,j]>max_val):
                
                max_val=obj_arr[i,j]
                max_P_init=P_init_arr[i]
                max_Pow=Pow_arr[j]
                
        #switch_0=1
        n=i%20

    fig_a, ax_a = plt.subplots()
    cs=ax_a.contourf(Pow_arr, P_init_arr*1e-05, obj_arr, levels=100)
    cbar = fig_a.colorbar(cs)
    cg1=ax_a.plot(T_v_lim_x, T_v_lim_y, '--', color='r', label='Max. chamber temperature')
    cg2=ax_a.plot(F_lim_x, F_lim_y, '--', color='k', label='Max. thrust')
    #ax_a.plot([pos_lim_P_init, pos_lim_Pow])
    ax_a.set_ylabel('Inlet pressure [bar]', fontsize=14)
    ax_a.set_xlabel('Input power [W]', fontsize=14)
    ax_a.set_xlim([1,9])
    ax_a.legend()
    fig_a.set_size_inches(8.5, 6.5)
    
    return max_val, max_P_init*1e-05, max_Pow
        
def HC_sens_plot():
    
    T_v_flat=np.loadtxt('Tv_flat_HC_opcond_4.dat') 
    P_v_flat=np.loadtxt('Pv_flat_HC_opcond_4.dat') 
    F_flat=np.loadtxt('F_flat_HC_opcond_4.dat')  
    Isp_flat=np.loadtxt('Isp_flat_HC_opcond_4.dat') 
    mdot_flat=np.loadtxt('mdot_flat_HC_opcond_4.dat') 
    
    F_to_pow_flat = np.zeros([400,9])
    
    for i in range(400):
        
        for j in range(9):
            
            F_to_pow_flat[i,j]=F_flat[i,j]/Pow_arr[j]
       
    fig_a, ax_a = plt.subplots()
    ax_a.contourf(P_init_arr*1e-05 ,T_0_arr ,T_v_arr_9 ,levels=100 , cmap='coolwarm')
    CS_a=ax_a.contour(P_init_arr*1e-05 ,T_0_arr ,T_v_arr_9 ,levels=[630, 640, 650, 670, 700, 730],colors='black')
    ax_a.set_xlabel('Inlet pressure [bar]', fontsize=14)
    ax_a.set_ylabel('Inlet temperature [K]' , fontsize=14)
    ax_a.clabel(CS_a, fontsize=14)
    
    fig_b, ax_b = plt.subplots()
    ax_b.contourf(P_init_arr*1e-05 ,T_0_arr ,P_v_arr_9*1e-05 ,levels=100 , cmap='rainbow')
    CS_b=ax_b.contour(P_init_arr*1e-05 ,T_0_arr ,P_v_arr_9*1e-05 ,levels=10, colors='black')
    ax_b.set_xlabel('Inlet pressure [bar]' , fontsize=14)
    ax_b.set_ylabel('Inlet temperature [K]', fontsize=14)
    ax_b.clabel(CS_b, fontsize=14)
    
    fig_c, ax_c = plt.subplots()
    ax_c.contourf(P_init_arr ,T_0_arr ,Isp_arr_9 ,levels=100 , cmap='rainbow')
    CS_c=ax_c.contour(P_init_arr ,T_0_arr ,Isp_arr_9 ,levels=[112, 115, 117.3, 119, 120, 122, 123.5, 125], colors='black')
    ax_c.set_xlabel('Inlet pressure [Pa]', fontsize=14)
    ax_c.set_ylabel('Inlet temperature [K]' , fontsize=14)
    ax_c.clabel(CS_c, fontsize=14)
    
    fig_d, ax_d = plt.subplots()
    ax_d.contourf(P_init_arr ,T_0_arr ,F_arr_9/1000 ,levels=100 , cmap='rainbow')
    CS_c=ax_d.contour(P_init_arr ,T_0_arr ,F_arr_9/1000 ,levels=10, colors='black')
    ax_d.set_xlabel('Inlet pressure [Pa]', fontsize=14)
    ax_d.set_ylabel('Inlet temperature [K]' , fontsize=14)
    ax_d.clabel(CS_c, fontsize=14)
    
    fig_e, ax_e = plt.subplots()
    ax_e.contourf(P_init_arr*1e-05 ,T_0_arr ,m_dot_arr_9 ,levels=100 , cmap='rainbow')
    CS_c=ax_e.contour(P_init_arr*1e-05 ,T_0_arr ,m_dot_arr_9 ,levels=10, colors='black')
    ax_e.set_xlabel('Inlet pressure [Pa]', fontsize=14)
    ax_e.set_ylabel('Inlet temperature [K]' , fontsize=14)
    ax_e.clabel(CS_c, fontsize=14)
    
    T_v_arr_8 = T_v_flat[:,7]
    T_v_arr_8=T_v_arr_8.reshape(20,20)
    T_v_arr_8=T_v_arr_8.T
    
    P_v_arr_9 = P_v_flat[:,8]
    P_v_arr_9=P_v_arr_9.reshape(20,20)
    P_v_arr_9=P_v_arr_9.T
    
    F_arr_9 = F_flat[:,8]
    F_arr_9=F_arr_9.reshape(20,20)
    F_arr_9=F_arr_9.T
    
    Isp_arr_9 = Isp_flat[:,8]
    Isp_arr_9=Isp_arr_9.reshape(20,20)
    Isp_arr_9=Isp_arr_9.T
    
    m_dot_arr_9 = mdot_flat[:,8]
    m_dot_arr_9=m_dot_arr_9.reshape(20,20)
    m_dot_arr_9=m_dot_arr_9.T