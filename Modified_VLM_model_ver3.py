from math import *
import numpy as np
import matplotlib.pyplot as plt
#from steam_model import *
import scipy.integrate as integrate
from Density_model import density
from Modified_VLM_Parameters import *


#----------------Constants for wall temperature model------------------#

C1_Tw=50
C2_Tw=1
C3_Tw=118 
K_Tw=0.1

'''
rho[0]=rho_v
u[0]=m_dot/(rho[0]*A_cs) #Initial velocity at inlet [m/s]
Re[0]=rho[0]*u[0]*D_h/mu_v

row_no=np.zeros(count)
column_no=np.zeros(count)

Cp[0]=(A_v+(B_v*(T[0,0]/1000))+(C_v*((T[0,0]/1000)**2))+(D_v*((T[0,0]/1000)**3))+(E_v/((T[0,0]/1000)**2)))/m_a_H2O 
Pr[0]=mu_v*Cp[0]/k_v
Nu[0]=0.023*(Re[0]**0.8)*(Pr[0]**0.4)
h_b[0]=Nu[0]*k_v/D_h
Q_dot[0]=h_b[0]*P_cs*delta_s*(T_w[0,0]-T[0,0])
#Q_dot[0]=rho[0]*A_cs*u[0]*L
'''

'''
for j in range(no_it-1):
    
    time_it=0 
    s=0
    T[j,0]=305 #T_s
    '''

def heating_chamber(power_inp, Tw_init, P_init, T_sat, T_0, m_dot, LS_geom_arr, fluid_prop_arr, Nu_arr): #w_channel, D_chan

    time_it=0 
    s=0
    switch=0
    
    D_channel, Asp_ch, A_cs, P_cs, D_h, L_tot, A_ht, L_lin, delta_s = LS_geom(LS_geom_arr)
    k_l, k_v, mu_l, mu_v, gamma_H2O, R_H2O = fluid_prop_arr
    A_Nu, B_Nu, C_Nu, D_Nu, E_Nu = Nu_arr
    
    '''
    n_ch=floor(3e-03/(2*D_chan))
    A_cs=n_ch*w_channel*D_chan
    P_cs=2*n_ch*(w_channel + D_chan)
    A_ht=P_cs*L_lin
    D_h=4*A_cs/P_cs
    Asp_ch=D_chan/w_channel
    '''
    
    #-------------------------Array declaration------------------------#
    
    T=np.zeros(count)
    T_w=np.zeros(count)
    u=np.zeros(count)
    x=np.zeros(count)
    x_TD=np.zeros(count)
    x_incip=np.zeros(count)
    x_crit=np.zeros(count)
    rho=np.full(count,1000,dtype=float)
    sigma_t=np.zeros(count)
    h_b=np.zeros(count)
    h_bnb=np.zeros(count)
    h_bcb=np.zeros(count)
    P=np.zeros(count)
    pos=np.zeros(count)
    Cp=np.zeros(count)
    #Cp_Si=np.zeros(count)
    mu_H2O=np.zeros(count)
    Re=np.zeros(count)
    Pr=np.zeros(count)
    La=np.zeros(count)
    Nu=np.zeros(count)
    Bo=np.zeros(count)
    We=np.zeros(count)
    La=np.zeros(count)
    Co=np.zeros(count)
    k=np.zeros(count)
    PR=np.zeros(count)
    X_tt=np.zeros(count)
    h_enth=np.zeros(count)
    Q_dot=np.zeros(count)
    Q_frac=np.zeros(count)
    Q_dot_tot=np.zeros(no_it)
    time_per_it=np.zeros(no_it) #Cumulative time completed in the previous iterations [s] 
    s_plot=np.zeros(count)
    T_plot=[]
    Tw_plot=[]
    
    #-----------------------Initial Conditions--------------------------#
    
    P[0]=P_init #5E+5 #Pressure at inlet [Pa]
    PR[0]=P[0]/P_crit #Initial pressure ratio
    T_s=T_sat #424.981 #Saturation temperature at 5 bar [K] (NIST steam tables)
    T[0]=T_0 #Initial temperature at inlet [K] (ARBITRARY)
    T_w[0]=Tw_init #Initial wall temperature [K]
    P_in=power_inp #Input power [W] (SHOULD BE BASED ON WALL TEMPERATURE)
    T_c=647.15 #Critcal temperature of water [K] (NIST)
    q_in=P_in*P_cs*delta_s/A_ht
    rho[0]=rho_l #Initial density of propellant [kg/m3]
    sigma_t[0]=(B_sigma*((T_c-T[0])/T_c)**mu_sigma)*(1+(b_sigma*((T_c-T[0])/T_c))) #Surface tension between fluid and wall [N/m]
    #m_dot=1e+6 #1.78794e-8 #Mass flow rate (Coupled with nozzle model)
    u[0]=m_dot/(rho[0]*A_cs) #Initial velocity at inlet [m/s]
    k[0]=k_l #Initial thermal conductivity 
    Cp[0]=(A_l+(B_l*(T[0]/1000))+(C_l*((T[0]/1000)**2))+(D_l*((T[0]/1000)**3))+(E_l/((T[0]/1000)**2)))/m_a_H2O #Shomate equation [J/kg/K]
    mu_H2O[0]=mu_l #A_mu*exp((B_mu/T[0])+(C_mu*T[0])+(D_mu*(T[0]**2))) #Viscosity model [Pa-s] (Reid & Pausintz)
    #mu_H2O_s=A_mu*exp((B_mu/T_w)+(C_mu*T_w)+(D_mu*(T_w**2))) #Dynamic viscosity of water at wall temperature [Pa-s] (Used in calculation of Nusselt number)
    #Cp_Si[0]=(A_Si+(B_Si*(T[0]/1000))+(C_Si*((T[0]/1000)**2))+(D_Si*((T[0]/1000)**3))+(E_Si/((T[0]/1000)**2)))/m_a_Si
    Re[0]=rho[0]*u[0]*D_h/mu_H2O[0] #Reynolds number
    Re_v=rho_v*X_rho*u[0]*D_h/mu_v #Reynolds number of vapor phase
    Pr[0]=mu_H2O[0]*Cp[0]/k[0] #Prandtl number
    Pr_v=mu_v*Cp_v/k_v #Prandl number of vapor phase
    La[0]=sigma_t[0]*rho[0]*D_h/mu_H2O[0] #Laplace number 
    Nu[0]=A_Nu*((D_h/w_channel)**B_Nu)*(Asp_ch**C_Nu)*(Re[0]**D_Nu)*(Pr[0]**E_Nu) #1.86*((Re[0]*Pr[0]*D_h/L_tot)**(1/3))*((mu_H2O[0]/mu_s)**0.14) #Nusselt number (Seide-Tate for laminar flow, from HMT book)
    Bo[0]=Q_dot[0]/(A_ht*rho[0]*u[0]*L) #Boiling number
    We[0]=(rho[0]**2)*(u[0]**2)*D_h/(sigma_t[0]*rho[0]) #Webber number
    Eo=g_0*(rho_l-rho_v)*(D_h**2)/sigma_t[0] #Eotvos number
    La[0]=sigma_t[0]*rho[0]*D_h/mu_H2O[0] #Laplace number
    Fr_lo=u[0]**2/(g_0*D_h) #Froude number for liquid phase
    h_b[0]=Nu[0]*k[0]/D_h #Initial heat transfer coefficient [W/m2/K] 
    m_l_val=rho_l*A_cs*delta_s #Mass of liquid propellant within delta_s of the chamber[kg]
    m_v=rho_v*A_cs*delta_s #Mass of vapor propellant within delta_s of the chamber[kg]
    m_l=np.full(count,m_l_val)
    X_tt[0]=1 #Martinelli parameter
    x_drymax=0.95 #Maximum limit of vapor quality
    h_b_DB = 0.023*k_v*D_h*(Re_v**0.8)*(Pr_v**1/3) #Dittus-Boleter superheated vepor corelation
    #h_enth_lsat=Cp_l*(T_s-T[0]) #Saturation enthalpy of liquid propellant [J/kg]
    x[0]=0 #Initial vapor quality
    
    x_tp=0
    x_vp=0
        
    for i in range(int(count-1)):
        
        if(T_0>T_sat and i==0): switch=1

        if(T[i]<T_s):
             
             T[i+1]=T[i]+(h_b[i]*P_cs*delta_s*(T_w[i]-T[i])/(rho[i]*A_cs*Cp[i]*u[i])) #Increment in temperature of liquid propellant
             Cp[i+1]=(A_l+(B_l*(T[i+1]/1000))+(C_l*((T[i+1]/1000)**2))+(D_l*((T[i+1]/1000)**3))+(E_l/((T[i+1]/1000)**2)))/m_a_H2O #Specific heat capacity based on Shomate equation [NIST]
             #Cp_Si[i]=(A_Si+(B_Si*(T_w[i]/1000))+(C_Si*((T_w[i]/1000)**2))+(D_Si*((T_w[i]/1000)**3))+(E_Si/((T_w[i]/1000)**2)))*1000/m_a_Si
             u[i+1]=u[i]
             rho[i+1]=rho[i]
             P[i+1]=P[i]
             Q_dot[i]=h_b[i]*P_cs*delta_s*(T_w[i]-T[i])
             Re[i+1]=rho[i+1]*u[i+1]*D_h/mu_l
             Pr[i+1]=mu_l*Cp[i+1]/k_l
             Nu[i+1]=A_Nu*((D_h/w_channel)**B_Nu)*(Asp_ch**C_Nu)*(Re[i+1]**D_Nu)*(Pr[i+1]**E_Nu)
             h_b[i+1]=Nu[i+1]*k_l/D_h
             #h_b[i+1]=h_b[i]
             #Q_dot[i+1]=rho[i+1]*A_cs*u[i+1]*Cp[i+1]*(T[i+1]-T[i])
             #q_dot[i+1]=Q_dot[i+1]/(P_cs*delta_s)
             #Q_dot2[i]=Q_dot[i]+(0.5*m_dot*((u[i+1]**2)-(u[i]**2)))
             #T_w[i+1]=T[i+1]+(Q_dot[i]/(h_b[i+1]*P_cs*delta_s)) 
             T_w[i+1]=T_w[i]+(q_in-Q_dot[i])/(Cp_Si*rho_Si*t_Si*w_channel*u[i])
    
                
        elif(T[i]>T_s):
            
            if(x[i]<1 and T_0<T_sat):
                    
                    if(switch==0):
                        
                        Q_dot[i]=Q_dot[i-1]
                        h_b[i]=h_b[i-1] #Arbitrary heat transfer coefficient for two phase flow [W/m2K]
                        u_s_v=rho_l*u[i-1]/rho_v
                        Re_s_v=rho_v*u_s_v*D_h/mu_v
                        Cp_s_v=(A_v+(B_v*(T[i]/1000))+(C_v*((T[i]/1000)**2))+(D_v*((T[i]/1000)**3))+(E_v/((T[i]/1000)**2)))/m_a_H2O #Shomate equation using constants for vapor phase
                        Cp_s_l=(A_l+(B_l*(T[i]/1000))+(C_l*((T[i]/1000)**2))+(D_l*((T[i]/1000)**3))+(E_l/((T[i]/1000)**2)))/m_a_H2O
                        Pr_s_v=mu_v*Cp_s_v/k_v
                        Nu_s_v=0.23*(Re_s_v**0.8)*(Pr_s_v**0.4)
                        h_b_sv=Nu_s_v*k_v/D_h
                        h_b_sl=h_b[i] 
                        x_tp=s/L_lin
                        
                    switch=1
                    switch2=1
                    
                    '''
                    u[i+1]=u[i]
                    n_u=0
                    x_current=0
                    
                    for n_u in range(50):
                        
                        x_current=(1-x[i])*((h_b_sl*P_cs*delta_s*(T_w[i]-T[i])))/(rho[i]*A_cs*u[i]*L)
                        
                        T[i+1]=T[i]+((x[i]*h_b[i]*P_cs*delta_s*(T_w[i]-T[i]))-(0.5*x[i]*rho[i]*A_cs*u[i]*((u[i+1]**2)-(u[i]**2))))/Cp[i]
                        #T[i+1]=(T_temp+(2*T[i]))/3
                        
                        rho[i+1]=density(T[i+1]-273.15,P[i]/1e+06)
                    
                        temp_u=u[i+1]
                        u[i+1]=u[i]*rho[i]/rho[i+1]
                        
                        x[i+1]=x[i]+x_current
                        
                        if(abs(temp_u-u[i+1])<1e-03):
                            print('u = ',u[i+1], ' for x = ', x[i+1])               
                            break
                        
                        n_u+=1
                    
                    '''
                    
                    T[i+1]=T[i]
                    
                    Q_frac[i+1]=min(1,abs(Q_frac[i]+(P_cs*delta_s*h_b[i]*(T_w[i]-T[i])/(m_dot*L)))) #Amount of power used to partially vaporise propellant
                    #x[i+1]=(1/(rho_l-rho_v))*(rho_l-((h_b[i]*A_ht*(T_w[i]-T[i]))/(A_cs*u[i]*L)))
                    rho[i+1]=(min(Q_frac[i+1],1)*rho_v)+(max((1-Q_frac[i+1]),0)*rho_l)
                    x[i+1]=(rho_l-rho[i+1])/(rho_l-rho_v)
                    u[i+1]=rho[i]*u[i]/rho[i+1] #Velocity increment based on mass conservation

                    Q_dot[i+1]=(Q_frac[i+1]-Q_frac[i])*rho[i+1]*A_cs*u[i+1]*L
                    #Q_dot[i+1]=((x[i]*h_b[i])+((1-x[i])*h_b_sl))*P_cs*delta_s*(T_w[i]-T[i])
                    Cp[i+1]=(min(Q_frac[i+1],1)*Cp_s_v)+(max((1-Q_frac[i+1]),0)*Cp_s_l)
                    #Cp[i+1]=(A_v+(B_v*(T[i+1]/1000))+(C_v*((T[i+1]/1000)**2))+(D_v*((T[i+1]/1000)**3))+(E_v/((T[i+1]/1000)**2)))/m_a_H2O
                    #Cp_Si[i]=(A_Si+(B_Si*(T_w[i]/1000))+(C_Si*((T_w[i]/1000)**2))+(D_Si*((T_w[i]/1000)**3))+(E_Si/((T_w[i]/1000)**2)))/m_a_Si
                    P[i+1]=P[i]-(rho[i]*u[i]*(u[i+1]-u[i]))
                    mu_H2O[i+1]=(min(x[i+1],1)*mu_v)+(max((1-x[i+1]),0)*mu_l) #Weighted average of viscosity  
                    k[i+1]=(min(x[i+1],1)*k_v)+(max((1-x[i+1]),0)*k_l) #Weighted average of thermal conductivity
                    
                    
                    Re[i+1]=rho[i+1]*u[i+1]*D_h/mu_H2O[i+1]
                    Pr[i+1]=mu_H2O[i+1]*Cp[i+1]/k[i+1]
                    Nu[i+1]=A_Nu*((D_h/w_channel)**B_Nu)*(Asp_ch**C_Nu)*(Re[i+1]**D_Nu)*(Pr[i+1]**E_Nu)
                    
                    h_b[i+1]=Nu[i+1]*k[i+1]/D_h #(min(Q_frac[i+1],1)*h_b_sv)+(max((1-Q_frac[i+1]),0)*h_b_sl)    
                    
                    T_w[i+1]=T_w[i]+(q_in-Q_dot[i+1])/(Cp_Si*rho_Si*t_Si*w_channel*u[i])
                
            elif(x[i]>=1 or T_0>T_sat):
                
                '''
                if(T[i]>=T_c): 
                    print('Critical temperature reached!')
                    break
                '''
                
                if(switch==1):
                    x_vp=(s+delta_s)/L_lin
                    T[i+1]=T[i] #Increment in temperature of vapor propellant
                    
                    rho[i+1]=density(T[i+1]-273.15,P[i]/1e+06)
                    
                    u[i+1]=u[i]*rho[i]/rho[i+1]
                    
                    P[i+1]=P[i]-(rho[i]*u[i]*(u[i+1]-u[i]))
                    #P[i+2]=P[i+1]
                    #print(u[i+1])
                    
                elif(switch==2):
                    
                    '''
                    T[i+1]=T[i]+(h_b[i]*P_cs*delta_s*(T_w[i]-T[i])/(rho[i]*A_cs*Cp[i]*u[i]))
                    rho[i+1]=density(T[i+1]-273.15,P[i+1]/1e+06)
                    u[i+1]=u[i]*rho[i]/rho[i+1]
                    P[i+1]=P[i]-(rho[i]*u[i]*(u[i+1]-u[i]))
                    
                    '''
                    u[i+1]=0
                
                    for n_u in range(50):
                        
                        T[i+1]=T[i]+(h_b[i]*P_cs*delta_s*(T_w[i]-T[i])/(rho[i]*A_cs*Cp[i]*u[i]))-(((u[i+1]**2)-u[i]**2)/(2*Cp[i])) #Increment in temperature of vapor propellant
                        
                        rho[i+1]=density(T[i+1]-273.15,P[i+1]/1e+06)
                        
                        temp_u=u[i+1]
                        u[i+1]=u[i]*rho[i]/rho[i+1]
                        
                        P[i+1]=P[i]-(rho[i]*u[i]*(u[i+1]-u[i]))
                        
                        if(abs(temp_u-u[i+1])<1e-03):
                            break
                        
                        n_u+=1
                    
                    
                switch=2
                
                Q_frac[i+1]=Q_frac[i]
                x[i+1]=1                 
                
                Cp[i+1]=(A_v+(B_v*(T[i+1]/1000))+(C_v*((T[i+1]/1000)**2))+(D_v*((T[i+1]/1000)**3))+(E_v/((T[i+1]/1000)**2)))/m_a_H2O #Shomate equation using constants for vapor phase
                    
                #Cp_Si[i]=(A_Si+(B_Si*(T_w[i]/1000))+(C_Si*((T_w[i]/1000)**2))+(D_Si*((T_w[i]/1000)**3))+(E_Si/((T_w[i]/1000)**2)))/m_a_Si
                                        
                #P[i+1]=P[i]-((rho[i+1]-rho[i])*u[i]*(u[i+1]-u[i])) #T[i+1]*P[i]/T[i]
                
                Re[i+1]=rho[i+1]*u[i+1]*D_h/mu_v
                Pr[i+1]=mu_v*Cp[i+1]/k_v
                Nu[i+1]=A_Nu*((D_h/w_channel)**B_Nu)*(Asp_ch**C_Nu)*(Re[i+1]**D_Nu)*(Pr[i+1]**E_Nu)
                h_b[i+1]=Nu[i+1]*k_v/D_h
                
                Q_dot[i]=rho[i]*u[i]*A_cs*Cp[i]*(T[i+1]-T[i])
                #q_dot[i+1]=Q_dot[i+1]/(P_cs*delta_s)
                #Q_dot[i+1]=(0.5*m_dot*((u[i+1]**2)-(u[i]**2)))
                #T_w[i+1]=T[i+1]+(Q_dot[i]/(h_b[i+1]*P_cs*delta_s))
                T_w[i+1]=T_w[i]+(q_in-Q_dot[i+1])/(Cp_Si*rho_Si*t_Si*w_channel*u[i])
            
        s=s+delta_s
        s_plot[i]=s/L_lin
        time_it=time_it+(delta_s/u[i+1])
        
    #time_per_it[j]=time_it
    
    #T_w[j+1,i] = 473 C1_Tw*exp((-s/L_tot)*(1/C2_Tw)) + C3_Tw*exp(-time_per_it[j]/tau) + T[j,0]
    
    #print('Wall temp for iteration ',j+1,' : ',T_w[0])
    
    #T_w[j+1,0]=T_w[j,0]+((L_tot/(np.sum(u)/count))*((K_wall*P_net[j])-(T_w[j,0]/tau_wall)))
    
    s_plot[count-1]=1
    
    #print('Volume of vapor : ', (1-x_vp)*L_lin*w_channel*2.7e-3)
    
    return T, T_w, P, x_tp, x_vp, x, h_b, s_plot #, Re, Pr, D_h, A_cs, A_ht, Asp_ch

'''
h_b_arr=[1e+4,2e+4,5e+4,1e+5,2e+5,4e+5,6e+5,8e+5,1e+6]
heating_c
T_out=np.zeros(9,dtype=float)
P_out=np.zeros(9,dtype=float)

for k in range(9):
    
    T_out[k],P_out[k]=heating_chamber(m_dot,h_b_arr[k])
'''
  
def plotting(x_tp,x_vp):
    
    fig1, ax1 = plt.subplots()
    ax1.plot(s_plot,P*1e-05,'g', linewidth=2)
    ax1.set_xlabel('Fractional length along the channel (s/L_tot)', fontsize=16)
    ax1.set_ylabel('Chamber pressure (in bar)', fontsize=16)
    ax1.vlines(x=x_tp, ymin=P[count-1]/1e+05, ymax=P[0]/1e+05, colors='black', ls='--', lw=1, label='Onset of two-phase flow')
    ax1.vlines(x=x_vp, ymin=P[count-1]/1e+05, ymax=P[0]/1e+05, colors='black', ls='-.', lw=1, label='Onset of vapor phase flow')
    ax1.legend()
    
    fig2, ax2 = plt.subplots()
    ax2.plot(s_plot,T_w,'g', label='Wall temperature', linewidth=2)
    ax2.plot(s_plot,T,'b', label='Chamber temperature', linewidth=2)
    ax2.set_xlabel('Fractional length along the channel (s/L_tot)', fontsize=16)
    ax2.set_ylabel('Temperature (in K)', fontsize=16)
    ax2.vlines(x=x_tp, ymin=T[0], ymax=1.01*max(T[count-1],T_w[count-1]), colors='black', ls='--', lw=1, label='Onset of two-phase flow')
    if(x_vp!=0):
        ax2.vlines(x=x_vp, ymin=T[0], ymax=1.01*max(T[count-1],T_w[count-1]), colors='black', ls='-.', lw=1, label='Onset of vapor phase flow')
    ax2.legend()
    
def Cp_l(T_var,A_l,B_l,C_l,D_l,E_l):
    return (A_l+(B_l*(T_var/1000))+(C_l*((T_var/1000)**2))+(D_l*((T_var/1000)**3))+(E_l/((T_var/1000)**2)))/m_a_H2O #Specific heaT_var capaciT_vary based on ShomaT_vare equaT_varion [NIST_var]


    
