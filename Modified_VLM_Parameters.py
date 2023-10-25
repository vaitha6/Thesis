import numpy as np
from math import *


g_0=9.81 #Gravitational acceleration [m/s^2]
R=8.314 #Universal gas constant [J/mol/K]

#--------------------Propellant properties----------------------#

m_a_H2O=1.8e-2 #Molecular weight of water [kg]
gamma_H2O=1.32712 #Specific heat ratio of liquid water (Calculated)
R_H2O=461.5 #Characteristic gas constant (M.D Silva)
#L=2.25E+6 #Latent heat of vaporisation of water [J] (from BYJUS)
L=2123.50e+3 #Specific enthalpy of saturated water at a pressure of 5 bar [J/kg] (NIST steam tables)
Cp_l=4187 #Specific heat of liquid phase [J/kg/K] (GOOGLE)
T_s=424.981 #Saturation temeprature of water [K] (Saturation temp. at 5 bar: NIST)
Cp_v=1900 #Specific heat of vapor phase [J/kg/K]
#Cp_m=0 #Specific heat of mixture [J/kg/K]
k_l=0.61 #Thermal conductivity of liquid phase [J/K/s] (from NIST)
k_v=0.033 #Thermal conductivity of vapor phase [J/K/s] (from NIST)
#k_m=np.zeros(count) #Thermal conductivity of mixture [J/K/s]
rho_l=919 #917.02 #Density of saturated liquid phase at 5 bar[kg/m3] (from NIST steam tables)
rho_v=2.5 #2.6436 #Density of saturated vapor phase at 5 bar[kg/m3] (from NIST steam tables)
X_rho=rho_l/rho_v #Liquid to vapor denisty ratio
mu_l=3.5e-4 #Dynamic viscosity of liquid phase [Pa-s]
mu_v=1.6e-5 #Dynamic viscosity of vapor phase[Pa-s]


'''
#-----------Properties for verifivation (D. Fontanarosa)------------#

m_a_H2O=1.8e-2 #Molecular weight of water [kg]
gamma_H2O=1.32712 #Specific heat ratio of liquid water (Calculated)
R_H2O=461.5 #Characteristic gas constant (M.D Silva)
#L=2.25E+6 #Latent heat of vaporisation of water [J] (from BYJUS)
L=2210.57e+03 #2123.50e+3 #Specific enthalpy of saturated water at a pressure of 5 bar [J/kg] (NIST steam tables)
Cp_l=4187 #Specific heat of liquid phase [J/kg/K] (GOOGLE)
T_s=396.4 #424 #Saturation temeprature of water [K] (Saturation temp. at 5 bar: NIST)
Cp_v=1900 #Specific heat of vapor phase [J/kg/K]
#Cp_m=0 #Specific heat of mixture [J/kg/K]
k_l=0.6683 #Thermal conductivity of liquid phase [J/K/s] (from NIST)
k_v=0.033 #Thermal conductivity of vapor phase [J/K/s] (from NIST)
#k_m=np.zeros(count) #Thermal conductivity of mixture [J/K/s]
rho_l=943.12 #917.02 #Density of saturated liquid phase at 5 bar[kg/m3] (from NIST steam tables)
rho_v=1.2283 #2.6436 #Density of saturated vapor phase at 5 bar[kg/m3] (from NIST steam tables)
X_rho=rho_l/rho_v #Liquid to vapor denisty ratio
mu_l=3.5e-4 #Dynamic viscosity of liquid phase [Pa-s]
mu_v=1.6e-5 #Dynamic viscosity of vapor phase[Pa-s]
'''

fluid_prop_arr=[k_l, k_v, mu_l, mu_v, gamma_H2O, R_H2O]


#Viscosity model constants (Validity range:273K-643K)
A_mu=1.856E-14 #[Pa-s]
B_mu=4209 #[K]
C_mu=0.04527 # [/K]
D_mu=-3.376E-5 #[/K2]

mu_s=1e-4 #Viscosity at saturation conditions (ARBITRARY)

P_atm=0 #5000 #Atmospheric pressure at exit of nozzle [Pa] (Vacuum)
P_s=3E+5 #Mean pressure of chamber [Pa] (M.D Silva)
P_crit=2.2064e+7 #Critical pressure [Pa] (NIST)

#Model constants for Antoine equation
A=8.14019
B=1810.94
C=244.485

#Nozzle model constants

#alpha_1=1/(C + B/(A - log(P_s)/log(10)))**(1/2) - B/(2*log(10)*(C + B/(A - log(P_s)/log(10)))**(3/2)*(A - log(P_s)/log(10))**2)
#beta_1=(P_s/sqrt(C+(B/(A-log10(P_s)))))-alpha_1*P_s

alpha_1=0.047589938 #Obtained from curve fitting
beta_1=542.7841431 

#Model constants for Shomate equation
A_v=30.092
B_v=6.832514
C_v=6.793435
D_v=-2.53448
E_v=0.082139

A_l=-203.6060
B_l=1523.290
C_l=-3196.413
D_l=2474.455
E_l=3.855326

#Model constants for surface tension (International Tables of Surface Tension)
B_sigma=235.8E-3
b_sigma=-0.625
mu_sigma=1.256

#Constants for Kandlikar correlation

#Convective boiling
C1_kkrcb=1.136
C2_kkrcb=-0.9
C3_kkrcb=667.2
C4_kkrcb=0.7
C5_kkrcb=0

#Nucleate boiling
C1_kkncb=0.6683
C2_kkncb=-0.2
C3_kkncb=10589
C4_kkncb=0.7
C5_kkncb=0

#Nusselt number expression constants
A_Nu=0.1165
B_Nu=0.81
C_Nu=0.79
D_Nu=0.62
E_Nu=0.63

Nu_arr=[A_Nu, B_Nu, C_Nu, D_Nu, E_Nu]


#--------------------Silicon wafer properties----------------------#

m_a_Si=28.0855 #Molecular weight of silicon [kg]
rho_Si=2329 #Density of silicon [kg/m3]
k_Si=130 #Thermal conductivity of silicon [W/m/K]
Cp_Si=820 #Specific heat of silicon [J/kg/K]

#Model constants for Shomate equation for solid-phase silicon
A_Si=22.81719
B_Si=3.89951
C_Si=-0.082885
D_Si=0.042111
E_Si=-0.354063

t_Si=3.6e-04 #Thickness of silicon wafer [m] !!!!!!!!!!!!!!!!!!!!!

#--------------------Heating chamber geometry----------------------#

delta_s=1e-6 #Position incremental step within heating chamber [m]
L_lin=9e-03 #Total length of heating chamber [m] !!!CHANGE
s=delta_s
count=int((L_lin-s)/delta_s)

'''
#----------Geometry for coupled model verification (D. Fontanarosa)---------#


n_chamber=9
L_tot=6e-03
D_channel=8e-5
w_channel=1.2e-4
Asp_ch=D_channel/w_channel #Aspect ratio of channel
A_cs=n_chamber*D_channel*w_channel #Cross-sectional area of heating chamber [m2]
P_cs=n_chamber*2*(D_channel+w_channel) #Perimeter of heating chamber cross-section [m]
D_h=4*A_cs/P_cs #Hydraulic diameter [m] (General case)
A_ht=(P_cs*L_tot) + 9.6e-07
'''

'''
A_cs_inlet=2.4e-07 #Cross-sectional area of inlet chamber [m2]
P_cs_inlet=0.00424 #Perimeter of inlet chamber cross-section [m2]
D_h_inlet=4*A_cs_inlet/P_cs_inlet #Hydraulic diameter of inlet chamber [m]
Asp_ch_inlet=16.667 #Aspect ratio of inlet chamber [m2]
'''
'''
#-----Small serpentine channels-----#

d1=7.66e-5 #outer radius of channel for small serpentine geometry (SHOULD BE CHANGED LATER!) (M.D Silva)
d2=8.2e-6 #inner radius of channel for small serpentine geometry (M.D Silva)

D_channel=d1-d2 #Height of channels within chamber [m] 
w_channel=1.4815e-04 #Width of channel (ARBITRARY)
Asp_ch=D_channel/w_channel #Aspect ratio of channel
P0_Asp=48*(1-(1.13553*Asp_ch)+(1.9467*(Asp_ch**2))-(1.7012*(Asp_ch**3))+(0.9564*(Asp_ch**4))-(0.2537*(Asp_ch**5))) #Constant for pressure drop in channel
n_bends=floor(L_lin/(2*(d1+d2+D_channel))) #No. of bends in channel
n_chamber=21 #No. of channels in chamber (ARBITRARY)
A_cs=n_chamber*D_channel*w_channel #Cross-sectional area of chamber [m2]
P_cs=n_chamber*2*(D_channel+w_channel) #Perimeter of chamber cross-section [m]
D_h=4*A_cs/P_cs #Hydraulic diameter [m] (General case)
L_tot=pi*(d1+d2)*n_bends*0.5 #Actual length of a single channel [m]
A_ht=P_cs*L_tot #n_chamber*n_bends*pi*(d1+d2)*P_cs*delta_s/L_tot #Area of heat transfer through the channels [m2] !!!!!!!!!!!!!!!!!!!!!!!
#A_ht=5.16e
'''
'''
#-----Large serpentine channels (Fabricated)-----#

d1=2.897E-04 #2.66e-4 #outer radius of channel for small serpentine geometry (SHOULD BE CHANGED LATER!) (M.D Silva)
d2=3.99E-05 #5.4e-5 #inner radius of channel for small serpentine geometry (M.D Silva)
w_channel=4.44e-5 #2.1e-4 #Width of channel (ARBITRARY)
n_bends=14 #floor(L_lin/(2*(d1+d2+D_channel))) #No. of bends in channel
n_chamber=5 #No. of channels in chamber (ARBITRARY)
'''

#-----Large serpentine channels (Actual design)-----#

d1=2.66e-4 #outer radius of channel for small serpentine geometry (SHOULD BE CHANGED LATER!) (M.D Silva)
d2=5.4e-5 #inner radius of channel for small serpentine geometry (M.D Silva)
w_channel=4.44e-5 #2.1e-4 #Width of channel (ARBITRARY)
n_bends=14 #floor(L_lin/(2*(d1+d2+D_channel))) #No. of bends in channel
n_chamber=5 #No. of channels in chamber (ARBITRARY)

'''
#-----Small diamond channels-----#

d1=1.03e-5 #1.6e-04 #length of diamond geometry
d2=2.8e-06 #4e-5 #width of diamond geometry
w_channel=4.44e-5 #Width of channel (ARBITRARY)
n_bends=floor(L_lin/(1.2*d2)) #No. of bends in channel
n_chamber=floor(3e-3/(2*d1)) #No. of channels in chamber (ARBITRARY)
#A_ht=2.27e-5
'''
'''
#-----Large diamond channels-----#

d1=5.477e-4 #5.8e-04 #length of diamond geometry
d2=1.443e-04 #1.6e-4 #width of diamond geometry
w_channel=4.44e-5 #Width of channel (ARBITRARY)
n_bends=floor(L_lin/(1.1*d1)) #No. of bends in channel
n_chamber=floor(3E-03/(2*d2))+2 #No. of channels in chamber (ARBITRARY)
'''

LS_geom_arr=[d1, d2, w_channel, n_bends, n_chamber, L_lin, delta_s]


def LS_geom(LS_geom_arr):
    
    '''
    #FOR DIAMOND
    D_channel=d1/2 #Height of channels within chamber [m] 
    Asp_ch=D_channel/w_channel #Aspect ratio of channel
    P0_Asp=48*(1-(1.13553*Asp_ch)+(1.9467*(Asp_ch**2))-(1.7012*(Asp_ch**3))+(0.9564*(Asp_ch**4))-(0.2537*(Asp_ch**5))) #Constant for pressure drop in channel
    L_tot=L_lin 
    A_cs=3E-03*w_channel #n_chamber*D_channel*w_channel #Cross-sectional area of chamber [m2]
    P_cs=2*(3E-03 + w_channel) #n_chamber*2*(D_channel+w_channel) #Perimeter of chamber cross-section [m]
    D_h=4*A_cs/P_cs #Hydraulic diameter [m] (General case)
    A_ht=n_chamber*n_bends*(2*t_noz*sqrt((d1**2)+(d2**2)))+((L_lin*3E-03)-(n_chamber*n_bends*d1*d2)) #Area of heat transfer through the channels [m2] !!!!!!!!!!!!!!!!!!!!!!!
    '''
    #FOR SERPENTINE
    D_channel=LS_geom_arr[0]-LS_geom_arr[1] #Height of channels within chamber [m] 
    Asp_ch=D_channel/LS_geom_arr[2] #Aspect ratio of channel
    P0_Asp=48*(1-(1.13553*Asp_ch)+(1.9467*(Asp_ch**2))-(1.7012*(Asp_ch**3))+(0.9564*(Asp_ch**4))-(0.2537*(Asp_ch**5))) #Constant for pressure drop in channel
    A_cs=LS_geom_arr[4]*D_channel*LS_geom_arr[2] #Cross-sectional area of chamber [m2]
    P_cs=LS_geom_arr[4]*2*(D_channel+LS_geom_arr[2]) #Perimeter of chamber cross-section [m]
    D_h=4*A_cs/P_cs #Hydraulic diameter [m] (General case)
    #D_h=2*sqrt(LS_geom_arr[2]*D_channel/pi)
    L_tot=pi*(LS_geom_arr[0]+LS_geom_arr[1])*LS_geom_arr[3]*0.5 #Actual length of a single channel [m]
    A_ht=P_cs*L_tot #LS_geom_arr[4]*LS_geom_arr[3]*pi*(LS_geom_arr[0]+LS_geom_arr[1])*P_cs*delta_s/L_tot #Area of heat transfer through the channels [m2] !!!!!!!!!!!!!!!!!!!!!!!
    #A_ht=5.16e-6
    #count=int((L_lin-s)/delta_s)

    
    return D_channel, Asp_ch, A_cs, P_cs, D_h, L_tot, A_ht, L_lin, delta_s #,count

#---------------------------Nozzle geometry------------------------#

'''
#-----Micronozzle (C. Ganani)-----#

A_t=2e-9 #Throat area [m2] (M.D Silva)
AR=16.971 #Area ratio 
W_nc=3E-3 #Width of converging section of nozzle [m] (M.D Silva)
W_nd=7.2E-4 #Width of diverging section of nozzle [m] (M.D Silva)
l_nd=2.6195e-04 #Length of divergent section of nozzle [m]
W_t=4.5E-5 #Width of nozzle throat [m] (M.D Silva)
'''
'''
#--------Micronozzle (D. Fontanarosa)----------#

A_t=1.8e-08
A_e=2.112e-7
AR=A_e/A_t
t_noz=1.2e-04
W_nc=1.12e-04
W_nd=1.76e-03
W_t=1.5e-4
l_nd=0.5*(W_nd-W_t)/tan(15*pi/180)
'''
'''
#-----Long nozzle/ Bell nozzle (Fabricated)-----#

#A_t=2e-9 #Throat area [m2] (M.D Silva)
AR=11 #Area ratio 
W_nc=2.979E-03 #3E-3 #Width of converging section of nozzle [m] (M.D Silva)
W_nd=4.89E-04 #5E-4 #Width of diverging section of nozzle [m] (M.D Silva)
l_nd=6.26E-04 #6.45e-04 #Length of divergent section of nozzle [m]
W_t=2.51E-5 #4.5E-5 #Width of nozzle throat [m] (M.D Silva)
'''

#-----Long nozzle/ Bell nozzle (Actual design)-----#

A_t=2e-9 #Throat area [m2] (M.D Silva)
AR=11 #Area ratio 
W_nc=3E-3 #Width of converging section of nozzle [m] (M.D Silva)
W_nd=5E-4 #Width of diverging section of nozzle [m] (M.D Silva)
l_nd=6.45e-04 #Length of divergent section of nozzle [m]
W_t=4.5E-5 #Width of nozzle throat [m] (M.D Silva)

'''
#-----Short nozzle-----#

A_t=2e-9 #Throat area [m2] (M.D Silva)
AR=17 #Area ratio 
A_e=A_t*AR
W_nc=2.9804E-03 #3E-3 #Width of converging section of nozzle [m] (M.D Silva)
W_nd=7.77E-04 #7.8E-4 #Width of diverging section of nozzle [m] (M.D Silva)
l_nd=6.436E-04 #6.6E-04 #Length of divergent section of nozzle [m]
W_t=2.6E-5 #4.5E-5 #Width of nozzle throat [m] (M.D Silva)
r_c=W_t/2 #Radius of curvature at throat [m] (ARBITRARY)
'''


t_noz=4.44e-05 #thickness of nozzle [m] (Calculated based on dimensions from M.D Silva)

'''
A_t=W_t*t_noz
A_e=W_nd*t_noz

AR=W_nd/W_t
'''

nozzle_geom_arr=[A_t, AR, W_nc, W_nd, l_nd, W_t, t_noz]

def nozzle_geom(nozzle_geom_arr):
    
    A_e=nozzle_geom_arr[0]*nozzle_geom_arr[1]
    r_c=nozzle_geom_arr[5]/2 #Radius of curvature at throat [m] (C. Ganani)
    L_slant=sqrt((((nozzle_geom_arr[3]/2)-(nozzle_geom_arr[5]/2))**2)+(nozzle_geom_arr[4]**2))

    return A_e, r_c, L_slant

#------Transient formulation for simulating heating chamber--------#

no_it=2 #No. of iterations for transient simulation of heating chamber

K_wall=28.5 #Constant determined for chip temperature model (M.D Silva)
tau_wall=119.5 #Time constant for chip temperature model (M.D Silva)


