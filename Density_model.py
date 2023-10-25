import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

coeff_poly=[0.598545543973339,-0.313493989303638,0.306467715627611,0.201221026824129,-0.184094624273186,-0.000623902597538823,-0.0817303083105097,0.0921628581430544,0.000196756729517733,-0.00345918376965264,0.0118214541372334,-0.0163383494390289,-0.000218494256306388,5.04564421253956e-05,0.00135069389431253]
   
def density(T_data,P_data):

    rho=(coeff_poly[0] + coeff_poly[1]*((T_data-557.6)/407.4) + coeff_poly[2]*((P_data-0.2336)/0.1158) + coeff_poly[3]*(((T_data-557.6)/407.4)**2) + coeff_poly[4]*((T_data-557.6)/407.4)*((P_data-0.2336)/0.1158) + coeff_poly[5]*(((P_data-0.2336)/0.1158)**2) + coeff_poly[6]*(((T_data-557.6)/407.4)**3) + coeff_poly[7]*(((T_data-557.6)/407.4)**2)*((P_data-0.2336)/0.1158) + coeff_poly[8]*((T_data-557.6)/407.4)*(((P_data-0.2336)/0.1158)**2) + coeff_poly[9]*(((P_data-0.2336)/0.1158)**3) + coeff_poly[10]*(((T_data-557.6)/407.4)**4) + coeff_poly[11]*(((T_data-557.6)/407.4)**3)*((P_data-0.2336)/0.1158) +  coeff_poly[12]*(((T_data-557.6)/407.4)**2)*(((P_data-0.2336)/0.1158)**2) +  coeff_poly[13]*((T_data-557.6)/407.4)*(((P_data-0.2336)/0.1158)**3) + (coeff_poly[14]*(((P_data-0.2336)/0.1158)**4)))  
    return rho

T_dat=np.arange(200,400,1)
P_dat=np.arange(10,2010,10)
rho=np.zeros([np.shape(T_dat)[0],np.shape(P_dat)[0]])

for i in range(np.shape(T_dat)[0]):
    
    for j in range(np.shape(P_dat)[0]):
    
        rho[i,j] = density(T_dat[i], P_dat[j])

'''        
fig = plt.figure()
ax = plt.axes(projection='3d')    

ax.plot3D(T_dat, P_dat, rho) 
'''