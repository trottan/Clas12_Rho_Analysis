import numpy as np
from array import array
import ROOT



splitSize = len(rhoAmp_p)/9
Np = np.array(np.split(np.array(rhoAmp_p),splitSize))
Nm = np.array(np.split(np.array(rhoAmp_m),splitSize))


Np_e =np.array(np.split( np.array(rhoAmpError_p),splitSize))
Nm_e =np.array(np.split(np.array(rhoAmpError_m),splitSize))

print(len(Nm),len(Np),len(Nm_e),len(Np_e))

Pb = 0.8692
print(Np,Nm)

#BSA = 1/Pb*(Np-Nm)/(Np+Nm)


BSA = 1/Pb*(Np-Nm)/(Nm+Np)



#eBSA = 2/Pb * np.sqrt(Np*Nm/(Np+Nm)**3)
eBSA = 1/Pb*np.sqrt(4*((Nm)**2/((Nm+Np)**4))*Np_e**2 + 4*((Np)**2/((Nm+Np)**4))*Nm_e**2)

if Np[1][0]== 0 and Nm[1][0]== 0:
    BSA[1][0] = 0
    BSA[1][8] = 0
    eBSA[1][0] = 0
    eBSA[1][8] = 0
    
if Np[2][0]== 0 and Nm[2][0]== 0:
    BSA[2][0] = 0
    BSA[2][8] = 0
    eBSA[2][0] = 0
    eBSA[2][8] = 0

print("BSA = ",BSA)
print("error BSA = ",eBSA)

trePhi = np.array([20,60,100,140,180,220,260,300,340])