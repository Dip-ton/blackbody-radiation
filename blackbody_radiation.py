#Planck's Spectral Distribution using Quantization of Energy: E=hnu

import os
import numpy as np
import matplotlib.pyplot as plt

TT=[1000,2000,3000,4000,5000]
nu=np.linspace (1e10,3e15,100000)

h=6.62607015e-34
c=299792458
k=1.380649e-23

T_list=[]
energy_list=[]
numax_list=[]

plt.figure()
for T in TT:
    T_list.append(T)
    I=(2*h*(nu**3))/(c**2)*(1/(np.exp(h*nu/(k*T))-1))
    area=np.trapz(I,nu)
    power=np.pi*area
    energy_list.append(power)
    idx=np.argmax(I)
    nu_max=nu[idx]
    numax_list.append(nu_max)
    plt.plot (nu,I, label=f'T={T}K')

plt.xlim(0,1.5e15)
plt.xlabel('Frequency(Hz)--->')
plt.ylabel('Spectral Radiance (W.m^-2.sr^-1.Hz^-1)--->')
plt.title("Planck's Spectral Distribution")
plt.grid (True)
plt. legend()
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")
os.makedirs(out_dir, exist_ok=True)
plt.savefig(os.path.join(out_dir, "planck_spectral_distribution.png"), dpi=300, bbox_inches="tight")
plt.show()

print ("Verification of Stefan-Boltzmann Law")
for i in range (len(T_list)):
    value=energy_list[i]/(T_list[i])**4
    print (f"T={T_list[i]}K,Stefan-Boltzmann Constant={value:.3e}")
plt.figure()
plt.loglog(T_list,energy_list,'o-')
plt.xlabel('Temperature (K)--->')
plt.ylabel('Total Energy (W.m^-2)--->')
plt.title("Verification of Stefan-Boltzmann law")
plt.grid (True)
plt.savefig(os.path.join(out_dir, "stefan_boltzmann_law.png"), dpi=300, bbox_inches="tight")
plt.show()
slope,intercept=np.polyfit(np.log(T_list),np.log(energy_list),1)
print ("Stefan-Boltzmann Exponent=",slope)
print ("Stefan-Boltzmann Constant=",np.exp(intercept))

print (f"Verification of Wien's Law")
for i in range (len(T_list)):
    value=numax_list[i]/T_list[i]
    print (f"T={T_list[i]}K,nu_max/T={value:.3e}")
plt.figure()
plt.plot(T_list,numax_list,'o-')
plt.xlabel('Temperature (K)--->')
plt.ylabel('Maximum Frequency (Hz)--->')
plt.title("Verification of Wien's Displacement Law")
plt.grid (True)
plt.savefig(os.path.join(out_dir, "wiens_displacement_law.png"), dpi=300, bbox_inches="tight")
plt.show()
slope,intercept=np.polyfit(T_list,numax_list,1)
print (f"Wien's displacement constant={slope}")