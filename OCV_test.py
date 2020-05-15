# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:48:34 2019

@author: Qiyang
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})


l = 1*1e-4 # electrolyte thickness in cm
n = 10000 # mesh density
temp = 800 # temperature in degree C
log_p_O2_ox = 0
log_p_O2_red = -21
x_L = np.linspace(0, 1, n)

def u_o2(log_p_O2, temp):
    return 8.617e-5*(temp+273)*np.log(float(10)**(log_p_O2))

def log_p_O2(u_O2, temp):
    return np.log10(np.exp(u_O2/(8.617e-5*(temp+273))))

def sigma_e(log_p_O2, temp):
    return 1.31 * 1e7 * np.exp (-3.88 /(8.617e-5*(temp+273))) * ((float(10)**(log_p_O2))**(-1/4))

def sigma_h(log_p_O2, temp):
    return 2.35 * 1e2 * np.exp (-1.67 /(8.617e-5*(temp+273))) * ((float(10)**(log_p_O2))**(1/4))

def sigma_el(log_p_O2, temp):
    sigma_e = 1.31 * 1e7 * np.exp (-3.88 /(8.617e-5*(temp+273))) * ((float(10)**(log_p_O2))**(-1/4))
    sigma_h = 2.35 * 1e2 * np.exp (-1.67 /(8.617e-5*(temp+273))) * ((float(10)**(log_p_O2))**(1/4))
    return sigma_e+sigma_h

def f(log_p_O2, j0):
    return 4*9.65e4*j0*(1/sigma_o+1/sigma_el(log_p_O2, temp))

sigma_o = 1.63e2 * np.exp (-0.79/(8.617e-5*(temp+273)))


j0 = 10**(-7)
j0_step = 0.1
u_O2 = np.zeros(n)
u1 = u_o2(log_p_O2_red, temp)
count = 0

while (np.abs(u1 - u_o2(log_p_O2_ox, temp)) > 0.1):
    count += 1
    u1 = u_o2(log_p_O2_red, temp)
    for i in range(n):
        u1 = f(log_p_O2(u1, temp), j0)*l/n + u1
        u_O2[i] = u1
    if u1 > u_o2(log_p_O2_ox, temp):
        j0 = j0 * (10**(-j0_step))
    else: 
        j0 = j0 * (10**(+j0_step))

# %%
figure1 = plt.figure(figsize = (6,7))
plt.plot(x_L, u_O2,'o', label = r'$\mu_{O2}$')
#plt.plot(u_O2/4,'o', label = r'$\tilde \mu_{e}$')
plt.ylabel(r'$\mu_{O2}$ or $\tilde \mu_{e}$')
plt.xlabel('x/L')
plt.legend()
plt.tight_layout()

# %%

figure2 = plt.figure(figsize = (6,7))
plt.semilogy(u_O2, sigma_e(log_p_O2(u_O2,temp), temp),'o', label = r'$\sigma_{e}$')
plt.semilogy(u_O2, sigma_h(log_p_O2(u_O2,temp), temp),'o', label = r'$\sigma_{h}$')
plt.semilogy(u_O2, sigma_el(log_p_O2(u_O2,temp), temp),'o', label = r'$\sigma_{e}+\sigma_{h}$')    

plt.ylabel(r'$\sigma_{el}$')
plt.xlabel('$\mu_{O2}$')
plt.legend()
plt.tight_layout()
    