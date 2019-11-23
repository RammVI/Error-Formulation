
'''
Saves values, constants, mesh info and obtained solutions from others .py files
'''

import os
import numpy as np

K     = 0.5 * 4. * np.pi * 332.064 
        
ep_m = 4.
ep_s = 80.
k = 0.125

# e_c : Carga del proton
# k_B : Constante de Stefan Boltzmann
# T   : Temperatura promedio, en grados Kelvin
# C   : Constante igual a e_c/(k_B*T). Para el caso se utilizara como 1 y se agregara una cte a la QoI
e_c = 1.60217662e-19 # [C] - proton charge
k_B = 1.38064852e-23 # [m2 kg s-2 K-1]
T   = 298. # [K]  
C = 1. # 
