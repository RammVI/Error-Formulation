import os
import numpy as np

class values():
    
    global ep_m , ep_s , k , e_c , k_B , T , C , K
    
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
    
class mesh_info():
    
    global mol_name , mesh_density , suffix , path , q , x_q , phi_space , phi_order , u_space , u_order
    
    mol_name = ''
    mesh_density = 2.0
    suffix   = '-0'
    path = os.path.join('Molecule',mol_name) 
    
    stern_thickness = 0
    
    q   = np.array((1.))
    x_q = np.array((0. , 0. , 0. ))
    
    phi_space , phi_order = 'P'  , 1
    u_space   , u_order   = 'DP' , 0
    

    
def run_pqr(mol_name):
    
    global q , x_q
    
    path = os.path.join('Molecule',mol_name) 
    
    q, x_q = np.empty(0), np.empty((0,3))
    pqr_file = os.path.join(path,mol_name+'.pqr')
    charges_file = open( pqr_file , 'r').read().split('\n')

    for line in charges_file:
        line = line.split()
        if len(line)==0: continue
        if line[0]!='ATOM': continue
        q = np.append( q, float(line[8]))
        x_q = np.vstack( ( x_q, np.array(line[5:8]).astype(float) ) )  

    return q , x_q