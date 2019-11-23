import os
import numpy as np
#global mol_name , mesh_density , suffix , path , q , x_q , phi_space , phi_order , u_space , u_order

def values():
    mol_name = ''
    mesh_density = 2.0
    suffix   = '-0'
    path = os.path.join('Molecule',mol_name) 

    stern_thickness = 0

    q   = np.array((1.))
    x_q = np.array((0. , 0. , 0. ))

    phi_space , phi_order = 'P'  , 1
    u_space   , u_order   = 'DP' , 0
    
