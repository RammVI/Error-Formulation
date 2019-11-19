
# 15.11.19
# Added a function to calculate volume integrals via Monte Carlo integration
# THIS RUTINE WILL DO ONLY MC INTREGRALS

# Acumulative 10.10.19 upgrades

import bempp.api, numpy as np
from math import pi
import os
import time

from constants import values
from constants import mesh_info
from constants import potential

from bempp.api.operators.potential import laplace as lp
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

from Grid_Maker_R2  import *
from Mesh_Ref_V2    import *
from quadrature     import *
from Potential_Solver import *
#from bem_parameters import *
from File_converter_python2 import *
import trimesh

global mol_name , mesh_density , suffix , path , q , x_q , phi_space , phi_order , u_space , u_order

def main_MC(name , dens , input_suffix ):
    
    mesh_info.mol_name     = name
    mesh_info.mesh_density = dens
    mesh_info.suffix       = input_suffix
    mesh_info.path         = os.path.join('Molecule' , mesh_info.mol_name)

    mesh_info.q , mesh_info.x_q = run_pqr(mesh_info.mol_name)

    mesh_info.u_space , mesh_info.u_order     = 'DP' , 0
    mesh_info.phi_space , mesh_info.phi_order = 'P' , 1
    mesh_info.u_s_space , mesh_info.u_s_order = 'P' , 1


    bempp.api.set_ipython_notebook_viewer()
    bempp.api.PLOT_BACKEND = "ipython_notebook"
       
    if input_suffix == '-0':
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix , GAMer = False)
    else:
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix )
    
    potential.dirichl_space_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order)
    potential.neumann_space_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order) 
    potential.dual_to_dir_s_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order)

    potential.U, potential.dU , operators_time , matrix_assembly_time , solv_time , U_it_count= U_tot(
        potential.dirichl_space_u  , potential.neumann_space_u , potential.dual_to_dir_s_u )
    
    U_R , dU_R = U_Reac(potential.U, potential.dU , potential.dirichl_space_u , potential.neumann_space_u )

    face_array = np.transpose(grid.leaf_view.elements)+1
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    mesh = trimesh.Trimesh(vertices = vert_array , faces = face_array-1)
    
    depreciated_term , used_points = scalar_times_laplacian(mesh , local_U_interior , 
                           local_U_Reac_interior , 40 , 0.001)   
    
    
    return depreciated_term , used_points
    #return dif[:,0]
    

if True:
    Resultados = open( 'Resultados_MC_15_11.txt' , 'w+' )

    Resultados.write( ' molecule & Density & Vol integral & N Points \n')

    for molecule in ('arg' , 'methanol'):
        
        for dens in ( 0.5 , 1.0 , 2.0 ):
            
            text = '{0} & {1}'.format( molecule , str(dens)  )

            depreciated_term , used_points = main_MC(molecule , dens , '-0' )
                
                
            text = text + ' & {0:.10f} & {1:d}'.format( depreciated_term , used_points )
            Resultados.write( text )
                            
    Resultados.write('Conditions\n')
    Resultados.write('Smooth and Use_Gamer: both Disabled \n')
    Resultados.write('Adjoint mesh with uniform refinement: Disabled \n')
    Resultados.write('Note: Used for Monte Carlo integration, but points are not strictly random, \n')
    Resultados.write('      points distance is constant in each direction, but not controlable. See quadrature \n')
    Resultados.write('      for more information \n')
    Resultados.close()
