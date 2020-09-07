
import bempp.api, numpy as np
from math import pi
import os
import time

from constants import values
from constants import mesh_info

from bempp.api.operators.potential import laplace as lp
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

from Grid_Maker_R2  import *
from Mesh_Ref_V2    import *
from quadrature     import *
from Potential_Solver import *
#from bem_parameters import *
from File_converter_python2 import *

global mol_name , mesh_density , suffix , path , q , x_q , phi_space , phi_order , u_space , u_order

def main_Adj(name , dens , input_suffix , output_suffix , percentaje , smooth = False , refine = True,
        Use_GAMer = False ):
    
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
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix , GAMer = Use_GAMer)
    else:
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix )
    
    init_time = time.time()
    
    dirichl_space_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order)
    neumann_space_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order) 
    dual_to_dir_s_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order)
    
    space_time = time.time() - init_time

    U, dU , operators_time , matrix_assembly_time , solv_time , U_it_count= U_tot(dirichl_space_u  , 
                                                                  neumann_space_u , dual_to_dir_s_u )
    
    S_trad = S_trad_calc_R( dirichl_space_u, neumann_space_u , U , dU )
    
    return S_trad

Resultados = open( 'Resultados/Resultados_Variando_S_Trad.txt' , 'w+' )

Resultados.write( 'molecule & density & Quad & Strad \n')

for molecule in ( 'methanol' , '1Br2Etano' ):
        
        for dens in ( 0.5 , 1.0 ):

            for i in np.arange(1,20):

                # bempp Parameters
                #print(bempp.api.global_parameters.hmat.coarsening)
                bempp.api.global_parameters.hmat.eps = 1e-8

                bempp.api.global_parameters.hmat.max_block_size = 2048
                bempp.api.global_parameters.hmat.min_block_size = 21
                bempp.api.global_parameters.hmat.max_rank = 30

                bempp.api.global_parameters.quadrature.double_singular = i ###

                bempp.api.global_parameters.quadrature.near.max_rel_dist = 2.0
                bempp.api.global_parameters.quadrature.near.single_order = i
                bempp.api.global_parameters.quadrature.near.double_order = i

                bempp.api.global_parameters.quadrature.medium.max_rel_dist = 4.0
                bempp.api.global_parameters.quadrature.medium.single_order = i
                bempp.api.global_parameters.quadrature.medium.double_order = i

                #bempp.api.global_parameters.quadrature.far.max_rel_dist = ?
                bempp.api.global_parameters.quadrature.far.single_order = i
                bempp.api.global_parameters.quadrature.far.double_order = i

                S_Trad = main_Adj( molecule , dens , '-0' , '-1' , 0.1,
                                                                  smooth=False , Use_GAMer= False)
                text = '{0} & {1} & {2:d}& {3:10.10f} \n'.format( molecule , str(dens) , i , S_Trad  )
                
                Resultados.write( text )
        
Resultados.write('Conditions\n')
Resultados.write('Smooth and Use_Gamer: Disabled \n')
Resultados.write('Adjoint mesh with %1% uniform refinement: Disabled \n')
Resultados.write('TESTING BEMPP QUADRATURE RULES PARAMETERS \n')
Resultados.close()
