
# RUTINE TO CALCULATE THE SOLVATION ENERGY TRAD, WITH U IN DP-0 AND WITH EXPONENTIAL FUNCTIONS

import bempp.api, numpy as np
from math import pi
import os
import time

from constants import values, mesh_info , potential

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

def fx(x,n,domain_index,result):
    result[:] =  x[0]
    
def fx2(x,n,domain_index,result):
    result[:] =  x[0]**2.

def main_MC(name , dens , input_suffix , N , N_ref):
    
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
    
    face_array = np.transpose(grid.leaf_view.elements)
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    lineal = np.arange(0 , len(face_array) , 1)
    
    
    #potential.U, potential.dU = U_tot_boundary(grid)
    
    const_space = bempp.api.function_space(grid,  "DP", 0)
    
    potential.U, potential.dU = U_tot_boundary(grid)
    
    U_R , dU_R = U_Reac(potential.U, potential.dU , potential.dirichl_space_u , potential.neumann_space_u )
    
    #U_R_c , dU_R_c = np.ones([len(face_array)])*0. ,  np.ones([len(face_array)])#bempp.api.GridFunction( const_space , fun =  fx)
    
    #U_R = bempp.api.GridFunction(const_space, fun=None, coefficients= U_R_c )  
    #dU_R= bempp.api.GridFunction(const_space, fun=None, coefficients= dU_R_c)#.coefficients.real  )  
    
    phi , dphi , adj_grid = phi_with_N_ref(name , grid , face_array , vert_array ,
                    dens , input_suffix , N_ref , return_grid = True , calculate_phi = True )
    
    space_ADJ = bempp.api.function_space(adj_grid,  "P", 1)
    
    dirichl_space = bempp.api.function_space(adj_grid,  'P' , 1)
    
    u_s  = bempp.api.GridFunction(dirichl_space , fun =  u_s_G)
    du_s = bempp.api.GridFunction(dirichl_space , fun = du_s_G)
    
    adj_face_array = np.transpose(adj_grid.leaf_view.elements)
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)
    
    #phi_c , dphi_c = bempp.api.GridFunction( space_ADJ , fun =  fx2) , np.ones([len(adj_vert_array)])
    
    #phi = bempp.api.GridFunction(space_ADJ, fun=None, coefficients= phi_c.coefficients.real )  
    
    # 
    #dphi= bempp.api.GridFunction(space_ADJ, fun=None, coefficients= dphi_c ) 
        
    
    
    S_Aprox, S_Aprox_i , relation = Aproximated_Sol_Adj_UDP0( U_R , dU_R , u_s , du_s ,
                                                             face_array , vert_array , 
                             adj_face_array , adj_vert_array , N , grid , adj_grid , N_ref ,
                                                return_relation=True)    
    
    
    #S_trad = S_trad_calc_R( potential.dirichl_space_phi, potential.neumann_space_phi , 
    #                        potential.phi,  potential.dphi ) 
    
    N_el_adj = len(adj_face_array)
    
    
    return N_el_adj , S_Aprox  #, S_trad

if True:
    Resultados = open( 'Resultados/Resultados_meth_3_13_07_2020.txt' , 'w+' )

    Resultados.write( ' molecule & Density & N_ref & N_El_adj & S App \n')

    #for molecule in ['1Br2Etano' , 'methanol' ]:
    if True:
        molecule = '1Br2Etano'
        #for dens in [0.5,1.0]:
        if True:
            dens = 2.0
            
            for N_ref in [0 , 1 ,3 ]:

                text = '{0} & {1} & {2} '.format( molecule , str(dens) , str(N_ref))

                N_el , S_app  = main_MC( molecule , dens , '-0' , 25 , N_ref)

                text = text + '& {0:d} & {1:.10f}  \n'.format( N_el , S_app )
                Resultados.write( text )
                            
    Resultados.write('Conditions: - \n'.format(dens))
    Resultados.write('Smooth and Use_Gamer: both Disabled \n')
    Resultados.write('Adjoint mesh with uniform refinement: Enabled with *3* consecutive refinements on the Adjoint mesh \n')
    Resultados.write('Notes:  \n')
    Resultados.write('       INTEGRATE OVER GAMMA FOR S APROX\n')
    Resultados.write('       !!!!!!!! \n')
    
    Resultados.close()
