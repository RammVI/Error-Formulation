
# 25.01.2020
# Routine that uses 4 uniform refinement for the adjoint
# Tested only for the 2.0 density for methanol

# Acumulative 10.10.19 upgrades

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
    
    
    potential.U, potential.dU = U_tot_boundary(grid)
    
    
    
    U_R , dU_R = U_Reac(potential.U, potential.dU , potential.dirichl_space_u , potential.neumann_space_u )
    
    U_c , dU_c = potential.U - U_R , potential.dU - dU_R
    
    phi , dphi , adj_grid = phi_with_N_ref(name , grid , face_array , vert_array ,
                    dens , input_suffix , N_ref , return_grid = True)
    
    phi_c  = bempp.api.GridFunction(potential.dirichl_space_phi, fun=u_s_G)
    dphi_c = bempp.api.GridFunction(potential.neumann_space_phi, fun=du_s_G)
    
    phi_r , dphi_r = phi - phi_c , dphi - dphi_c
    
    adj_face_array = np.transpose(adj_grid.leaf_view.elements)
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)
    
    S_Aprox   , _ , relation = Aproximated_Sol_Adj_UDP0( U_R , dU_R , phi , dphi , face_array , vert_array , 
                             adj_face_array , adj_vert_array , 1 , grid , adj_grid , N_ref ,
                                                return_relation=True)   
    
    S_Aprox_E1, _ , relation = Aproximated_Sol_Adj_UDP0( U_R , dU_R , phi_r , dphi_r , face_array , vert_array , 
                             adj_face_array , adj_vert_array , 1 , grid , adj_grid , N_ref ,
                                                return_relation=True)   
    
    S_Aprox_E2, _ , relation = Aproximated_Sol_Adj_UDP0( U_c , dU_c , phi_c , dphi_c , face_array , vert_array , 
                             adj_face_array , adj_vert_array , 1 , grid , adj_grid , N_ref ,
                                                return_relation=True)   
    
    S_trad = S_trad_calc_R( potential.dirichl_space_u, potential.neumann_space_u , 
                           potential.U,  potential.dU ) 
    #print(relation , N_ref)
    
    ### TEST MODE
    
    #const_space = bempp.api.function_space(grid,  "DP", 0)
    #const_space_ADJ = bempp.api.function_space(adj_grid,  "DP", 0)
    
    #rel_ADJ = bempp.api.GridFunction(const_space_ADJ, fun=None, coefficients= relation )   
    
    #bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_ADJ_ESP_{2:d}_relation.vtk'.format(
    #                                                                    dens,input_suffix , N_ref)
    #                 , grid_function = rel_ADJ , data_type = 'element')
    
    #rel = bempp.api.GridFunction(const_space, fun=None, coefficients= np.arange(0,len(face_array)) )   
    
    #bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_ADJ_ESP_{2:d}_original_relation.vtk'.format(
    #                                                                    dens,input_suffix , N_ref)
    #                 , grid_function = rel , data_type = 'element')
    
    
    ####
    
    
    #S_Ex , S_Ex_i = Zeb_aproach_with_u_s_Teo( adj_face_array , adj_vert_array , phi , dphi , N ,
    #                                        grid_relation=relation )
    
    #const_space = bempp.api.function_space(grid,  "DP", 0)

    #S_Aprox_i_BEMPP = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Aprox_i[:,0])
    #S_Ex_i_BEMPP    = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ex_i[:,0])
    
    #dif =S_Aprox_i_BEMPP-S_Ex_i_BEMPP

    #dif_F = bempp.api.GridFunction(const_space, fun=None, coefficients=np.abs(dif.coefficients.real) )   
    
    #bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_ADJ_ESP_{2:d}_RRR.vtk'.format( dens,input_suffix , N_ref)
    #                 , grid_function = S_Aprox_i_BEMPP , data_type = 'element')
    
    N_el = len(adj_face_array)
    
    ### NEW, ONLY FOR S_TRAD CALCULATED IN THE NEW MESH - REMEMBER THAT PHI EQUALS u (potential)
    
    #mesh_info.phi_space , mesh_info.phi_order = 'DP' , 0
    #mesh_info.u_s_space , mesh_info.u_s_order = 'DP' , 0
    
    #phi , dphi , adj_grid = phi_with_N_ref(name , grid , face_array , vert_array ,
    #                dens , input_suffix , N_ref , return_grid = True)
    
    #S_trad = S_trad_calc_R( potential.dirichl_space_phi, potential.neumann_space_phi , 
    #                       potential.phi,  potential.dphi ) 
    
    return N_el , S_trad , S_Aprox , S_Aprox_E1 , S_Aprox_E2 # , dif.coefficients.real  , S_trad , S_Aprox , S_Ex

if True:
    Resultados = open( 'Resultados/Resultados_meth_19_07_2020.txt' , 'w+' )

    Resultados.write( ' molecule & Density & N_ref & N_El_adj & S_T & S_Aprox & S_A_E1 & S_A_E2 \n')

    for molecule in ['1Br2Etano' , 'methanol' ]:
    
        #for dens in [0.5,1.0]:
        if True:
            dens = 0.5
            
            for N_ref in [0 , 1 , 2]:#, 3]:

                text = '{0} & {1} & {2} '.format( molecule , str(dens) , str(N_ref))

                N_el , S_T , S_A , S_A_E1 , S_A_E2 = main_MC( molecule , dens , '-0' , 25 , N_ref)

                text = text + '& {0:d} & {1:.10f} & {2:.10f} & {3:.10f} & {4:.10f} \n'.format( N_el , S_T , S_A ,
                                                                    S_A_E1 , S_A_E2)
                Resultados.write( text )
                            
    Resultados.write('Conditions: Special integrals case \n'.format(dens))
    Resultados.write('Smooth and Use_Gamer: both Disabled \n')
    Resultados.write('Adjoint mesh with uniform refinement: Enabled with *4* refinements on the Adjoint mesh \n')
    Resultados.write('Notes:  \n')
    Resultados.write('      U is from DP-0 and phi from p-1 \n')
    
    Resultados.close()
