
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
    
    phi , dphi , adj_grid = phi_with_N_ref(name , grid , face_array , vert_array ,
                    dens , input_suffix , N_ref , return_grid = True)
    
    adj_face_array = np.transpose(adj_grid.leaf_view.elements)
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)
    
    S_Aprox, S_Aprox_i , relation = Aproximated_Sol_Adj_UDP0(U_R , dU_R , phi , dphi , face_array , vert_array , 
                             adj_face_array , adj_vert_array , N , grid , adj_grid , N_ref ,
                                                return_relation=True)    
    
    print(relation , N_ref)
    
    S_Ex , S_Ex_i = Zeb_aproach_with_u_s_Teo( adj_face_array , adj_vert_array , phi , dphi , N ,
                                            grid_relation=relation )
    
    S_trad = S_trad_calc_R( potential.dirichl_space_phi, potential.neumann_space_phi , 
                            potential.phi, potential.dphi )
    
    const_space = bempp.api.function_space(grid,  "DP", 0)

    S_Aprox_i_BEMPP = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Aprox_i[:,0])
    S_Ex_i_BEMPP    = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ex_i[:,0])
    
    dif =S_Aprox_i_BEMPP-S_Ex_i_BEMPP

    dif_F = bempp.api.GridFunction(const_space, fun=None, coefficients=np.abs(dif.coefficients.real) )   
    
    bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_ADJ_ESP_{2:d}.vtk'.format( dens,input_suffix , N_ref)
                     , grid_function = dif_F , data_type = 'element')
    
    N_el = len(adj_face_array)
    
    return N_el , dif.coefficients.real , S_trad , S_Aprox , S_Ex 

if True:
    Resultados = open( 'Resultados/Resultados_meth_16_04.txt' , 'w+' )

    Resultados.write( ' molecule & Density & N_ref & N_El_adj & Trad Solvation & Aprox Solvation & Exact Solvation \n')

    for molecule in ['1Br2Etano' , 'methanol' ]:
        
        for dens in [0.5,1.0]:

            for N_ref in [0 , 1 , 2 , 3]:

                text = '{0} & {1} & {2} '.format( molecule , str(dens) , str(N_ref))

                N_el , Dif , S_T , S_A , S_E = main_MC( molecule , dens , '-0' , 25 , N_ref)

                text = text + '& {3:d} & {0:.10f} & {1:.10f} & {2:.10f} \n'.format( S_T , S_A , S_E , N_el)
                Resultados.write( text )
                            
    Resultados.write('Conditions: mesh density = {0} \n'.format(dens))
    Resultados.write('Smooth and Use_Gamer: both Disabled \n')
    Resultados.write('Adjoint mesh with uniform refinement: Enabled with *4* refinements on the Adjoint mesh \n')
    Resultados.write('Notes:  \n')
    Resultados.write('      U is from DP-0 and phi from P-1 \n')
    
    Resultados.close()
