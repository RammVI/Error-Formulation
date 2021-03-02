
# Calculates the local error by using 2 different solutions U in DP-0 and u in P-1 from a uniformly refined mesh
# 21-07-2020


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

def main_Adj(name , dens , input_suffix , percentaje , N ,  N_ref , smooth = False , refine = True,
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

    U, dU = U_tot_boundary(grid , save_plot=True)
    
    S_trad = S_trad_calc_R( potential.dirichl_space_u, potential.neumann_space_u , U , dU )  

    S_Ap , S_Ap_i = delta_G_exact( grid , U , dU , mesh_info.u_space , mesh_info.u_order  , N ,
                                                    save_energy_plot=True                    )
    
    S_Ex , S_Ex_i , _ , N_El_adj = S_Exact_in_Adjoint_Mesh_with_N_Ref(name , grid , dens , input_suffix , N ,
                                            N_ref , save_energy_plot=True  , test_mode = True )

    const_space = bempp.api.function_space(grid,  "DP", 0)

    S_Ap_bempp = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ap_i[:,0])
    S_Ex_bempp    = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ex_i[:,0])

    dif =S_Ap_i-S_Ex_i

    dif_F = bempp.api.GridFunction(const_space, fun=None, coefficients=np.abs(dif[:,0] ) )
    
    bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_N_ref_{2:d}.vtk'.format( dens,
                                                        input_suffix , N_ref )
                     , grid_function = dif_F , data_type = 'element')
    
    if True:
        
        
        face_array = np.transpose(grid.leaf_view.elements) + 1        
        status = value_assignor_starter(face_array , np.abs(dif[:,0]) , percentaje)
        const_space = bempp.api.function_space(grid,  "DP", 0)
        Status    = bempp.api.GridFunction(const_space, fun=None, coefficients=status)
        bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_Marked_elements_{2}.vtk'.format( 
                                        dens, input_suffix , N_ref )
                     , grid_function = Status , data_type = 'element')
    
    if False: #Change for refine
        new_face_array , new_vert_array = mesh_refiner(face_array , vert_array , np.abs(dif[:,0]) , percentaje )

        if smooth:
            fine_vert_array = text_to_list(name , '_40.0-0' , '.vert' , info_type=float)                
            new_vert_array = smoothing_vertex( new_vert_array , fine_vert_array )
            if Use_GAMer:
                face_array , vert_array = Improve_Mesh(new_face_array , new_vert_array , path , 
                                                      mesh_info.mol_name+ '_' + str(dens) + output_suffix )
                new_face_array , new_vert_array = face_array.copy() , vert_array.copy()
                
        vert_and_face_arrays_to_text_and_mesh( name , new_vert_array ,
                                            new_face_array.astype(int)[:] , output_suffix, dens , Self_build=True)

        grid = Grid_loader( name , dens , output_suffix )
        #print('New mesh:')
        #grid.plot()
    
    N_elements = len(face_array)
    
    return ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj )
    #return dif[:,0]
    

txt_name = 'Resultados_22_07_2020.txt'
    
if True:
    Resultados = open( 'Resultados/'+txt_name , 'w+' )

    Resultados.write( ' density  & Times-Refinated & Elem & Adj Elem & S_Trad & S_Ap & S_Ex \n')
    
    Resultados.close()

    if True:
        molecule = 'methanol'
        
        for dens in [0.5 , 1.0 ]:
            
            for N_ref in [0 , 1 , 2 , 3 , 4]: 
                
                text = '{0} & {1}'.format( str(dens) , N_ref  )

                ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj )  = main_Adj( 
                    molecule , dens , '-0' , 0.1, 25 , N_ref , refine = False, 
                    smooth=False , Use_GAMer= False)
                
                Resultados = open( 'Resultados/'+txt_name , 'a' )                
                
                text = text + ' & {0:d} & {1:d} & {2:.10f} & {3:.10f} & {4:.10f} \n'.format( 
                     N_elements , N_El_adj  , S_trad , S_Ap , S_Ex  )
                                                                                        
                Resultados.write( text )
                
                Resultados.close()
                
    Resultados = open( 'Resultados/'+txt_name , 'a' ) 
    
    Resultados.write('Conditions\n')
    Resultados.write('Smooth and Use_Gamer: Disabled \n')
    Resultados.write('Adjoint mesh with uniform refinement: Enabled \n')
    Resultados.close()


