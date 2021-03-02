
# 28.04.19
# New run for new mesh densities

# 10.10.19
# Now different times are measured, and N_it are also saved.
# Changed a potential function in Potential_Solver

# 11.09.19
# Added an adjoint mesh for S_Zeb


# 30.08.19
# Updates : Use of reaction potential instead of using the harmonic and regular component
# Vertices are exploded to the real geometry and the starting geometry is repaired with GAMer
# and built with MSMS


# RUN
# export LD_LIBRARY_PATH="/home/{user_name}/lib"
# In order to use FETK toolkit
# Be aware that the python shell won't give any warnings, but the cmd screen will.

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
    
    endt_time = time.time()
    fin_time = endt_time - init_time
    
    U_R , dU_R = U_Reac(U, dU , dirichl_space_u , neumann_space_u )
    
    ref_init_time = time.time()

    dirichl_space_phi = bempp.api.function_space(grid,  mesh_info.phi_space , mesh_info.phi_order)
    neumann_space_phi = bempp.api.function_space(grid,  mesh_info.phi_space , mesh_info.phi_order) 
    dual_to_dir_s_phi = bempp.api.function_space(grid,  mesh_info.phi_space , mesh_info.phi_order)
    
    phi , dphi , _ = adjoint_equation( dirichl_space_phi , neumann_space_phi , dual_to_dir_s_phi)

    aux_path = '_'+str(mesh_info.mesh_density)+ mesh_info.suffix

    face_array = np.transpose(grid.leaf_view.elements)+1
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    S_Cooper , S_Cooper_i = S_Cooper_calc( face_array , vert_array , phi , dphi , U_R , dU_R , 25)
    
    dirichl_space_u_s  = bempp.api.function_space(grid,  mesh_info.u_s_space , mesh_info.u_s_order )
    neumann_space_du_s = bempp.api.function_space(grid,  mesh_info.u_s_space , mesh_info.u_s_order )
 
    
    S_Zeb    , S_Zeb_i  , it_adj_count  = S_Zeb_in_Adjoint_Mesh(mesh_info.mol_name , face_array , vert_array ,
                                                                dens , input_suffix , 25)
    
    const_space = bempp.api.function_space(grid,  "DP", 0)

    Solv_i_Cooper_Func = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Cooper_i[:,0])
    Solv_i_Zeb_Func    = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Zeb_i[:,0])

    dif =S_Cooper_i-S_Zeb_i

    dif_F = bempp.api.GridFunction(const_space, fun=None, coefficients=np.abs(dif[:,0] ) )
    
    bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}.vtk'.format( dens,input_suffix )
                     , grid_function = dif_F , data_type = 'element')
    
    bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}.msh'.format( dens,input_suffix )
                     , grid_function = dif_F , data_type = 'element')
    
    if refine:
        new_face_array , new_vert_array = mesh_refiner(face_array , vert_array , np.abs(dif[:,0]) , percentaje )

        if smooth:
            fine_vert_array = text_to_list(name , '_40.0-0' , '.vert' , info_type=float)                
            new_vert_array = smoothing_vertex( new_vert_array , fine_vert_array )
            if Use_GAMer:
                face_array , vert_array = Improve_Mesh(new_face_array , new_vert_array , mesh_info.path  , 
                                                      mesh_info.mol_name+ '_' + str(dens) + output_suffix )
                new_face_array , new_vert_array = face_array.copy() , vert_array.copy()
                
        vert_and_face_arrays_to_text_and_mesh( name , new_vert_array ,
                                            new_face_array.astype(int)[:] , output_suffix, dens , Self_build=True)

        grid = Grid_loader( name , dens , output_suffix )
        print('New mesh:')
        #grid.plot()
    
    N_elements = len(face_array)
    
    ref_end_time = time.time()
    ref_fin_time = ref_end_time - ref_init_time
    
    return (S_trad , S_Cooper , S_Zeb , N_elements , fin_time , ref_fin_time , operators_time ,
            matrix_assembly_time , solv_time , it_adj_count , U_it_count)
    #return dif[:,0]
    

if True:
    Resultados = open( 'Resultados/Resultados_Adjoint_Mesh.txt' , 'w+' )

    Resultados.write( 'molecule & density  & Times-Refinated & Num of Elem  & Strad & SCooper & SZeb & calc_time '+
                    '& ref_time & op_time & Matrix_Time & Solv_Time & U_it_Count & Phi_it_Count \n')

    for molecule in ( 'methanol' , 'arg' ):
        
        for dens in ( 0.5 , 1.0 ):
            c_r = 0
            instance = ('-0' , '-1' , '-2', '-3', '-4', '-5' , '-6')#, '-7')
            for suffix in instance[:-1]:

                text = '{0} & {1} & {2}'.format( molecule , str(dens) , suffix  )

                (S_trad , S_Cooper , S_Zeb , N_elements , fin_time , ref_time ,
                 operators_time , matrix_assembly_time , solv_time ,
                 it_adj_count , U_it_count)  = main_Adj( molecule , dens , instance[c_r] , instance[c_r+1] , 0.1,
                                                                                  smooth=True , Use_GAMer= True)
                
                
                
                text = text + ' & {0:d} & {1:.10f} & {2:.10f} & {3:.10f} & {4:.10f} & {5:.10f}'.format( 
                    N_elements , S_trad , S_Cooper , S_Zeb , fin_time, ref_time  
                ) +' & {0:.10f} & {1:.10f} & {2:.10f} & {3:d} & {4:d} \n'.format(  
                        operators_time , matrix_assembly_time , solv_time , it_adj_count , U_it_count  )
                Resultados.write( text )
                c_r+=1
    Resultados.write('Conditions\n')
    Resultados.write('Smooth and Use_Gamer: Enabled \n')
    Resultados.write('Adjoint mesh with %1% uniform refinement: Enabled \n')
    Resultados.close()
