
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
from sphere_geom    import *
#from bem_parameters import *
from File_converter_python2 import *

from analytical     import *

global mol_name , mesh_density , suffix , path , q , x_q , phi_space , phi_order , u_space , u_order

def main_Adj(name , dens , input_suffix , output_suffix , percentaje ,  N ,  N_ref ,
             smooth = True , refine = True, Use_GAMer = True , sphere=False):
    
    mesh_info.mol_name     = name
    mesh_info.mesh_density = dens
    mesh_info.suffix       = input_suffix
    mesh_info.path         = os.path.join('Molecule' , mesh_info.mol_name)

    mesh_info.q , mesh_info.x_q = run_pqr(mesh_info.mol_name)

    mesh_info.u_space , mesh_info.u_order     = 'DP' , 0
    mesh_info.phi_space , mesh_info.phi_order = 'P' , 1
    mesh_info.u_s_space , mesh_info.u_s_order = 'P' , 1
       
    if input_suffix == '-0' and not sphere:
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix , GAMer = Use_GAMer)
    elif input_suffix!='-0':
        print('Loading previus mesh')
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix )
        
    face_array = np.transpose(grid.leaf_view.elements)
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    init_time = time.time()
    U, dU , operators_time_U , assembly_time_U , solving_time_U , it_count_U = U_tot_boundary(grid
                                                , return_time = True , save_plot=True, tolerance = 1e-8)
    total_solving_time = time.time() - init_time
    ################
    
    phi , dphi , adj_grid = phi_with_N_ref(name , grid , face_array , vert_array ,
                    dens , input_suffix , N_ref , return_grid = True)
    
    U_R , dU_R = U_Reac( U, dU , potential.dirichl_space_u , potential.neumann_space_u )
    
    ################
    
    init_time = time.time()
    S_trad = S_trad_calc_R( potential.dirichl_space_u, potential.neumann_space_u , U , dU ) 
    S_trad_time = time.time()-init_time
    
    ################
    
    adj_face_array = np.transpose(adj_grid.leaf_view.elements)
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)
    
    init_time = time.time()
    S_Ap , S_Ap_i , relation = Aproximated_Sol_Adj_UDP0( U_R , dU_R , phi , dphi , face_array , vert_array , 
                             adj_face_array , adj_vert_array , 1 , grid , adj_grid , N_ref ,
                                                return_relation=True)
    S_Ap_time = time.time() - init_time
    
    #################
    
    init_time = time.time()
    S_Ex , S_Ex_i , _ , N_El_adj = S_Exact_in_Adjoint_Mesh_with_N_Ref(name , grid , dens , input_suffix , N ,
                                            N_ref , save_energy_plot=True  , test_mode = True )
    S_Ex_time = time.time() - init_time

    ################
    
    const_space = bempp.api.function_space(grid,  "DP", 0)

    S_Ap_bempp = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ap_i[:,0])
    S_Ex_bempp    = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ex_i[:,0])

    dif =S_Ap_i-S_Ex_i

    dif_F = bempp.api.GridFunction(const_space, fun=None, coefficients=np.abs(dif[:,0] ) )
    
    bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_N_ref_{2:d}.vtk'.format( dens,
                                                        input_suffix , N_ref )
                     , grid_function = dif_F , data_type = 'element')
    
    if True: #Marked elements        
        face_array = np.transpose(grid.leaf_view.elements) + 1        
        status = value_assignor_starter(face_array , np.abs(dif[:,0]) , percentaje)
        const_space = bempp.api.function_space(grid,  "DP", 0)
        Status    = bempp.api.GridFunction(const_space, fun=None, coefficients=status)
        bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_Marked_elements_{2}.vtk'.format( 
                                        dens, input_suffix , N_ref )
                     , grid_function = Status , data_type = 'element')
    
    face_array = np.transpose(grid.leaf_view.elements)+1
    vert_array = np.transpose(grid.leaf_view.vertices)
            
    N_elements = len(face_array)
    
    if refine:

        new_face_array , new_vert_array = mesh_refiner(face_array , vert_array , np.abs(dif[:,0]) , percentaje )

        if smooth:

            fine_vert_array = text_to_list(name , '_40.0-0' , '.vert' , info_type=float)                
            aux_vert_array  = smoothing_vertex( new_vert_array , fine_vert_array )
            
        elif not smooth:
        
            aux_vert_array = new_vert_array.copy()
        
        if Use_GAMer:
            
            new_face_array , aux_vert_array = Improve_Mesh(new_face_array , aux_vert_array , mesh_info.path , 
                                                          mesh_info.mol_name+ '_' + str(dens) + output_suffix )

            vert_and_face_arrays_to_text_and_mesh( name , aux_vert_array , new_face_array.astype(int)[:] ,
                                                   output_suffix, dens , Self_build=True)

            grid = Grid_loader( name , dens , output_suffix )
    
    return ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , total_solving_time , S_trad_time , S_Ap_time ,
             S_Ex_time , operators_time_U , assembly_time_U , solving_time_U , it_count_U )
    #return dif[:,0]
    

txt_name = 'Resultados_27_08_2020_methanol_E_phi_3_Adaptative.txt'
percentaje , s = 0.1 , True
N_it     = 10
    
if True:
    Resultados = open( 'Resultados/'+txt_name , 'w+' )

    Resultados.write( 'Molecule & density  & Times-Refinated & Elem & Adj Elem & S_Trad & S_Ap & S_Ex ' 
                    + ' solving_time & S_trad_time & S_Ap_time & S_Ex_time & Operators_Build_Time_U & Assembly_Time_U & Linear_Sys_Solv_time_U & Iterations_U \n')
    
    Resultados.close()    

    if True:
        molecule = 'methanol_E_phi'
        
        for dens in [0.5,1.0]:
                         
            suffixes = np.array(['-0','-1','-2','-3','-4','-5','-6','-7'])
                        
            c = 0
            for i in suffixes[:-1]:
                
                text = '{0} & {1:f} & {2:d}'.format( molecule, dens , 1  )

                ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , tot_solving_time , S_trad_time , S_Ap_time ,
                  S_Ex_time , op_time_U , assem_time_U , solv_time_U , it_count_U )  = main_Adj( 
                    molecule , dens , suffixes[c]  , suffixes[c+1] , percentaje ,
                    N=25 , N_ref=3 , smooth = True , refine = True, Use_GAMer = True)
                
                Resultados = open( 'Resultados/'+txt_name , 'a' )                
                
                text = text + ' & {0:d} & {1:d} & {2:.10f} & {3:.10f} & {4:.10f} & {5:.10f} & {6:.10f} & {7:.10f} & {8:.10f} & {9:.10f} & {10:.10f} & {11:.10f} & {12:.10f} \n'.format( 
                     N_elements , N_El_adj  , S_trad , S_Ap , S_Ex , tot_solving_time , S_trad_time , S_Ap_time , S_Ex_time , op_time_U , assem_time_U , solv_time_U , it_count_U)
                                                                                        
                Resultados.write( text )
                
                c+=1                
                Resultados.close()
                
    Resultados = open( 'Resultados/'+txt_name , 'a' ) 
    
    Resultados.write('Conditions\n')
    Resultados.write('Smooth and Use_Gamer: Disabled \n')
    Resultados.write('3 refinements on phi!      \n')
    Resultados.write('Using E_phi!!!!!!!   \n')
    Resultados.close()
