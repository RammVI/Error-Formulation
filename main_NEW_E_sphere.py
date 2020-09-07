
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

def main_Adj(name , dens , input_suffix , output_suffix , percentaje ,  r , x_q , q , N ,  N_ref ,
             smooth = False , refine = True, Use_GAMer = False , sphere=True):
    
    mesh_info.mol_name     = name
    mesh_info.mesh_density = dens
    mesh_info.suffix       = input_suffix
    mesh_info.path         = os.path.join('Molecule' , mesh_info.mol_name)
    
    print('{1:.0f} {0} '.format(mesh_info.suffix , percentaje * 100.) )

    mesh_info.q , mesh_info.x_q = q , x_q

    mesh_info.u_space , mesh_info.u_order     = 'DP' , 0
    mesh_info.phi_space , mesh_info.phi_order = 'P' , 1
    mesh_info.u_s_space , mesh_info.u_s_order = 'P' , 1


    bempp.api.set_ipython_notebook_viewer()
    bempp.api.PLOT_BACKEND = "ipython_notebook"
       
    if input_suffix == '-0' and not sphere:
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix , GAMer = Use_GAMer)
    else:
        grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix )
    
    init_time = time.time()
    U, dU , operators_time_U , assembly_time_U , solving_time_U , it_count_U = U_tot_boundary(grid
                                                , return_time = True , save_plot=True)
    total_solving_time = time.time() - init_time
    
    ################
    
    init_time = time.time()
    S_trad = S_trad_calc_R( potential.dirichl_space_u, potential.neumann_space_u , U , dU ) 
    S_trad_time = time.time()-init_time
    
    ################
    
    init_time = time.time()
    S_Ap , S_Ap_i = delta_G_exact( grid , U , dU , mesh_info.u_space , mesh_info.u_order  , N ,
                                                    save_energy_plot=True                    )
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
    
    if True:
        
        
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

            aux_vert_array = np.zeros(( len(new_vert_array),3 ))

            c=0
            for vert in new_vert_array:
                aux_vert_array[c] = smothing_func(vert , r) 
                c+=1
            
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
    

txt_name = '_Resultados_Uniform.txt'
percentaje , s =2.5 , True
N_it     = 5

dens = 0

name = 'off_centered_nref_0'

r = 1.
x_q = np.array( [[  1.E-12 ,  1.E-12 ,   r/2. ]]  )
q = np.array( [1.] )
#x_q = np.array( [[  1.E-12 ,  1.E-12 ,  0.62 ],
#                 [  1.E-12 ,  0.62*np.cos(np.pi*1.5 + 5.*np.pi/180.) ,
#                                                      0.62*np.sin(np.pi*1.5 + 5.*np.pi/180. ) ] ,
#                 [  1.E-12 ,  0.62*np.cos(np.pi*1.5 - 5.*np.pi/180.) ,
#                                                      0.62*np.sin(np.pi*1.5 - 5.*np.pi/180. )  ]
#                       ] )
#q = np.array( [1. , 1. , -1.]) 



results_name = 'Resultados/'+ name + txt_name 
    
if True:
    Resultados = open( results_name , 'w+' )
    
    a   = 1.
    N   = 20
    
    f_ex = solution(q, x_q, ep_m, ep_s, r, k , a, N)
    print('Exact Solution: {0:5f}--------------------------------------------'.format(f_ex))
    Resultados.write('Exact solution: {0:10f} \n'.format(f_ex))

    Resultados.write( ' density  & Times-Refinated & Elem & Adj Elem & S_Trad & S_Ap & S_Ex &' 
                    + ' solving_time & S_trad_time & S_Ap_time & S_Ex_time & Operators_Build_Time_U & Assembly_Time_U & Linear_Sys_Solv_time_U & Iterations_U \n')
    
    Resultados.close()
    
    path = os.path.join( 'Molecule' , name , str( int(percentaje*100)) , str(s) )            
    os.system('mkdir -p {0}'.format( path ) )
    
    dec_v_array , dec_f_array = decahedron_points(r)
    vert_and_face_arrays_to_text_and_mesh( name, dec_v_array , dec_f_array+1 , '-s0' 
                                                      , dens=0 , Self_build=True)
    # Refines to obtain a better aproach of a sphere
    main_Adj(name , dens , '-s0' , '-s1' , percentaje=2.5 , r=r , x_q=x_q , q=q , N=25, N_ref=3 , smooth = True ,
               refine = True, Use_GAMer = True)
    
    if True:
        
        
        if True:
                         
            suffixes = suffix_names(N_it)
                        
            c = 1
            for i in suffixes[1:-1]:
                
                text = '{0} & {1}'.format( str(dens) , 0  )

                ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , tot_solving_time , S_trad_time , S_Ap_time ,
                  S_Ex_time , op_time_U , assem_time_U , solv_time_U , it_count_U )  = main_Adj( 
                    name , dens , suffixes[c]  , suffixes[c+1] , percentaje , r=r , x_q=x_q , q=q , 
                    N=25 , N_ref=0, smooth = True , refine = True, Use_GAMer = True , sphere=True)
                
                Resultados = open( results_name , 'a' )                
                
                text = text + ' & {0:d} & {1:d} & {2:.10f} & {3:.10f} & {4:.10f} & {5:.10f} & {6:.10f} & {7:.10f} & {8:.10f} & {9:.10f} & {10:.10f} & {11:.10f} & {12:d} \n'.format( 
                     N_elements , N_El_adj  , S_trad , S_Ap , S_Ex , tot_solving_time , S_trad_time , S_Ap_time , S_Ex_time , op_time_U , assem_time_U , solv_time_U , it_count_U)
                                                                                          
                Resultados.write( text )
                
                c+=1                
                Resultados.close()
                
    Resultados = open( results_name , 'a' ) 
    
    Resultados.write('Conditions\n')
    Resultados.write('Smooth and Use_Gamer: Enabled \n')
    Resultados.write('Adjoint mesh with uniform refinement: Enabled \n')
    Resultados.write('1 refinement on phi grid \n')
    Resultados.close()

    moving_files(name , percentaje , s)
