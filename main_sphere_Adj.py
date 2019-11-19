
# October 10 - 2019
# This rutine uses the three term splitting formulation

# Same as the rutine defined above, uses an uniform refined mesh for the Adjoint

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


#global mol_name , mesh_density , suffix , path , q , x_q , phi_space , phi_order , u_space , u_order

def main(name , dens , input_suffix , output_suffix , percentaje , r , x_q
         , q , smooth = False , refine = True  ):
    
    init_time = time.time()

    mesh_info.mol_name     = name
    mesh_info.mesh_density = dens
    mesh_info.suffix       = input_suffix
    mesh_info.path         = os.path.join('Molecule' , mesh_info.mol_name)
    
    print('{1:.0f} {0} '.format(mesh_info.suffix , percentaje * 100.) )

    mesh_info.q , mesh_info.x_q = q , x_q

    mesh_info.u_space , mesh_info.u_order     = 'DP' , 0
    mesh_info.phi_space , mesh_info.phi_order =  'P' , 1
    mesh_info.u_s_space , mesh_info.u_s_order =  'P' , 1


    bempp.api.set_ipython_notebook_viewer()
    bempp.api.PLOT_BACKEND = "ipython_notebook"

    grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix )
        
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

    u_s  = bempp.api.GridFunction(dirichl_space_u_s , fun=u_s_G )
    du_s = bempp.api.GridFunction(neumann_space_du_s, fun=du_s_G)      

    S_Zeb    , S_Zeb_i  , it_adj_count  = S_Zeb_in_Adjoint_Mesh(mesh_info.mol_name , face_array , vert_array , dens , input_suffix , 25)

    const_space = bempp.api.function_space(grid,  "DP", 0)

    Solv_i_Cooper_Func = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Cooper_i[:,0])
    Solv_i_Zeb_Func    = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Zeb_i[:,0])

    dif =S_Cooper_i-S_Zeb_i

    dif_F = bempp.api.GridFunction(const_space, fun=None, coefficients=np.abs(dif[:,0] ) )
    #dif_F.plot()    

    bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}.vtk'.format( dens,input_suffix )
                         , grid_function = dif_F , data_type = 'element')

    #if False:
    #    status = value_assignor_starter(face_array , np.abs(dif[:,0]) , percentaje) #final_status(face_array , dif[:,0] , 2 )
    #    const_space = bempp.api.function_space(grid,  "DP", 0)
    #    Status    = bempp.api.GridFunction(const_space, fun=None, coefficients=status)
    #    Status.plot()

    if refine:

        new_face_array , new_vert_array = mesh_refiner(face_array , vert_array , np.abs(dif[:,0]) , percentaje )

        if smooth:

            aux_vert_array = np.zeros(( len(new_vert_array),3 ))

            c=0
            for vert in new_vert_array:
                aux_vert_array[c] = smothing_func(vert , r) 
                c+=1

            #print(aux_vert_array)
            
        elif not smooth:
        
            aux_vert_array = new_vert_array.copy()
            

    
    # Run over ImproveSurfMesh  --> FETK Toolkit
    # Place this program path in constants.py  !
    
    new_face_array , aux_vert_array = Improve_Mesh(new_face_array , aux_vert_array , mesh_info.path , 
                                                  mesh_info.mol_name+ '_' + str(dens) + output_suffix )
    
    vert_and_face_arrays_to_text_and_mesh( name , aux_vert_array ,
                                                new_face_array.astype(int)[:] , output_suffix, dens , Self_build=True)

    grid = Grid_loader( name , dens , output_suffix )
    
    #grid.plot()

    N_elements = len(face_array)
    
    print(N_elements)
    
    ref_end_time = time.time()
    ref_fin_time = ref_end_time - ref_init_time
    
        
    return (S_trad , S_Cooper , S_Zeb , N_elements , fin_time , ref_fin_time , operators_time ,
            matrix_assembly_time , solv_time , it_adj_count , U_it_count)

def rutine(name , N_it, percentajes , r = 1.):
    
    if name == 'sphere_cent':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,  1.E-12 ]]  )
        q = np.array( [1.] )
    
    if name == 'sphere_excent':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,   r/2. ]]  )
        q = np.array( [1.] )
    
    if name == 'sphere_dip_C5':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,   1.E-12 ] ,
                         [  1.E-12 ,  1.E-12 ,   0.5    ] ])
        q = np.array( [1. , 1. ] )
    
    if name == 'sphere_dip': 
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,    0.5 ] ,
                         [  1.E-12 ,  1.E-12 ,   -0.5    ] ])
        q = np.array( [1. , 1. ] )
        
    if name == 'sphere_dip_12':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,    0.5    ] ,
                         [  1.E-12 ,  1.E-12 ,   -0.5    ] ])
        q = np.array( [1. , 2. ] )
        
    if name == 'special':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,  0.62 ],
                             [  1.E-12 ,  0.62*np.cos(np.pi*1.5 + 5.*np.pi/180.) ,
                                                      0.62*np.sin(np.pi*1.5 + 5.*np.pi/180. ) ] ,
                             [  1.E-12 ,  0.62*np.cos(np.pi*1.5 - 5.*np.pi/180.) ,
                                                      0.62*np.sin(np.pi*1.5 - 5.*np.pi/180. )  ]
                            ])
        q = np.array( [1. , 1. , -1.]) 
    #for percentaje in percentajes:
    if True:
        
        percentaje = percentajes
        
        #for s in (True,False):
        if True:
            
            s = True
            suffixes = suffix_names(N_it)

            #x_q = np.array( [[  1.E-12 ,  1.E-12 ,   1.E-12 ]]  ) #,
                           #  [  1.E-12 ,  1.E-12 ,  -0.5    ] ])
            #q = np.array( [1.] ) #, 2. ]) #, 1. ] )

            #r = 1.
            #name = 'sphere_cent'
            dens = 0

            path = os.path.join( 'Molecule' , name , str( int(percentaje*100)) , str(s) )            
            os.system('mkdir -p {0}'.format( path ) )

            bempp.api.set_ipython_notebook_viewer()
            bempp.api.PLOT_BACKEND = "ipython_notebook"

            # Creates Initial Decahedron
            dec_v_array , dec_f_array = decahedron_points(r)
            vert_and_face_arrays_to_text_and_mesh( name, dec_v_array , dec_f_array+1 , '-s0' 
                                                      , dens=0 , Self_build=True)
            # Refines to obtain a better aproach of a sphere
            main(name , dens , '-s0' , '-s1' , percentaje=2.5 , r=r , x_q=x_q , q=q , smooth = True)

            Resultados = open('Molecule/{0}/Resultados_{0}.txt'.format( name ) , 'w+')

            Resultados.write('x_q = \n')
            Resultados.write(str(x_q) )
            Resultados.write('\n q =   \n')
            Resultados.write(str(  q) )
            Resultados.write('  \n')


            Resultados.write('Percentaje & Number of elements & Strad & SAprox & SEx & smooth & time \n')
            Resultados.close()

            c=1
            for i in suffixes[1:-1]:

                (S_trad , S_Cooper , S_Zeb , N_elements , fin_time , ref_time ,
                 operators_time , matrix_assembly_time , solv_time ,
                 it_adj_count , U_it_count) = main(name , dens , suffixes[c]  , suffixes[c+1] , percentaje , r , x_q , q , smooth = s )
                
                
                print(c , suffixes[c] , suffixes[c+1] )
                
                Resultados = open('Molecule/{0}/Resultados_{0}.txt'.format( name ) , 'a')
                
                text = '- & '+str(i) + ' & {0:d} & {1:.10f} & {2:.10f} & {3:.10f} & {4:.10f} & {5:.10f}'.format( 
                    N_elements , S_trad , S_Cooper , S_Zeb , fin_time, ref_time  
                ) +' & {0:.10f} & {1:.10f} & {2:.10f} & {3:d} & {4:d} \n'.format(  
                        operators_time , matrix_assembly_time , solv_time , it_adj_count , U_it_count  )
                Resultados.write(text)
                #Resultados.write('{4:.2f} & {3:d} & {0:.10f} & {1:.10f} & {2:.10f} & {5} & {6:.10f}\n'.format(
                #            line[0] , line[1], line[2] , line[3] , percentaje, str(s) , line[4]))
                Resultados.close()

                #res = np.vstack( (res , np.array( (name , line[0] , line[1], line[2]
                #                                  , line[3] , percentaje, str(s) , line[4])) ) )
                c+=1
                
                #os.system( "python3.6 File_converter_python3.py" )

                if c==5 and percentaje == 1.5: break

            #R   = 1.

            a   = 1.
            N   = 20

            f_ex = solution(q, x_q, ep_m, ep_s, r, k , a, N)
            print('Exact Solution: {0:5f}--------------------------------------------'.format(f_ex))
            #print(solution_2(q, x_q, ep_m, ep_s, r, kappa, a, N))
            Resultados = open('Molecule/{0}/Resultados_{0}.txt'.format( name ) , 'a')
            Resultados.write('Exact solution: {0:10f}'.format(f_ex))
            
            Resultados.write('Conditions\n')
            
            Resultados.write('Smooth and Use_Gamer: Enabled \n')
            Resultados.write('Adjoint mesh with uniform refinement: Enabled \n')
            Resultados.write('Note: None \n')
            
            Resultados.close()

            moving_files(name , percentaje , s)
            
        
            
    return None

N_it = 31
percentajes = 0.1

name = 'special'
rutine(name , N_it, percentajes , r = 1.)#
rutine(name , 5   , 1.5         , r = 1.)

#name = 'sphere_excent'
#rutine(name , N_it, percentajes , r = 1.)
#rutine(name , 5   , 1.5         , r = 1.)
