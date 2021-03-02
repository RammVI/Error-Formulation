
# checked 19.09
# 11-09 added adjoint mesh
from pathos.multiprocessing import ProcessingPool as Pool
import bempp.api
import numpy as np
import os
import time

from bempp.api.operators.potential import laplace as lp
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

from constants import mesh_info
from constants import values
from constants import *

from quadrature import *

from Mesh_Ref_V2 import *
from Grid_Maker_R2 import *

# for debuging
from random import randint
from constants import potential

# --------------------------------------------------------------------------------

def zero_i(x, n, domain_index, result):
    result[:] = 0

def u_s_G(x,n,domain_index,result):
    global ep_m 
    result[:] =  1. / (4.*np.pi*ep_m)  * np.sum( mesh_info.q / np.linalg.norm( x - mesh_info.x_q , axis = 1 )  )


def du_s_G(x,n,domain_index,result):
    global ep_m
    result[:] = -1./(4.*np.pi*ep_m)  * np.sum( np.dot( x-
                            mesh_info.x_q , n  )  * mesh_info.q / np.linalg.norm( x - mesh_info.x_q , axis = 1  )**3  )

def harmonic_component(dirichl_space , neumann_space , dual_to_dir_s , u_s , du_s):
    
    global ep_m , ep_s, k

    # Se crean los operadores asociados al sistema, que dependen de los espacios de las soluciones
    # identity : I  : matriz identidad
    # dlp_in : K_in : Double-Layer operator para la region interior
    # slp_in : V_in : Single-Layer operator para la region interior
    # _out : Mismos operadores pero para la region exterior, con k=kappa=0.125
    identity = sparse.identity(     dirichl_space, dirichl_space, dual_to_dir_s)
    slp_in   = laplace.single_layer(neumann_space, dirichl_space, dual_to_dir_s)
    dlp_in   = laplace.double_layer(dirichl_space, dirichl_space, dual_to_dir_s)

    # V_in du_s = (1/2+K_in)u_h = -(1/2+K_in)u_s (BC)
    sol, info,it_count = bempp.api.linalg.gmres( slp_in, -(dlp_in+0.5*identity)*u_s , return_iteration_count=True, tol=1e-8)
    print("The linear system for du_h was solved in {0} iterations".format(it_count))



    u_h = -u_s
    du_h = sol
    
    return u_h , du_h

def regular_component(dirichl_space , neumann_space , dual_to_dir_s , du_s , du_h):
    
    global ep_m , ep_s, k
    
    identity = sparse.identity(     dirichl_space, dirichl_space, dual_to_dir_s)
    slp_in   = laplace.single_layer(neumann_space, dirichl_space, dual_to_dir_s)
    dlp_in   = laplace.double_layer(dirichl_space, dirichl_space, dual_to_dir_s)
    slp_out  = modified_helmholtz.single_layer(neumann_space, dirichl_space, dual_to_dir_s, k)
    dlp_out  = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dual_to_dir_s, k)
    
    # Se crea la matriz / Lado izquierdo de la ecuacion
    # | ( I/2 + K_L-in  )     (      -V_L-in     ) |  u_r  = 0
    # | ( I/2 - K_Y-out )     ( ep_m/ep_s V_Y-out) | du_r  = ep_m/ep_s V_Y-out*(du_s+du_h)  (BC)
    blocked = bempp.api.BlockedOperator(2, 2)
    blocked[0, 0] = 0.5*identity + dlp_in
    blocked[0, 1] = -slp_in
    blocked[1, 0] = 0.5*identity - dlp_out
    blocked[1, 1] = (ep_m/ep_s)*slp_out

    # Se crea el lado derecho de la ecuacion 
    zero = bempp.api.GridFunction(dirichl_space, fun=zero_i)
    rhs = [ zero ,  -slp_out *(ep_m/ep_s)* (du_s+du_h)]

    # Y Finalmente se resuelve para u_r y du_r
    sol, info,it_count = bempp.api.linalg.gmres( blocked , rhs, return_iteration_count=True, tol=1e-8)
    print("The linear system for u_r and du_r was solved in {0} iterations".format(it_count))
    u_r , du_r = sol
    
    return u_r , du_r

# --------------------------------------------------------------------------------

def carga_i(x, n, domain_index, result):
    global ep_m

    # Right side of the eqn, with the Green function convolution
    result[:] = np.sum(mesh_info.q/np.linalg.norm( x - mesh_info.x_q, axis=1 ))/(4.*np.pi*ep_m)

def U_tot_boundary(grid , return_time =False , save_plot=False , tolerance = 1e-8):
    '''
    Returns the total electrostatic potential in the boundary.
    Parameters:
    grid   : Bempp object
    '''
    init_time_spaces = time.time()
    potential.dirichl_space_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order)
    potential.neumann_space_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order) 
    potential.dual_to_dir_s_u = bempp.api.function_space(grid,  mesh_info.u_space, mesh_info.u_order)
    spaces_time      = time.time()-init_time_spaces
    U, dU , operators_time , assembly_time , GMRES_time , UdU_time , it_count= U_tot(potential.dirichl_space_u ,
                        potential.neumann_space_u , potential.dual_to_dir_s_u , tolerance )
    if save_plot:
        
        U_func  = bempp.api.GridFunction(potential.dirichl_space_u , fun=None, coefficients= U.coefficients.real)
        dU_func = bempp.api.GridFunction(potential.neumann_space_u , fun=None, coefficients=dU.coefficients.real)
        
        bempp.api.export('Molecule/' + mesh_info.mol_name +'/' + mesh_info.mol_name + '_{0}{1}_U.vtk'.format(
                                      mesh_info.mesh_density,   mesh_info.suffix )
                     , grid_function = U_func , data_type = 'element')
        
        bempp.api.export('Molecule/' + mesh_info.mol_name +'/' + mesh_info.mol_name + '_{0}{1}_dU.vtk'.format(
                                      mesh_info.mesh_density,   mesh_info.suffix )
                     , grid_function = dU_func , data_type = 'element')
    
    if return_time:
        return U, dU , spaces_time , operators_time , assembly_time , GMRES_time , UdU_time , it_count
    
    return U , dU


def U_tot(dirichl_space , neumann_space , dual_to_dir_s , tolerance):
    '''
    Computes the Total electrostatic mean potential on the boundary.
    Params:
    dirichl_space
    neumann_space
    dual_to_dir_s
    tolerance     : Tolerance of the solver
    '''
    
    init_time_operators = time.time()

    identity = sparse.identity(     dirichl_space, dirichl_space, dual_to_dir_s)
    slp_in   = laplace.single_layer(neumann_space, dirichl_space, dual_to_dir_s)
    dlp_in   = laplace.double_layer(dirichl_space, dirichl_space, dual_to_dir_s)

    slp_out  = modified_helmholtz.single_layer(neumann_space, dirichl_space, dual_to_dir_s, k)
    dlp_out  = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dual_to_dir_s, k)
    
    charged_grid_fun = bempp.api.GridFunction(dirichl_space, fun=q_times_G_L)
    zero_grid_fun    = bempp.api.GridFunction(neumann_space, fun=zero_i     )
    
    operators_time = time.time() - init_time_operators

    init_time_matrix = time.time()
    blocked = bempp.api.BlockedOperator(2, 2)
    blocked[0, 0] = 0.5*identity + dlp_in
    blocked[0, 1] = -slp_in
    blocked[1, 0] = 0.5*identity - dlp_out
    blocked[1, 1] = (ep_m/ep_s)*slp_out

    rhs = [charged_grid_fun, zero_grid_fun]
    
    assembly_time = time.time() - init_time_matrix
    
    print('GMRES Tolerance = {0}'.format(str(tolerance)))
    init_time_GMRES = time.time()
    sol, info,it_count = bempp.api.linalg.gmres( blocked, rhs , return_iteration_count=True , tol=tolerance ,
                                               use_strong_form=True)
    GMRES_time = time.time()-init_time_GMRES
    print("The linear system for U_tot was solved in {0} iterations".format(it_count))
    init_time_UdU = time.time()
    U , dU = sol
    UdU_time      = time.time()-init_time_UdU
    
           
    return U, dU , operators_time , assembly_time , GMRES_time , UdU_time , it_count

def U_Reac(U, dU , dirichl_space , neumann_space ):
    
    U_s  = bempp.api.GridFunction(dirichl_space , fun =  u_s_G)
    dU_s = bempp.api.GridFunction(neumann_space , fun = du_s_G)
    
    U_R  =  U -  U_s
    dU_R = dU - dU_s
    
    return U_R , dU_R

def S_trad_calc_R( dirichl_space, neumann_space , U , dU ):

    # Se definen los operadores
    slp_in_O = lp.single_layer(neumann_space, mesh_info.x_q.transpose()) 
    dlp_in_O = lp.double_layer(dirichl_space, mesh_info.x_q.transpose())

    # Y con la solucion de las fronteras se fabrica el potencial evaluada en la posicion de cada carga
    U_R_O = slp_in_O * dU  -  dlp_in_O * U

    # Donde agregando algunas constantes podemos calcular la energia de solvatacion S
    
    S     = K * np.sum(mesh_info.q * U_R_O).real
    print("Solvation Energy : {:7.8f} [kcal/mol] ".format(S) )
    
    return S


# --------------------------------------------------------------------------------




def phi_in_Adjoint_Mesh(mol_name , face_array , vert_array , dens , input_suffix , return_grid = False ):
    '''
    Finds and create the adjoint mesh.
    mol_name : Molecule/Ion name, only for files saving.
    face_array : array of faces 
    vert_array : array of vertices
    dens       : used mesh density
    input_suffix : Normaly related to a number of iterations. If doing for a mesh obtained via MSMS/NanoShaper
                   use "-0".
    return_grid : Boolean 
    '''
    
    adj_face , adj_vertex = mesh_refiner(face_array , vert_array , np.ones((len(face_array[0:,]))) , 1.5 )

    vert_and_face_arrays_to_text_and_mesh( mol_name , adj_vertex , adj_face.astype(int) , input_suffix + '_adj' ,
                                          dens=dens, Self_build=True)
    
    adj_grid = Grid_loader( mol_name , dens , input_suffix + '_adj' )

    adj_face_array = np.transpose(adj_grid.leaf_view.elements) + 1
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)

    adj_el_pos = elements_position_in_normal_grid(adj_face_array , adj_vert_array , face_array , vert_array )
    
    potential.dirichl_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space ,
                                                           mesh_info.phi_order)
    potential.neumann_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space ,
                                                           mesh_info.phi_order) 
    potential.dual_to_dir_s_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space ,
                                                           mesh_info.phi_order)
    
    phi , dphi , it_count = adjoint_equation( potential.dirichl_space_phi , 
                        potential.neumann_space_phi , potential.dual_to_dir_s_phi)
    
    if return_grid:
        return phi , dphi , adj_grid
    
    
    return phi , dphi


def S_Exact_in_Adjoint_Mesh_with_N_Ref(mol_name , grid  , dens , input_suffix , N , N_ref ,
                                        save_energy_plot=False  , test_mode = False , return_times = False):
    
    face_array = np.transpose(grid.leaf_view.elements)
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    aux_face = face_array.copy()
    aux_vert = vert_array.copy()
    
    if N_ref == 0:
        adj_grid = grid
        
    elif N_ref>=1:
        i=1
        while i <= N_ref:

            new_face , new_vertex = mesh_refiner(  aux_face +1 , aux_vert , np.ones((len(aux_face[0:,]))) , 1.5 )

            vert_and_face_arrays_to_text_and_mesh( mol_name , new_vertex , new_face.astype(int), input_suffix +
                                                  '_adj_'+ str(i), dens=dens, Self_build=True)

            aux_face , aux_vert = new_face.copy()- 1 , new_vertex.copy()
            i+=1

        adj_grid = Grid_loader( mol_name , dens , input_suffix + '_adj_' + str(N_ref) )       
        
    
    adj_face_array = np.transpose(adj_grid.leaf_view.elements) + 1
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)

    adj_el_pos = check_contained_triangles__(grid , adj_grid)
    
    #print(len(adj_face_array))
    #print(adj_el_pos)
    init_time_spaces_adj  = time.time()
    dirichl_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order)
    neumann_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order) 
    dual_to_dir_s_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order)
    spaces_time_adj       = time.time()-init_time_spaces_adj
    
    phi , dphi , operators_time_adj , matrix_time_adj , GMRES_time_adj , phidphi_time , it_count = adjoint_equation( dirichl_space_phi , neumann_space_phi , dual_to_dir_s_phi , 
                                            save_plot = True , suffix = '_'+str(N_ref) )
    init_time_S_Ex = time.time()
    S_Ex    , rearange_S_Ex_i  , S_Ex_i , _= Exact_aproach_with_u_s_Teo( adj_face_array , adj_vert_array , phi , dphi , N , 
                                                 grid_relation = adj_el_pos , return_values_on_Adj_mesh = True)
    S_Ex_time      = time.time()-init_time_S_Ex
    N_el_adjoint = len(adj_face_array)
    
    
    if test_mode:
        
        const_space = bempp.api.function_space(adj_grid,  "DP", 0)
        counter     = bempp.api.GridFunction(const_space, fun=None, coefficients=adj_el_pos )
        bempp.api.export('Molecule/' + mol_name +'/' + mol_name + '_{0}{1}_Counter_{2}.msh'.format( 
                                        dens, input_suffix , N_ref )
                     , grid_function = counter , data_type = 'element')
    
    if save_energy_plot:
        const_space = bempp.api.function_space(adj_grid,  "DP", 0)
        S_Ex_BEMPP  = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ex_i[:,0])
        bempp.api.export('Molecule/' + mol_name +'/' + mol_name + '_{0}{1}_S_Exact_{2}.msh'.format( 
                                        dens, input_suffix , N_ref )
                     , grid_function = S_Ex_BEMPP , data_type = 'element')
    
    return (S_Ex , rearange_S_Ex_i , it_count , N_el_adjoint , spaces_time_adj , operators_time_adj , 
                matrix_time_adj , GMRES_time_adj , phidphi_time , S_Ex_time )


    
def S_Exact_in_Adjoint_Mesh(mol_name , face_array , vert_array , dens , input_suffix , N):
    
    adj_face , adj_vertex = mesh_refiner(face_array , vert_array , np.ones((len(face_array[0:,]))) , 1.5 )

    vert_and_face_arrays_to_text_and_mesh( mol_name , adj_vertex , adj_face.astype(int) , input_suffix + '_adj' ,
                                          dens=dens, Self_build=True)
    
    adj_grid = Grid_loader( mol_name , dens , input_suffix + '_adj' )
    
    adj_face_array = np.transpose(adj_grid.leaf_view.elements) + 1
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)

    adj_el_pos = elements_position_in_normal_grid(adj_face_array , adj_vert_array , face_array , vert_array )
    
    dirichl_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order)
    neumann_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order) 
    dual_to_dir_s_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order)
    
    phi , dphi , it_count = adjoint_equation( dirichl_space_phi , neumann_space_phi , dual_to_dir_s_phi)
    
    S_Exact    , S_Exact_i = Exact_aproach_with_u_s_Teo( adj_face_array , adj_vert_array , phi , dphi , N ,
                                                 grid_relation=adj_el_pos)
    
    rearange_S_Exact_i = np.zeros((len(face_array),1))
    
    c=0
    for G_i in S_Exact_i:
        rearange_S_Exact_i[int(adj_el_pos[c])] += G_i
        c+=1
        
    #print(rearange_S_Exact_i)
    
    return S_Exact , rearange_S_Exact_i , it_count


# ----------------------------------------------------------------
    

def q_times_G_L(x, n, domain_index, result):
    global ep_m
    result[:] = 1. / (4.*np.pi*ep_m)  * np.sum( mesh_info.q  / np.linalg.norm( x - mesh_info.x_q, axis=1 ) )

def adjoint_equation( dirichl_space , neumann_space , dual_to_dir_s , save_plot = True ,suffix = ''):
    
    global ep_m , ep_s , k
    
    init_time_operators = time.time()
    identity = sparse.identity(     dirichl_space, dirichl_space, dual_to_dir_s)
    slp_in   = laplace.single_layer(neumann_space, dirichl_space, dual_to_dir_s)
    dlp_in   = laplace.double_layer(dirichl_space, dirichl_space, dual_to_dir_s)
    slp_out  = modified_helmholtz.single_layer(neumann_space, dirichl_space, dual_to_dir_s, k)
    dlp_out  = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dual_to_dir_s, k)
    operators_time_adj = time.time() - init_time_operators
    
    init_time_operators = time.time()
    blocked = bempp.api.BlockedOperator(2, 2)
    blocked[0, 0] = 0.5*identity + dlp_in
    blocked[0, 1] = -slp_in
    blocked[1, 0] = 0.5*identity - dlp_out
    blocked[1, 1] = (ep_m/ep_s)*slp_out
    
    
    zero = bempp.api.GridFunction(dirichl_space , fun=zero_i)
    P_GL = bempp.api.GridFunction(dirichl_space, fun=q_times_G_L)
    rs_r = [P_GL , zero]
    matrix_time_adj = time.time() - init_time_operators
    
    init_time_GMRES_adj = time.time()
    sol_r, info,it_count = bempp.api.linalg.gmres( blocked, rs_r , return_iteration_count=True, tol=1e-5,
                                                  use_strong_form=True )
    GMRES_time_adj = time.time()-init_time_GMRES_adj
    print("The linear system for phi was solved in {0} iterations".format(it_count))
    phidphi_time = time.time()
    phi_r , dphi_r = sol_r
    phidphi_time = time.time()-phidphi_time 
    
    if save_plot:
        
        phi_func  = bempp.api.GridFunction( dirichl_space , fun=None, coefficients= phi_r.coefficients.real)
        dphi_func = bempp.api.GridFunction( neumann_space , fun=None, coefficients=dphi_r.coefficients.real)
        
        bempp.api.export('Molecule/' + mesh_info.mol_name +'/' + mesh_info.mol_name + '_{0}{1}_phi_{2}.msh'.format(
                                      mesh_info.mesh_density,   mesh_info.suffix , suffix )
                     , grid_function = phi_func )
        
        bempp.api.export('Molecule/' + mesh_info.mol_name +'/' + mesh_info.mol_name + '_{0}{1}_dphi_{2}.vtk'.format(
                                      mesh_info.mesh_density,   mesh_info.suffix , suffix )
                     , grid_function = dphi_func , data_type = 'node')
    
    return phi_r , dphi_r  , operators_time_adj , matrix_time_adj , GMRES_time_adj , phidphi_time , it_count

def S_trad_calc( dirichl_space, neumann_space , u_h , du_h , u_r , du_r):
    
    # En base a los puntos donde se encuentran las cargas, calculemos el potencial u_r y u_h
    # Esto luego de que podemos escribir la energia de solvatacion como
    # G_solv = Sum_i q_i *u_reac = Sum_i q_i * (u_h+u_r)           evaluado en cada carga.

    # Se definen los operadores
    slp_in_O = lp.single_layer(neumann_space, mesh_info.x_q.transpose()) 
    dlp_in_O = lp.double_layer(dirichl_space, mesh_info.x_q.transpose())

    # Y con la solucion de las fronteras se fabrica el potencial evaluada en la posicion de cada carga
    u_r_O = slp_in_O * du_r  -  dlp_in_O * u_r
    u_h_O = slp_in_O * du_h  -  dlp_in_O * u_h

    terms =  u_r_O + u_h_O

    # Donde agregando algunas constantes podemos calcular la energia de solvatacion S
    
    S     = K * np.sum(mesh_info.q * terms).real
    print("Solvation Energy : {:7.8f} [kcal/mol] ".format(S) )
    
    return S

def S_Aproximate_calc( face_array , vert_array , phi_r , dphi_r , U_Reac , dU_Reac , N):
    
    Solv_Aproximate = np.zeros((len(face_array),1))
    c = 0
    for face in face_array:

        I1 = int_calc_i( face , face_array , vert_array , phi_r , mesh_info.phi_space 
                        , mesh_info.phi_order , ep_m*dU_Reac , mesh_info.u_space , mesh_info.u_order  , N)
        
        I2 = int_calc_i( face , face_array , vert_array , ep_m * dphi_r , mesh_info.phi_space
                        , mesh_info.phi_order , U_Reac , mesh_info.u_space , mesh_info.u_order  , N)
        Solv_Aproximate[c] = I1-I2

        c+=1

    S_Aproximate_i = K*Solv_Aproximate
    S_Aproximate = np.sum(S_Aproximate_i )
    print('Aproximate Solv = {0:10f} '.format(S_Aproximate)) 
    
    return S_Aproximate , S_Aproximate_i


def Aproximated_Sol_Adj_UDP0(U_R , dU_R , phi , dphi , face_array , vert_array , 
                             face_array_adj , vert_array_adj , N , coarse_grid , adj_grid , N_ref ,
                            return_relation=False):
    '''
    Returns the integral over Gamma for phi and dphi per element.
    U MUST BE IN DP0 IN ORDER TO WORK.
    Params: U_R , dU_R , phi , dphi , face_array , vert_array , face_array_adj , vert_array_adj , N 
    face_array , vert_array         : Array of faces and vertices from the coarse mesh
    face_array_adj , vert_array_adj : Array of faces and vertices from the adjoint mesh
    N: number of quadrature points used (Only for phi), as U is in DP0 means is constant per element.
    '''
    
    phi_array , dphi_array = integral_per_element(phi , dphi , face_array_adj , vert_array_adj ,
                                                  mesh_info.phi_space , mesh_info.phi_order ,
                                                  mesh_info.phi_space , mesh_info.phi_order , N , adj_grid)
    
    
    relationship = check_contained_triangles__( coarse_grid , adj_grid )
        
    rearange_S_Aprox_i = np.zeros((len(face_array),1))
    
    c_adj=0
    for c in relationship:
        rearange_S_Aprox_i[c] += ep_m *(  dU_R.coefficients.real[c]*phi_array[c_adj] - 
                                           U_R.coefficients.real[c]*dphi_array[c_adj]   )
        c_adj+=1
        
    rearange_S_Aprox_i = K*rearange_S_Aprox_i.copy()

    S_Aprox = np.sum(rearange_S_Aprox_i )#[0]
    
    print('Aproximated Solvation Energy : {0:10f}'.format(S_Aprox))
        
    if return_relation==True:
        return S_Aprox , rearange_S_Aprox_i , relationship
    
    
    
    return S_Aprox , rearange_S_Aprox_i
       
       
def S_Exact_calc( face_array , vert_array , phi , dphi , u_s , du_s , N):
    
    Solv_Exact = np.zeros((len(face_array),1))
    c = 0
    print(mesh_info.u_s_order)
    for face in face_array:

        I1 = int_calc_i( face , face_array , vert_array , phi , mesh_info.phi_space 
                        , mesh_info.phi_order , ep_m*(du_s) , mesh_info.u_s_space , mesh_info.u_s_order  , N)
        
        I2 = int_calc_i( face , face_array , vert_array ,  dphi , mesh_info.phi_space
                         , mesh_info.phi_order , ep_m*(u_s) , mesh_info.u_s_space , mesh_info.u_s_order , N)
        Solv_Exact[c] = I2-I1

        c+=1
    Solv_Exact_i = Solv_Exact
    S_Exact = K*np.sum(Solv_Exact )
    print('Exact Solv = {0:10f} '.format(S_Exact)) 
    

    
    return S_Exact , Solv_Exact_i
       

def delta_G_exact(grid , U , dU , U_space , U_order , N , save_energy_plot = False):
    '''
    Calculates the exact solvation energy using theorical u_c and du_c
    Params:
    face_array : Array of faces
    vert_array : Array of vertices
    U          : BEMPP object
    dU         : BEMPP object
    N          : Number of points for the gauss cuadrature
    save_energy_plots = True or False : If energy plots wants to be saved in the mesh directory
    Returns 
    S**{ex} , np.array(S**{ex} per element) 
    '''
    
    face_array = np.transpose(grid.leaf_view.elements) + 1
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    normals = normals_to_element( face_array , vert_array )
    
    Solv_Ex = np.zeros((len(face_array),1))
    c = 0
    
    for face in face_array:
        
        f1 , f2 , f3 = face-1
        v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]
        
        normal = normals[c]
        
        Area = 0.5 * np.linalg.norm( np.cross(v2-v1 , v3-v1) )
        
        if U_order != 0:
            A = matrix_lineal_transform( v1 , v2 , v3 )
        elif U_order == 0:
            A = 1.
        
        X_K , W = evaluation_points_and_weights(v1,v2,v3 , N)
        
        U_aux  = unpack_info( face , face_array, vert_array , U  , U_space , U_order )
        dU_aux = unpack_info( face , face_array, vert_array , dU , U_space , U_order )
        
        I1 , I2 = 0. , 0. 
        
        point_count = 0 
        
        for x in X_K:
            
            U_local  = local_f( x , A , U_aux  , U_order)
            dU_local = local_f( x , A , dU_aux , U_order)
            
            u_s_local  = u_s_Teo( x )
            du_s_local = du_s_Teo( x , normal )
            
            I1 += ep_m * dU_local * u_s_local * W[point_count] 
            
            I2 += ep_m * U_local * du_s_local * W[point_count]
            
            point_count+=1
            
        Solv_Ex[c] = (I1-I2)*Area

        c+=1
        
    Solv_Ex_i = K * Solv_Ex
    
    Solv_Ex_total = np.sum(Solv_Ex_i )
    print('Estimated Solvation Energy : {0:10f} '.format(Solv_Ex_total)) 
    
    if save_energy_plot:
        const_space = bempp.api.function_space(grid,  "DP", 0)
        S_Ap_BEMPP  = bempp.api.GridFunction(const_space, fun=None, coefficients=Solv_Ex_i[:,0])
        bempp.api.export('Molecule/' + mesh_info.mol_name +'/' + mesh_info.mol_name +
                        '_{0}{1}_S_Aprox.vtk'.format( mesh_info.mesh_density, mesh_info.suffix)
                     , grid_function = S_Ap_BEMPP , data_type = 'element')
        
    
    return Solv_Ex_total , Solv_Ex_i

def delta_G_exact_new(grid , U , dU , U_space , U_order , N , save_energy_plot = False):
    '''
    Calculates the exact solvation energy using theorical u_c and du_c
    Params:
    face_array : Array of faces
    vert_array : Array of vertices
    U          : BEMPP object
    dU         : BEMPP object
    N          : Number of points for the gauss cuadrature
    save_energy_plots = True or False : If energy plots wants to be saved in the mesh directory
    Returns 
    S**{ex} , np.array(S**{ex} per element) 
    '''
    
    face_array = np.transpose(grid.leaf_view.elements) 
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    normals = normals_to_element( face_array , vert_array )
    
    Solv_Ex = np.zeros((len(face_array),1))
    c = 0
    
    if U_order == 0 and U_space == 'DP':
        
        for face in face_array:

            f1 , f2 , f3 = face
            v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]

            normal = normals[c]

            Area = 0.5 * np.linalg.norm( np.cross(v2-v1 , v3-v1) )

            X_K , W = evaluation_points_and_weights(v1,v2,v3 , N)

            u_c_func  = lambda x_i : u_s_Teo(x_i)
            du_c_func = lambda x_i : du_s_Teo(x_i, normal)
            u_c_face  = np.array([ u_c_func(x_i) for x_i in X_K]) 
            du_c_face = np.array([du_c_func(x_i) for x_i in X_K]) 

            U_local , dU_local = U.coefficients.real[c] , dU.coefficients.real[c]

            I2 =  np.sum( U_local  * du_c_face * W )
            I1 =  np.sum( dU_local * u_c_face  * W )

            Solv_Ex[c] = (I1-I2)*Area         

            c+=1

        Solv_Ex_i = K * ep_m * Solv_Ex

        Solv_Ex_total = np.sum(Solv_Ex_i )
        print('Estimated Solvation Energy : {0:10f} '.format(Solv_Ex_total)) 

        if save_energy_plot:
            const_space = bempp.api.function_space(grid,  "DP", 0)
            S_Ap_BEMPP  = bempp.api.GridFunction(const_space, fun=None, coefficients=Solv_Ex_i[:,0])
            bempp.api.export('Molecule/' + mesh_info.mol_name +'/' + mesh_info.mol_name +
                            '_{0}{1}_S_Aprox.vtk'.format( mesh_info.mesh_density, mesh_info.suffix)
                         , grid_function = S_Ap_BEMPP , data_type = 'element')


    return Solv_Ex_total , Solv_Ex_i

def delta_G_tent( grid , U , dU , U_space , U_order , N ):
    
    face_array = np.transpose(grid.leaf_view.elements) 
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    if U_space == 'DP' and U_order == 0:
    
        Solv_Ex = np.zeros([len(face_array) , 1])


        mesh = trimesh.Trimesh(vertices=vert_array , faces=face_array)
        normals = mesh.face_normals
        Areas = mesh.area_faces

        quadrule = quadratureRule_fine(N)
        X_K , W  = quadrule[0].reshape(-1,3) , quadrule[1]

        c=0
        for face in face_array:
                        
            v1 , v2 , v3 = vert_array[face[0]] , vert_array[face[1]] , vert_array[face[2]]
            
            X = evaluation_points_and_weights_new(v1,v2,v3 , N , X_K , W)
            
            u_c_func  = lambda x_i : u_s_Teo(x_i )
            du_c_func = lambda x_i : du_s_Teo(x_i, normals[c] )
            u_c_face  = np.array([ u_c_func(x_i) for x_i in X]) 
            du_c_face = np.array([du_c_func(x_i) for x_i in X]) 


            U_local , dU_local = U.coefficients.real[c] , dU.coefficients.real[c]

            I2 =  np.sum( U_local  * du_c_face * W )
            I1 =  np.sum( dU_local * u_c_face  * W )

            Solv_Ex[c] = (I1-I2)*Areas[c]       

            c+=1
    
    Solv_Ex_i = K * ep_m * Solv_Ex
    
    print('Aproximate solvation energy: {0:.6f}[kcal/kmol]'.format(np.sum(Solv_Ex_i)))
    return np.sum(Solv_Ex_i) , Solv_Ex_i


def Exact_aproach_with_u_s_Teo( face_array , vert_array , phi , dphi , N , grid_relation=None ,
                            return_values_on_Adj_mesh = False ):
    
    if np.min(face_array)==0:
        face_array += 1
    
    normals = normals_to_element( face_array , vert_array )
    
    Solv_Exact = np.zeros((len(face_array),1))
    c = 0
    
    for face in face_array:
        
        f1 , f2 , f3 = face-1
        v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]
        
        normal = normals[c]
        
        
        A = matrix_lineal_transform( v1 , v2 , v3 )
        
        Area = 0.5 * np.linalg.norm( np.cross(v2-v1 , v3-v1) )
        
        X_K , W = evaluation_points_and_weights(v1,v2,v3 , N)
        
        phi_a  = unpack_info( face , face_array, vert_array , phi  , mesh_info.phi_space , mesh_info.phi_order)
        dphi_a = unpack_info( face , face_array, vert_array , dphi , mesh_info.phi_space , mesh_info.phi_order)
        
        I1 , I2 = 0. , 0. 
        
        point_count = 0        
        for x in X_K:
            
            phi_local  = local_f( x , A , phi_a  , mesh_info.phi_order)
            dphi_local = local_f( x , A , dphi_a , mesh_info.phi_order)
            
            u_s_local  = u_s_Teo( x )
            du_s_local = du_s_Teo( x , normal )
            
            I1 += ep_m * dphi_local * u_s_local * W[point_count] 
            
            I2 += ep_m * phi_local * du_s_local * W[point_count]
            
            point_count+=1
            
        Solv_Exact[c] = (I1-I2)*Area

        c+=1
        
    Solv_Exact_i = K * Solv_Exact
    
    S_Exact = np.sum(Solv_Exact_i )
    print('Estimated Exact Solvation Energy : {0:10f} '.format(S_Exact)) 
    
    if np.any(type(grid_relation)!=None):
        c=0
        
        rearange_S_Ex_i = np.zeros(( np.max(grid_relation.astype(int))+1,1))
        
        for G_i in Solv_Exact_i:
            rearange_S_Ex_i[int(grid_relation[c])] += G_i
            c+=1
            
        if return_values_on_Adj_mesh:
            return S_Exact , rearange_S_Ex_i , Solv_Exact_i , grid_relation
        
    elif type(grid_relation)==None:
        print('Cant relate elements from the adj with the original mesh. Continuee....')
        return S_Exact
        
    
    return S_Exact , rearange_S_Ex_i

def u_s_Teo( x ):
    
    return (1. / (4.*np.pi*ep_m) ) * np.sum( mesh_info.q / np.linalg.norm( x - mesh_info.x_q, axis=1 ) )

def du_s_Teo(x,n):
    
    return -1./(4.*np.pi*ep_m)  * np.sum( np.dot( x-
                            mesh_info.x_q , n)  * mesh_info.q / np.linalg.norm( x - mesh_info.x_q, axis=1 )**3 )

def normals_to_element( face_array , vert_array ):
    '''
    Calculates the normals pointing outwards the inner domain
    Inputs
    face_array  : Array containing vertices position in vert_array for each element.
    vert_array  : Array of vertices [N,3]
    '''
    normals = np.empty((0,3))
    element_cent = np.empty((0,3))
    
    for face in face_array:
        
        f1,f2,f3 = face-1
        v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]
        n = np.cross( v2-v1 , v3-v1 ) 
        normals = np.vstack((normals , n/np.linalg.norm(n) )) 
        element_cent = np.vstack((element_cent, (v1+v2+v3)/3. ))

    return normals

def z_value(x, n, domain_index, result):
    result[:] = x[0]
    

    
def local_U_Reac_interior( x ):
    '''
    Computes the local reaction potential value for a given position and the boundary solution
    U_boundary: TOTAL potential on the boundary
    dirichl_space : Space of U
    neumann_space : Space of DU/Dn
    x         : Position or array of points
    '''
    
    aux_x = np.array([x]).transpose()
    
    slp_in_O = lp.single_layer(potential.neumann_space_u, aux_x ) 
    dlp_in_O = lp.double_layer(potential.dirichl_space_u, aux_x )
    
    U_reac   =  slp_in_O * potential.dU.real  -  dlp_in_O * potential.U.real
    
    return U_reac

def local_U_interior( x ):
    '''
    Calculates the total electrostatic potential in the molecular region.
    U_reac  : Reaction potential at x.
    x       : Position
    '''
    U_S = u_s_Teo( x )
    
    return local_U_Reac_interior( x ) + U_S


def phi_with_N_ref(mol_name , coarse_grid , face_array , vert_array , dens ,
                     input_suffix , N_ref , return_grid = False , calculate_phi = True ):
    '''
    Finds and creates the adjoint mesh using N_ref cycles of UNIFORM refinement.
    mol_name : Molecule/Ion name, only for files saving.
    face_array : array of faces 
    vert_array : array of vertices
    dens       : used mesh density
    input_suffix : Normaly related to a number of iterations. If doing for a mesh obtained via MSMS/NanoShaper
                   use "-0".
    return_grid : Boolean
    '''
    
    aux_face = face_array.copy()
    aux_vert = vert_array.copy().astype(float)
    
    
    if N_ref == 0:
        adj_grid = coarse_grid
        
    elif N_ref>=1:
    
        for i in range(1,N_ref+1):

            new_face , new_vertex = mesh_refiner(aux_face +1 , aux_vert , np.ones((len(aux_face[0:,]))) , 1.5 )

            vert_and_face_arrays_to_text_and_mesh( mol_name , new_vertex , new_face.astype(int), input_suffix +
                                                  '_adj_'+ str(i), dens=dens, Self_build=True)

            aux_face , aux_vert = new_face.copy()- 1 , new_vertex.copy()

        adj_grid = Grid_loader( mol_name , dens , input_suffix + '_adj_' + str(N_ref) )
        
    if not calculate_phi:
        phi , dphi = 0. , 0. 
        
        return phi , dphi , adj_grid

    adj_face_array = np.transpose(adj_grid.leaf_view.elements)
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)

    #adj_el_pos = check_contained_triangles_alternative_2(coarse_grid , adj_grid , N_ref )
    
    potential.dirichl_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space ,
                                                           mesh_info.phi_order)
    potential.neumann_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space ,
                                                           mesh_info.phi_order) 
    potential.dual_to_dir_s_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space ,
                                                           mesh_info.phi_order)
        
    phi , dphi , it_count = adjoint_equation( potential.dirichl_space_phi , 
                        potential.neumann_space_phi , potential.dual_to_dir_s_phi)
    
    potential.phi , potential.dphi = phi , dphi
        
    if return_grid:
        return phi , dphi , adj_grid
    
    
    return phi , dphi

def main( name , dens , input_suffix , output_suffix , percentaje ,  N ,  N_ref ,
             smooth = True , Mallador = 'NanoShaper',refine = True, Use_GAMer = True , sphere=False , Estimator = 'E_u',
             x_q = None, q=None , r = np.nan , fine_vert_array = None):
    '''
    Calculates the solvation energy and refines the mesh.
    Params:
    name  : Name of the molecule. Must have a Molecule/{name} folder with the .pqr if not running sphere cases
    dens  : Mesh density. Set this to 0 if using sphere_cases
    input_suffix  : Mesh to be refined. If this is equal to "-0" the mesh will be build using MSMS
    output_suffix : suffix of the refined mesh.
    percentaje    : Percentaje of elements to be refined, which absolute error contribution is
                    less than percentaje*sum of the total absoulute error
    N             : Number of points used in the Gauss cuadrature. Can be {1,7,13,17,19,25,37,48,52,61,79}
    N_ref         : Number of UNIFORM refinements used to calculate phi. 
    smooth        : Smooths the mesh using a 40.0 [el/A2] grid, must be created before with MSMS
    refine        : True if the mesh is going to be refined, False if not
    Mallador      : Can be MSMS or NanoShaper
    Use_GAMer     : True if using GAMer or False if not. Read the README to use this.
    sphere        : True if the mesh is a spherical grid, and False if not
    Estimator     : E_phi or E_u
    x_q           : For the SPHERE case, can be defined in np.array([N_q , 3]) format
    q             : For the SPHERE case, can be defined in np.array([N_q , 1]) format
    r             : For the SPHERE case, this is the radius of the sphere
    
    This function gives the following output
    S_trad : Solvation energy using bempp potential operators
    S_Ap   : I
    S_Ex   : II
    N_elements : Number of elements of the input_suffix grid
    N_El_adj   : Number of elements of the adjoint grid
    total_solving_time : Time needed to create the operators and solving.
    S_trad_time        : Time needed to calculate S_trad, not counting total_solving_time
    S_Ap_time          : Time needed to calculate I, not counting total_solving_time
    S_Ex_time          : Time needed to calculate II
    operators_time_U   : Time needed to build the associated U operators [Not used]
    assembly_time_U    : Time needed to assembly the blocked operator    [Not used]
    solving_time_U     : Time needed to solve the U system
    it_count_U         : Number of iterations to solve the U system
    '''
    
    if sphere:
        
        x_q , q = saved_sphere_distributions(name , r)
        
        if x_q == None or q==None or r == np.nan:
            print('x_q, q or sphere radius where not defined.')
        
        if Estimator ==   'E_u':
            
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
            print('Total solving time for U and dU: {0:.4f}'.format(total_solving_time))
            init_time = time.time()
            S_trad = S_trad_calc_R( potential.dirichl_space_u, potential.neumann_space_u , U , dU ) 
            S_trad_time = time.time()-init_time
            print('Measured time to obtain S_trad : {0:.4f}'.format(S_trad_time))
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
                    print(np.min(aux_face_array))

                    vert_and_face_arrays_to_text_and_mesh( name , aux_vert_array , new_face_array.astype(int)[:] ,
                                                           output_suffix, dens , Self_build=True)

                    grid = Grid_loader( name , dens , output_suffix , 'Self')

            return ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , total_solving_time , S_trad_time , S_Ap_time ,
                     S_Ex_time , operators_time_U , assembly_time_U , solving_time_U , it_count_U )
            
        elif Estimator == 'E_phi':
            
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

            face_array = np.transpose(grid.leaf_view.elements)
            vert_array = np.transpose(grid.leaf_view.vertices)

            init_time = time.time()
            U, dU , operators_time_U , assembly_time_U , solving_time_U , it_count_U = U_tot_boundary(grid
                                                        , return_time = True , save_plot=True)
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

                    grid = Grid_loader( name , dens , output_suffix , 'Self')

            return ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , total_solving_time , S_trad_time , S_Ap_time ,
                     S_Ex_time , operators_time_U , assembly_time_U , solving_time_U , it_count_U )
        
    else:
        
        if Estimator ==   'E_u':
            
            mesh_info.mol_name     = name
            mesh_info.mesh_density = dens
            mesh_info.suffix       = input_suffix
            mesh_info.path         = os.path.join('Molecule' , mesh_info.mol_name)

            mesh_info.q , mesh_info.x_q = run_pqr(mesh_info.mol_name)
            
            init_time_0 = time.time()
            mesh_info.u_space , mesh_info.u_order     = 'DP' , 0
            mesh_info.phi_space , mesh_info.phi_order = 'P' , 1
            mesh_info.u_s_space , mesh_info.u_s_order = 'P' , 1
            
            if input_suffix == '-0' and not sphere:
                grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density , mesh_info.suffix , Mallador, GAMer = False)
            else:
                grid = Grid_loader( mesh_info.mol_name , mesh_info.mesh_density ,  mesh_info.suffix , 'Self')
            t_0 = time.time() - init_time_0
            print('Total time to create and load the mesh {0:.6f} [s]'.format(t_0))

            init_time_solving_out = time.time()
            U, dU , spaces_time_U , operators_time_U , assembly_time_U , GMRES_time_U , UdU_time, it_count_U = U_tot_boundary(grid
                                                        , return_time = True , save_plot=True , tolerance = 1e-5)
            total_solving_time = time.time() - init_time_solving_out
            print('Total solving time for U and dU: {0:.4f}'.format(total_solving_time))
            ################

            init_time_S_trad = time.time()
            S_trad   = S_trad_calc_R( potential.dirichl_space_u, potential.neumann_space_u , U , dU ) 
            t_S_trad = time.time()-init_time_S_trad
            print('Measured time to obtain S_trad : {0:.4f}'.format(t_S_trad))
            ################

            init_time_S_Ap = time.time()
            S_Ap , S_Ap_i = delta_G_tent_Pool( grid , U.coefficients.real , dU.coefficients.real , mesh_info.u_space ,
                                           mesh_info.u_order  , N )
            S_Ap_time = time.time() - init_time_S_Ap
            print('Time to calculate S_ap: {0:.2f}'.format(S_Ap_time))

            #################

            [S_Ex , S_Ex_i , it_count_phi , N_El_adj , flat_ref_time_adj , spaces_time_adj , operators_time_adj , 
             #(mol_name , grid  , dens , input_suffix , N ,save_energy_plot=False  , test_mode = False , return_times = False)
            matrix_time_adj , GMRES_time_adj , phidphi_time , S_Ex_time]  = S_Exact_in_Adjoint_Mesh_with_N_Ref_Pool(name , grid , dens , input_suffix , N , N_ref 
                                                    , save_energy_plot=False , test_mode=True , return_times = True)
            total_phi_time = matrix_time_adj + GMRES_time_adj + phidphi_time
            print('Time to solve the adjoint and dependencies: {0:.2f}'.format(total_phi_time))
            print('Time to calculate S_ex: {0:.2f}'.format(S_Ex_time))
            ################
            
            total_phi_time = matrix_time_adj + GMRES_time_adj + phidphi_time
            print('Time to solve the adjoint and dependencies: {0:.2f}'.format(total_phi_time))
            
            init_time_E = time.time()
            const_space = bempp.api.function_space(grid,  "DP", 0)            
            S_Ap_bempp = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ap_i[:,0])
            S_Ex_bempp    = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ex_i)

            dif =S_Ap_i[:,0]-S_Ex_i

            dif_F = bempp.api.GridFunction(const_space, fun=None, coefficients=np.abs(dif) )

            bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_N_ref_{2:d}.vtk'.format( dens,
                                                                input_suffix , N_ref )
                             , grid_function = dif_F , data_type = 'element')
            E_time = time.time()-init_time_E
      
            init_time_ref = time.time()
            if True: #Marked elements      
                init_time_status = time.time()
                face_array = np.transpose(grid.leaf_view.elements) + 1    
                status = value_assignor_starter(face_array , np.abs(dif) , percentaje)
                const_space = bempp.api.function_space(grid,  "DP", 0)
                Status    = bempp.api.GridFunction(const_space, fun=None, coefficients=status)
                status_time = time.time()-init_time_status
                bempp.api.export('Molecule/' + name +'/' + name + '_{0}{1}_Marked_elements_{2}.vtk'.format( 
                                                dens, input_suffix , N_ref )
                             , grid_function = Status , data_type = 'element')

            face_array = np.transpose(grid.leaf_view.elements)+1
            vert_array = np.transpose(grid.leaf_view.vertices)

            N_elements = len(face_array)

            if refine:
                init_time_mesh_refiner = time.time()
                new_face_array , new_vert_array = mesh_refiner(face_array , vert_array , np.abs(dif) , percentaje )
                mesh_refiner_time  = time.time()-init_time_mesh_refiner
                if smooth:
                    init_time_smoothing = time.time()
                    #fine_vert_array = text_to_list(name , '_40.0-0' , '.vert' , info_type=float)                
                    aux_vert_array  = smoothing_vertex( new_vert_array , fine_vert_array )
                    smoothing_time      = time.time()-init_time_smoothing
                elif not smooth:

                    aux_vert_array = new_vert_array.copy()

                if Use_GAMer:
                    init_time_GAMer = time.time()
                    new_face_array , aux_vert_array = Improve_Mesh(new_face_array , aux_vert_array , mesh_info.path , 
                                                                  mesh_info.mol_name+ '_' + str(dens) + output_suffix )

                    vert_and_face_arrays_to_text_and_mesh( name , aux_vert_array , new_face_array.astype(int)[:] ,
                                                           output_suffix, dens , Self_build=True)

   
            
                    grid = Grid_loader( name , dens , output_suffix ,'Self')
                    GAMer_time      = time.time()-init_time_GAMer
            t_ref =time.time()- init_time_ref
            
            times = np.array([ t_0 , spaces_time_U , operators_time_U , assembly_time_U ,
                              GMRES_time_U , UdU_time, t_S_trad, S_Ap_time, flat_ref_time_adj ,spaces_time_adj , 
                              operators_time_adj , matrix_time_adj , GMRES_time_adj , phidphi_time , 
                              S_Ex_time , E_time , status_time , mesh_refiner_time , smoothing_time  ,
                              GAMer_time ])
            
            return ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , it_count_U  , times )
            #return ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , total_solving_time , S_trad_time , S_Ap_time ,
            #         S_Ex_time , operators_time_U , assembly_time_U , solving_time_U , it_count_U , t_ref)
        
        elif Estimator == 'E_phi':
            
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
                                                        , return_time = True , save_plot=True, tolerance = 1e-5)
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
            init_time = time.time()
            
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

                    #fine_vert_array = text_to_list(name , '_40.0-0' , '.vert' , info_type=float)                
                    aux_vert_array  = smoothing_vertex( new_vert_array , vert_array, fine_vert_array )

                elif not smooth:

                    aux_vert_array = new_vert_array.copy()

                if Use_GAMer:

                    new_face_array , aux_vert_array = Improve_Mesh(new_face_array , aux_vert_array , mesh_info.path , 
                                                                  mesh_info.mol_name+ '_' + str(dens) + output_suffix )

                    vert_and_face_arrays_to_text_and_mesh( name , aux_vert_array , new_face_array.astype(int)[:] ,
                                                           output_suffix, dens , Self_build=True)

                    grid = Grid_loader( name , dens , output_suffix ,'Self')
                
            t_ref =time.time()- init_time

            return ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , total_solving_time , S_trad_time , S_Ap_time ,
                     S_Ex_time , operators_time_U , assembly_time_U , solving_time_U , it_count_U , t_ref)
            

def suffixes(N_ref):
    '''
    Defines the suffixes for the .msh/.vert/.off/.face/.... for the molecule
    Params
    N_ref   : Number of refinements
    Returns
    An array with N_ref+1 elements
    '''
    suf = np.array(['-0'])
    
    for i in (np.arange(N_ref+1)[1:]*(-1)).astype(str):
        suf = np.vstack([suf,i])
        
    suf = np.reshape(suf,[1,N_ref+1])[0]
    
    return suf

def create_txt_file(txt_name , headings = 'Molecule & density  & Times-Refinated & Elem & Adj Elem & S_Trad & S_Ap & S_Ex ' 
                    + '& Iterations_U & ' ,times = True):

    # Headings can also be
    # 'Molecule & density  & Times-Refinated & Elem & Adj Elem & S_Trad & S_Ap & S_Ex ' 
    #               + ' solving_time & S_trad_time & S_Ap_time & S_Ex_time & Operators_Build_Time_U & Assembly_Time_U & '
    #               + 'Linear_Sys_Solv_time_U & Iterations_U % t_ref'
    '''
    Creates the txt file in Results/.
    Input
    txt_name  : Name of the text file
    headings  : Headings of the text file, separated by  & .
    Output
    None
    '''
    if times:
        headings= headings +   't_0 & spaces_time_U & operators_time_U & assembly_time_U & ' 
        headings= headings +'GMRES_time_U & UdU_time & S_trad_time & S_Ap_time & Flat ref time & spaces_time_adj & '
        headings= headings +'operators_time_adj & matrix_time_adj & GMRES_time_adj & phidphi_time & ' 
        headings= headings +'S_Ex_time , E_time , status_time , mesh_refiner_time , smoothing_time  & GAMer_time '
        
    if txt_name[-4:]!='.txt':
        print(txt_name[-4:])
        txt_name = txt_name+'.txt'
    
    Text_file = open( 'Results/'+txt_name , 'w+' )
    Text_file.write( headings + '\n')
    Text_file.close()   

    return None

def append_txt_file(txt_name , text):
    '''
    Add info to Results/txt_name.txt
    Input 
    txt_name : Name of the text file
    text     : Line to be added to the txt file
    Output
    None
    '''
    if text[-2:]!='\n':
        text=text + '\n'
        
    if txt_name[-4:]!='.txt':
        print(txt_name[-4:])
        txt_name = txt_name+'.txt'
    
    Text_file = open( 'Results/'+txt_name , 'a' )  
    Text_file.write( text )             
    Text_file.close()
    
    return None

def parse_aux_func():
    
    import argparse
    from datetime import datetime
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-Molecule" , "-M"   , help="Molecule name in Molecule/mol_name/.")
    parser.add_argument("-Density"  , "-D"   , help="Starting mesh density", type=float)
    parser.add_argument("-N_it"     , "-n"   , help="Number of iterations/consecutive mesh refinements",type=int)
    parser.add_argument("-E"        , help="Error estimator. Can be E_u or E_phi",type=str)

    #Optional input here
    parser.add_argument("-N_ref"      , help="[Optional] Number of elements for the phi mesh. Default=0" ,
                        default=0,type=int)
    parser.add_argument("-Percentaje" , "-P" , help="[Optional] Percentaje of the elements to be refined, acording to their error contribution",
                                        default=0.1, type=float)
    parser.add_argument("-sphere"            , help="[Optional] True if the geometry is a sphere, and false if not" ,
                                        action="store_true", default=False)
    parser.add_argument("-Results_name"      , help="[Optional] Result text file name in Results/." ,
                                        default="Results_{0}_{1}_{2}.txt".format(datetime.now().year ,
                                        datetime.now().month, datetime.now().day), type=str)
    parser.add_argument("-radius"      , help="[Optional] Radius of the spherical geometry" ,
                                        default= 1.0 , type=float)
    parser.add_argument("-N_Gauss"     , help="[Optional] Number of points for the Gauss Quadrature" ,
                                        default= 25 , type=int)

    args = parser.parse_args()
    
    return args

def saved_sphere_distributions(name , r):
    '''
    Useful when running spherical grids.
    Inputs
    name  : Name of the distribution. Can be 'sphere_cent', 'sphere_offcent' or 'charge-dipole'.
    r     : Sphere radius
    Returns x_q , q 
    '''
        
    if name == 'sphere_cent':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,  1.E-12 ]]  )
        q = np.array( [1.] )
    
    if name == 'sphere_offcent':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,   r/2. ]]  )
        q = np.array( [1.] )
        
    if name == 'charge-dipole':
        x_q = np.array( [[  1.E-12 ,  1.E-12 ,  0.62 ],
                 [  1.E-12 ,  0.62*np.cos(np.pi*1.5 + 5.*np.pi/180.) ,
                                                      0.62*np.sin(np.pi*1.5 + 5.*np.pi/180. ) ] ,
                 [  1.E-12 ,  0.62*np.cos(np.pi*1.5 - 5.*np.pi/180.) ,
                                                      0.62*np.sin(np.pi*1.5 - 5.*np.pi/180. )  ]
                       ] )
        q = np.array( [1. , 1. , -1.]) 
    
    return x_q , q

def Delta_G_tent_int(c , face_array , vert_array , normals , Areas , U , dU , N , X_K , W):
    
    face = face_array[c]
    
    v1 , v2 , v3 = vert_array[face[0]] , vert_array[face[1]] , vert_array[face[2]]
            
    X = evaluation_points_and_weights_new(v1,v2,v3 , N , X_K , W)
                
    u_c_face_2  = np.array(list(map( lambda x_i : u_s_Teo(x_i)              , X )))
    du_c_face_2 = np.array(list(map( lambda x_i : du_s_Teo(x_i, normals[c]) , X )))

    U_local , dU_local = U[c] , dU[c]

    I2 =  np.sum( U_local  * du_c_face_2 * W )
    I1 =  np.sum( dU_local * u_c_face_2  * W )

    return (I1-I2)* Areas[c]

import time


def delta_G_tent_2(grid , U , dU , U_space , U_order , N):
    # Works good
    face_array = np.transpose(grid.leaf_view.elements) 
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    if U_space == 'DP' and U_order == 0:

        mesh = trimesh.Trimesh(vertices=vert_array , faces=face_array)
        normals = mesh.face_normals
        Areas = mesh.area_faces

        quadrule = quadratureRule_fine(N)
        X_K , W  = quadrule[0].reshape(-1,3) , quadrule[1]
        
        #Integral_func = lambda c : Delta_G_tent_int(c , face_array , vert_array , normals , Areas
        #                                            ,U , dU , N , X_K , W)
        init_time_integral = time.time()
        Integral      = np.array(list(map( lambda c : Delta_G_tent_int(c , face_array ,
                vert_array , normals , Areas ,U , dU , N , X_K , W) , np.arange(len(face_array)) )))
        
        

    Solv_Ex_i = K * ep_m * Integral 
    
    print('Aproximated Solvation energy {0:.6f}'.format(np.sum(Solv_Ex_i)))
    
    return np.sum(Solv_Ex_i) , np.reshape(Solv_Ex_i , (-1,1))

def delta_G_tent_3(grid , U , dU , U_space , U_order , N):
    
    face_array = np.transpose(grid.leaf_view.elements) 
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    init_time_x_i = time.time()
    
    quadrule = quadratureRule_fine(N)
    X_K , W  = quadrule[0].reshape(-1,3) , quadrule[1]
    
    x_i = np.concatenate(np.array( [evaluation_points_and_weights_new(vert_array[face[0]] ,
                                                 vert_array[face[1]] ,
                                                 vert_array[face[2]] , N , X_K , W) for face in face_array ])
                        )
    print('Time to create grid points = {0:.4f}'.format(time.time()-init_time_x_i))
                         
    if U_space == 'DP' and U_order == 0:

        mesh = trimesh.Trimesh(vertices=vert_array , faces=face_array)
        normals = mesh.face_normals
        Areas = mesh.area_faces

        
        
        #Integral_func = lambda c : Delta_G_tent_int(c , face_array , vert_array , normals , Areas
        #                                            ,U , dU , N , X_K , W)
        
        
        
        
        u_c_init_time = time.time()
        u_c_func  = lambda c : u_s_Teo(x_i[c])
        u_c_face  = np.array([ u_c_func(c) for c in range(len(x_i))])
        du_c_func = lambda c : du_s_Teo(x_i[c] , normals[c//N] )
        u_c_face  = np.array([ u_c_func(c) for c in range(len(x_i))]).reshape((N,-1)) 
        du_c_face = np.array([du_c_func(c) for c in range(len(x_i))]).reshape((N,-1)) 
        print('u_c time {0:.4f}'.format(time.time()-u_c_init_time))
        
        integrals_init_time = time.time()
        I1 = np.array([ np.sum( u_c_face[:,i] * W) for i in range(len(face_array))  ]) * dU
        I2 = np.array([ np.sum(du_c_face[:,i] * W) for i in range(len(face_array))  ]) *  U
        
        Integral = (I1-I2)*Areas     
        print('Integrals time {0:.4f}'.format(time.time()-integrals_init_time))

    Solv_Ex_i = K * ep_m * Integral 
    
    print('Aproximated Solvation energy {0:.6f}'.format(np.sum(Solv_Ex_i)))
    
    return np.sum(Solv_Ex_i) , np.reshape(Solv_Ex_i , (-1,1))



def delta_G_tent_Pool(grid , U , dU , U_space , U_order , N):
    
    face_array = np.transpose(grid.leaf_view.elements) 
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    if U_space == 'DP' and U_order == 0:

        mesh = trimesh.Trimesh(vertices=vert_array , faces=face_array)
        normals = mesh.face_normals
        Areas = mesh.area_faces

        quadrule = quadratureRule_fine(N)
        X_K , W  = quadrule[0].reshape(-1,3) , quadrule[1]
        
        #Integral_func = lambda c : Delta_G_tent_int(c , face_array , vert_array , normals , Areas
        #                                            ,U , dU , N , X_K , W)
        def func(c):
            return Delta_G_tent_int(c , face_array ,
                        vert_array , normals , Areas ,U , dU , N , X_K , W)
        
        Integral      = np.array(Pool().map( func , np.arange(len(face_array)) ))
        #Pool().terminate()
        #Pool().join()
        

    Solv_Ex_i = K * ep_m * Integral 
    #Pool().clear()
    #print('Aproximated Solvation energy {0:.6f}'.format(np.sum(Solv_Ex_i)))
    
    return np.sum(Solv_Ex_i) , np.reshape(Solv_Ex_i , (-1,1))

def S_Exact_in_Adjoint_Mesh_with_no_Ref_Pool(mol_name , grid  , dens , input_suffix , N ,
                                            save_energy_plot=False  , return_times = True):
    
    face_array = np.transpose(grid.leaf_view.elements)
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    init_time_spaces_adj  = time.time()
    dirichl_space_phi = bempp.api.function_space(grid,  mesh_info.phi_space , mesh_info.phi_order)
    neumann_space_phi = bempp.api.function_space(grid,  mesh_info.phi_space , mesh_info.phi_order) 
    dual_to_dir_s_phi = bempp.api.function_space(grid,  mesh_info.phi_space , mesh_info.phi_order)
    spaces_time_adj       = time.time()-init_time_spaces_adj
    
    phi , dphi , operators_time_adj , matrix_time_adj , GMRES_time_adj , phidphi_time , it_count = adjoint_equation( 
                    dirichl_space_phi , neumann_space_phi , dual_to_dir_s_phi , save_plot = True )
    init_time_S_Ex = time.time()
    S_Ex , S_Ex_i= Exact_aproach_with_u_s_Teo_Pool( face_array , vert_array ,
                                                   phi.coefficients.real , dphi.coefficients.real , N )
    S_Ex_time      = time.time()-init_time_S_Ex
    print('Exact solvation energy solved in {0:.5f} [s]'.format(S_Ex_time))
    N_el_adjoint   = len(face_array)
        
    print('Exact solvation energy {0:.5f} [kcal/kmol]'.format(S_Ex))
    
    return (S_Ex , S_Ex_i.reshape((-1,1)) , it_count , N_el_adjoint , spaces_time_adj , operators_time_adj , 
                matrix_time_adj , GMRES_time_adj , phidphi_time , S_Ex_time )

def S_ex_integrate_face(c , face_array , vert_array , normals, phi , dphi , X_K , W  , N):
    
    f1 , f2 , f3 = face_array[c]
    v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]          
        
    A = matrix_lineal_transform( v1 , v2 , v3 )    
    
    x_i = evaluation_points_and_weights_new(v1,v2,v3 , N , X_K , W)
        
    phi_a  = np.array([  phi[f1] ,  phi[f2] ,  phi[f3] ])
    dphi_a = np.array([ dphi[f1] , dphi[f2] , dphi[f3] ])
    
    phi_local  = np.array(list(map( lambda x_g : local_f( x_g , A ,  phi_a  , 1), x_i)))
    dphi_local = np.array(list(map( lambda x_g : local_f( x_g , A , dphi_a  , 1), x_i)))
            
    u_s_local  = np.array(list(map( lambda x_g : u_s_Teo( x_g )                 , x_i)))
    du_s_local = np.array(list(map( lambda x_g : du_s_Teo( x_g , normals[c] )   , x_i)))
            
    Integrate = np.sum( (dphi_local * u_s_local  -  phi_local * du_s_local) * W )    
                
    return Integrate


def Exact_aproach_with_u_s_Teo_Pool( face_array , vert_array , phi , dphi , N ):
    
    mesh    = trimesh.Trimesh(vertices=vert_array , faces=face_array)
    normals = mesh.face_normals
    Areas   = mesh.area_faces
    
    quadrule = quadratureRule_fine(N)
    X_K , W  = quadrule[0].reshape(-1,3) , quadrule[1]
    
    def integrate_i(c):
        return S_ex_integrate_face(c , face_array , vert_array  , normals, phi ,
                                   dphi, X_K , W , N)
    #def testirijillo(i):
    #    return 1./(i+10.)**5
    
    Integrates = np.array(Pool().map( integrate_i , np.arange(len(face_array)) )) 
        
    Solv_Exact_i = K * Integrates * ep_m * Areas
            
    
    return np.sum(Solv_Exact_i) , Solv_Exact_i


def S_Exact_in_Adjoint_Mesh_with_N_Ref_Pool(mol_name , grid  , dens , input_suffix , N , N_ref ,
                                        save_energy_plot=False  , test_mode = False , return_times = False):
    
    face_array = np.transpose(grid.leaf_view.elements)
    vert_array = np.transpose(grid.leaf_view.vertices)
    
    aux_face = face_array.copy()
    aux_vert = vert_array.copy()
    print(np.shape(face_array))
    
    
    flat_ref_time_init = time.time()
    if N_ref == 0:
        adj_grid = grid
        flat_ref_time = time.time()-flat_ref_time_init
        
    elif N_ref>=1:
        flat_ref_time_init = time.time()
        
        aux_grid = grid
        i=1
        while i <= N_ref:

            #new_face , new_vertex = mesh_refiner(  aux_face +1 , aux_vert , np.ones((len(aux_face[0:,]))) , 1.5 )
            
            #aux_face , aux_vert   = new_face.copy() , new_vertex.copy()
            
            #aux_grid =  aux_grid.refine()
            
            mesh = trimesh.Trimesh(vertices=aux_vert , faces=aux_face)
            aux_mesh = mesh.subdivide()
            aux_vert , aux_face = aux_mesh.vertices , aux_mesh.faces

            i+=1
        
        flat_ref_time = time.time()-flat_ref_time_init

    if N_ref>=1:
        print('The flat refinement was done in {0:.2f} seconds'.format(flat_ref_time))
    
    
    saving_refined_mesh_init = time.time()
    vert_and_face_arrays_to_text_and_mesh( mol_name , aux_vert , aux_face+1 ,
                                                  input_suffix +
                                                  '_adj_'+ str(i-1), dens=dens, Self_build=True)
    
    saving_refined_mesh_time = time.time() - saving_refined_mesh_init
    #adj_grid = aux_grid
    adj_grid  = Grid_loader( mol_name , dens , input_suffix + '_adj_' + str(N_ref) , 'Self')       
        
    print('The grid was uniformly refined in {0:.2f} seconds'.format(flat_ref_time))
    
    adj_face_array = np.transpose(adj_grid.leaf_view.elements) 
    adj_vert_array = np.transpose(adj_grid.leaf_view.vertices)
    print(np.shape(adj_face_array))

    #adj_el_pos = check_contained_triangles__(grid , adj_grid)
    
    #print(len(adj_face_array))
    #print(adj_el_pos)
    init_time_spaces_adj  = time.time()
    dirichl_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order)
    neumann_space_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order) 
    dual_to_dir_s_phi = bempp.api.function_space(adj_grid,  mesh_info.phi_space , mesh_info.phi_order)
    spaces_time_adj   = time.time()-init_time_spaces_adj
    
    phi , dphi , operators_time_adj , matrix_time_adj , GMRES_time_adj , phidphi_time , it_count = adjoint_equation( dirichl_space_phi , neumann_space_phi , dual_to_dir_s_phi , 
                                            save_plot = True , suffix = '_'+str(N_ref) )
    init_time_S_Ex = time.time()
    
    S_Ex , S_Ex_j = Exact_aproach_with_u_s_Teo_Pool( adj_face_array , adj_vert_array ,
                                                    phi.coefficients.real , dphi.coefficients.real , N )
    #S_Ex    , rearange_S_Ex_i  , S_Ex_i , _= Exact_aproach_with_u_s_Teo( adj_face_array , adj_vert_array , phi , dphi , N , 
    #                                             grid_relation = adj_el_pos , return_values_on_Adj_mesh = True)
        
    N_el_adjoint = len(adj_face_array)
    
    
    if test_mode:
        
        adj_el_pos = (np.arange(len(adj_face_array))/4).astype(int)
        
        const_space = bempp.api.function_space(adj_grid,  "DP", 0)
        counter     = bempp.api.GridFunction(const_space, fun=None, coefficients=adj_el_pos )
        bempp.api.export('Molecule/' + mol_name +'/' + mol_name + '_{0}{1}_Counter_{2}.msh'.format( 
                                        dens, input_suffix , N_ref )
                     , grid_function = counter , data_type = 'element')
        
        const_space_u = bempp.api.function_space(grid,  "DP", 0)
        counter_u     = bempp.api.GridFunction(const_space_u, fun=None, coefficients=np.arange(len(face_array)) )
        bempp.api.export('Molecule/' + mol_name +'/' + mol_name + '_{0}{1}_Counter_original_{2}.msh'.format( 
                                        dens, input_suffix , N_ref )
                     , grid_function = counter_u , data_type = 'element')
    
    if save_energy_plot:
        const_space = bempp.api.function_space(adj_grid,  "DP", 0)
        S_Ex_BEMPP  = bempp.api.GridFunction(const_space, fun=None, coefficients=S_Ex_i)
        bempp.api.export('Molecule/' + mol_name +'/' + mol_name + '_{0}{1}_S_Exact_{2}.msh'.format( 
                                        dens, input_suffix , N_ref )
                     , grid_function = S_Ex_BEMPP , data_type = 'element')
    
    
    print('Exact solvation energy {0:.5f} [kcal/kmol]'.format(S_Ex))
    rearange_S_Ex_i = np.sum( np.reshape(S_Ex_j , (-1,4**N_ref) )  , axis = 1)
    S_Ex_time      = time.time()-init_time_S_Ex
    
    
    return (S_Ex , rearange_S_Ex_i , it_count , N_el_adjoint , flat_ref_time , spaces_time_adj , operators_time_adj , 
                matrix_time_adj , GMRES_time_adj , phidphi_time , S_Ex_time )
