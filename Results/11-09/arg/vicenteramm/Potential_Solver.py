
# checked 19.09
# 11-09 added adjoint mesh

import bempp.api
import numpy as np
import os

from bempp.api.operators.potential import laplace as lp
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

from constants import mesh_info
from constants import values
from constants import *

from quadrature import *

from Mesh_Ref_V2 import *
from Grid_Maker_R2 import *

# --------------------------------------------------------------------------------

def zero_i(x, n, domain_index, result):
    result[:] = 0

def u_s_G(x,n,domain_index,result):
    global ep_m 
    result[:] =  1. / (4.*np.pi*ep_m)  * np.sum( mesh_info.q / np.linalg.norm( x - mesh_info.x_q, axis=1 ) )

def du_s_G(x,n,domain_index,result):
    global ep_m
    result[:] = -1./(4.*np.pi*ep_m)  * np.sum( np.dot( x-
                            mesh_info.x_q , n  )  * mesh_info.q / np.linalg.norm( x - mesh_info.x_q, axis=1 )**3 )

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



def U_and_U_Reac(dirichl_space , neumann_space , dual_to_dir_s):
    '''
    Computes reaction potential for a given space with the definition: U = U^Reac + U^Coulomb.    
    '''
    global ep_m , ep_s , k

    identity = sparse.identity(     dirichl_space, dirichl_space, dual_to_dir_s)
    slp_in   = laplace.single_layer(neumann_space, dirichl_space, dual_to_dir_s)
    dlp_in   = laplace.double_layer(dirichl_space, dirichl_space, dual_to_dir_s)

    slp_out  = modified_helmholtz.single_layer(neumann_space, dirichl_space, dual_to_dir_s, k)
    dlp_out  = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dual_to_dir_s, k)

    charged_grid_fun = bempp.api.GridFunction(dirichl_space, fun=q_times_G_L)
    zero_grid_fun    = bempp.api.GridFunction(neumann_space, fun=zero_i     )

    blocked = bempp.api.BlockedOperator(2, 2)
    blocked[0, 0] = 0.5*identity + dlp_in
    blocked[0, 1] = -slp_in
    blocked[1, 0] = 0.5*identity - dlp_out
    blocked[1, 1] = (ep_m/ep_s)*slp_out

    rhs = [charged_grid_fun, zero_grid_fun]

    sol, info,it_count = bempp.api.linalg.gmres( blocked, rhs , return_iteration_count=True , tol=1e-8)
    print("The linear system for U_tot was solved in {0} iterations".format(it_count))
    U , dU = sol
    
    U_s  = bempp.api.GridFunction(dirichl_space , fun =  u_s_G)
    dU_s = bempp.api.GridFunction(neumann_space , fun = du_s_G)
    
    U_R  =  U -  U_s
    dU_R = dU - dU_s
    
    return U, dU , U_R , dU_R

def S_trad_calc_R( dirichl_space, neumann_space , U , dU ):
    
    # En base a los puntos donde se encuentran las cargas, calculemos el potencial u_r y u_h
    # Esto luego de que podemos escribir la energia de solvatacion como
    # G_solv = Sum_i q_i *u_reac = Sum_i q_i * (u_h+u_r)           evaluado en cada carga.

    # Se definen los operadores
    slp_in_O = lp.single_layer(neumann_space, mesh_info.x_q.transpose()) 
    dlp_in_O = lp.double_layer(dirichl_space, mesh_info.x_q.transpose())

    # Y con la solucion de las fronteras se fabrica el potencial evaluada en la posicion de cada carga
    U_R_O = slp_in_O * dU  -  dlp_in_O * U

    # Donde agregando algunas constantes podemos calcular la energia de solvatacion S
    
    S     = K * np.sum(mesh_info.q * U_R_O).real
    print("Three Term Splitting Solvation Energy : {:7.8f} [kCal/mol] ".format(S) )
    
    return S


# --------------------------------------------------------------------------------

def S_Zeb_in_Adjoint_Mesh(mol_name , face_array , vert_array , dens , input_suffix , N):
    
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
    
    phi , dphi = adjoint_equation( dirichl_space_phi , neumann_space_phi , dual_to_dir_s_phi)
    
    S_Zeb    , S_Zeb_i = Zeb_aproach_with_u_s_Teo( adj_face_array , adj_vert_array , phi , dphi , N)
    
    rearange_S_Zeb_i = np.zeros((len(face_array),1))
    
    c=0
    for G_i in S_Zeb_i:
        rearange_S_Zeb_i[int(adj_el_pos[c])] += G_i
        c+=1
        
    #print(rearange_S_Zeb_i)
    
    return S_Zeb , rearange_S_Zeb_i


# ----------------------------------------------------------------
    

def q_times_G_L(x, n, domain_index, result):
    global ep_m
    result[:] = 1. / (4.*np.pi*ep_m)  * np.sum( mesh_info.q  / np.linalg.norm( x - mesh_info.x_q, axis=1 ) )

def adjoint_equation( dirichl_space , neumann_space , dual_to_dir_s):
    
    global ep_m , ep_s , k
    
    identity = sparse.identity(     dirichl_space, dirichl_space, dual_to_dir_s)
    slp_in   = laplace.single_layer(neumann_space, dirichl_space, dual_to_dir_s)
    dlp_in   = laplace.double_layer(dirichl_space, dirichl_space, dual_to_dir_s)
    slp_out  = modified_helmholtz.single_layer(neumann_space, dirichl_space, dual_to_dir_s, k)
    dlp_out  = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dual_to_dir_s, k)

    blocked = bempp.api.BlockedOperator(2, 2)
    blocked[0, 0] = 0.5*identity + dlp_in
    blocked[0, 1] = -slp_in
    blocked[1, 0] = 0.5*identity - dlp_out
    blocked[1, 1] = (ep_m/ep_s)*slp_out

    zero = bempp.api.GridFunction(dirichl_space , fun=zero_i)
    P_GL = bempp.api.GridFunction(dirichl_space, fun=q_times_G_L)
    rs_r = [P_GL , zero]

    sol_r, info,it_count = bempp.api.linalg.gmres( blocked, rs_r , return_iteration_count=True, tol=1e-8)
    print("The linear system for phi was solved in {0} iterations".format(it_count))
    phi_r , dphi_r = sol_r
    
    return phi_r , dphi_r

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
    print("Three Term Splitting Solvation Energy : {:7.8f} [kCal/mol] ".format(S) )
    
    return S

def S_Cooper_calc( face_array , vert_array , phi_r , dphi_r , U_Reac , dU_Reac , N):
    
    Solv_Cooper = np.zeros((len(face_array),1))
    c = 0
    for face in face_array:

        I1 = int_calc_i( face , face_array , vert_array , phi_r , mesh_info.phi_space 
                        , mesh_info.phi_order , ep_m*dU_Reac , mesh_info.u_space , mesh_info.u_order  , N)
        
        I2 = int_calc_i( face , face_array , vert_array , ep_m * dphi_r , mesh_info.phi_space
                        , mesh_info.phi_order , U_Reac , mesh_info.u_space , mesh_info.u_order  , N)
        Solv_Cooper[c] = I1-I2

        c+=1

    S_Cooper_i = Solv_Cooper
    S_Cooper = K*np.sum(Solv_Cooper )
    print('Cooper Solv = {0:10f} '.format(S_Cooper)) 
    
    return S_Cooper , S_Cooper_i

def S_Zeb_calc( face_array , vert_array , phi , dphi , u_s , du_s , N):
    
    Solv_Zeb = np.zeros((len(face_array),1))
    c = 0
    print(mesh_info.u_s_order)
    for face in face_array:

        I1 = int_calc_i( face , face_array , vert_array , phi , mesh_info.phi_space 
                        , mesh_info.phi_order , ep_m*(du_s) , mesh_info.u_s_space , mesh_info.u_s_order  , N)
        
        I2 = int_calc_i( face , face_array , vert_array ,  dphi , mesh_info.phi_space
                         , mesh_info.phi_order , ep_m*(u_s) , mesh_info.u_s_space , mesh_info.u_s_order , N)
        Solv_Zeb[c] = I2-I1

        c+=1
    Solv_Zeb_i = Solv_Zeb
    S_Zeb = K*np.sum(Solv_Zeb )
    print('Zeb Solv = {0:10f} '.format(S_Zeb)) 
    

    
    return S_Zeb , Solv_Zeb_i


def Zeb_aproach_with_u_s_Teo( face_array , vert_array , phi , dphi , N):
    
    normals = normals_to_element( face_array , vert_array )
    
    Solv_Zeb = np.zeros((len(face_array),1))
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
            
        Solv_Zeb[c] = (I1-I2)*Area

        c+=1
    Solv_Zeb_i = Solv_Zeb
    S_Zeb = K*np.sum(Solv_Zeb )
    print('Zeb Solv = {0:10f} '.format(S_Zeb)) 
    
    return S_Zeb , Solv_Zeb_i

def u_s_Teo( x ):
    
    return (1. / (4.*np.pi*ep_m) ) * np.sum( mesh_info.q / np.linalg.norm( x - mesh_info.x_q, axis=1 ) )
    
    #result[:] =  C / (4.*np.pi*ep_m)  * np.sum( mesh_info.q / np.linalg.norm( x - mesh_info.x_q, axis=1 ) )

def du_s_Teo(x,n):
    
    return -1./(4.*np.pi*ep_m)  * np.sum( np.dot( x-
                            mesh_info.x_q , n)  * mesh_info.q / np.linalg.norm( x - mesh_info.x_q, axis=1 )**3 )

def normals_to_element( face_array , vert_array ):

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
