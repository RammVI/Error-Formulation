
import bempp.api, numpy as np
from math import pi
import os

from Grid_Maker import *

# 1) Se elabora y carga la malla

# 1------------------------------------------------------------------------------------------------
# Se definen las variables generales de la malla, vale decir
# mol_name       : Nombre/Sigla de la molécula
# stern_tickness : Espesor de stern_layer
# probe_radius   : Espesor de la probeta para construir la SES
# mesh_density   : Densidad de malla en el/A^2 
# min_area       : Área mínima, agregar este en el caso de usar MSMS mejora la solución

mol_name = 'arg' ; stern_thickness = 0 ; probe_radius = 1.4
mesh_density = 8.0 ; min_area = 0 

# 1------------------------------------------------------------------------------------------------
# Aquí se utilizan las funciones de la librería Grid_Maker para elaborar la malla, se trabaja
# a partir de archivos en formato .pdb

pdb_to_pqr(mol_name , stern_thickness)
pqr_to_xyzr(mol_name , stern_thickness , method = 'amber' )
xyzr_to_msh(mol_name , mesh_density , probe_radius , stern_thickness , min_area , Mallador = 'MSMS')

# 1------------------------------------------------------------------------------------------------
# Esta rutina crea una lista en formato de texto para las áreas de cada elemento, utilizadas más adelante.

path = 'Molecule/'+mol_name
triangle_areas( path , mol_name , mesh_density )

# 1------------------------------------------------------------------------------------------------
# Se especifica el visualizador de las soluciones
bempp.api.set_ipython_notebook_viewer()
bempp.api.PLOT_BACKEND = "ipython_notebook"


# 1------------------------------------------------------------------------------------------------
# Aquí se carga la malla al kernel
path = 'Molecule/'+mol_name+'/'
grid_name_File =  os.path.join(path,mol_name + '_'+str(mesh_density)+'.msh')
grid = bempp.api.import_grid(grid_name_File)

# 1------------------------------------------------------------------------------------------------
# Grafica de la malla
# grid.plot()

# 2) Parámetros del sistema
# 2------------------------------------------------------------------------------------------------
# Se definen las variables físicas del problema
# ep_m : Permisividad del medio en la molécula
# ep_s : Permisividad del medio en el solvente
# k    : Kappa, es la inversa del largo de Debye.
ep_m = 4.
ep_s = 80.
k = 0.125

# e_c : Carga del protón
# k_B : Constante de Stefan Boltzmann
# T   : Temperatura promedio, en grados Kelvin
# C   : Constante igual a e_c/(k_B*T). Para el caso se utilizará como 1 y se agregará una cte a la QoI
e_c = 1.60217662e-19 # [C] - proton charge
k_B = 1.38064852e-23 # [m2 kg s-2 K-1]
T   = 298. # [K]  
C = 1. # 

# 2------------------------------------------------------------------------------------------------
# Se crea una lista vacía para almacenar las cargas desde el archivo .pqr y se almacenan los valores 
# y posiciones en q y x_q respectivamente
q, x_q = np.empty(0), np.empty((0,3))
pqr_file = os.path.join(path,mol_name+'.pqr')
charges_file = open( pqr_file , 'r').read().split('\n')

for line in charges_file:
    line = line.split()
    if len(line)==0: continue
    if line[0]!='ATOM': continue
    q = np.append( q, float(line[8]))
    x_q = np.vstack( ( x_q, np.array(line[5:8]).astype(float) ) )   
    
# 3) Definicion de u_s y su derivada normal

def zero_i(x, n, domain_index, result):
    result[:] = 0
    
def u_s_G(x,n,domain_index,result):
    global q,x_q,ep_m,C
    result[:] = C / (4.*np.pi*ep_m)  * np.sum( q / np.linalg.norm( x - x_q, axis=1 ) )
    
def du_s_G(x,n,domain_index,result):
    global q,x_q,ep_m,C
    result[:] = -C/(4.*np.pi*ep_m)  * np.sum( np.dot( x-
                            x_q , n)  * q / np.linalg.norm( x - x_q, axis=1 )**3 )
    
# 4) Resolución del problema

# Se cargan los operadores respectivos
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

# Se crean los espacios de la solución
# DP - i : Polinomio de orden i discontinuo entre elementos
#  P - i : Polinomio de orden i continuo entre elementos
# Para nuestro caso utilizaremos DP-0, es decir, valores constantes por elemento
dirichl_space = bempp.api.function_space(grid,  "DP", 0)
neumann_space = bempp.api.function_space(grid,  "DP", 0) 
dual_to_dir_s = bempp.api.function_space(grid,  "DP", 0)

# Se crea u_s y su derivada normal du_s en la frontera
u_s  = bempp.api.GridFunction(dirichl_space, fun=u_s_G)
du_s = bempp.api.GridFunction(neumann_space, fun=du_s_G) 

# Se crean los operadores asociados al sistema, que dependen de los espacios de las soluciones
# identity : I  : matriz identidad
# dlp_in : K_in : Double-Layer operator para la region interior
# slp_in : V_in : Single-Layer operator para la region interior
# _out : Mismos operadores pero para la región exterior, con k=kappa=0.125
identity = sparse.identity(     dirichl_space, dirichl_space, dual_to_dir_s)
slp_in   = laplace.single_layer(neumann_space, dirichl_space, dual_to_dir_s)
dlp_in   = laplace.double_layer(dirichl_space, dirichl_space, dual_to_dir_s)
slp_out  = modified_helmholtz.single_layer(neumann_space, dirichl_space, dual_to_dir_s, k)
dlp_out  = modified_helmholtz.double_layer(dirichl_space, dirichl_space, dual_to_dir_s, k)


# 4.1 Solución armónica ---------------------------------------------------------------

# Dada por V_in du_s = (1/2+K_in)u_h = -(1/2+K_in)u_s (BC)
sol, info,it_count = bempp.api.linalg.gmres( slp_in, -(dlp_in+0.5*identity)*u_s , return_iteration_count=True, tol=1e-4)
print("The linear system for du_h was solved in {0} iterations".format(it_count))

u_h = -u_s
du_h = sol



# 4.2 Solución regular -----------------------------------------------------------------
k_p = k # 0.125

# Se crea la matriz / Lado izquierdo de la ecuación
# | ( I/2 + K_L-in  )     (      -V_L-in     ) |  u_r  = 0
# | ( I/2 - K_Y-out )     ( ep_m/ep_s V_Y-out) | du_r  = ep_m/ep_s V_Y-out*(du_s+du_h)  (BC)
blocked = bempp.api.BlockedOperator(2, 2)
blocked[0, 0] = 0.5*identity + dlp_in
blocked[0, 1] = -slp_in
blocked[1, 0] = 0.5*identity - dlp_out
blocked[1, 1] = (ep_m/ep_s)*slp_out

# Se crea el lado derecho de la ecuación 
zero = bempp.api.GridFunction(dirichl_space, fun=zero_i)
rhs = [ zero ,  -slp_out *(ep_m/ep_s)* (du_s+du_h)]

# Y Finalmente se resuelve para u_r y du_r
sol, info,it_count = bempp.api.linalg.gmres( blocked , rhs, return_iteration_count=True, tol=1e-4)
print("The linear system for u_r and du_r was solved in {0} iterations".format(it_count))
u_r , du_r = sol

# 5) Calculo de la energía de solvatación

# Se importan los operadores para el potencial en el dominio de la molecula
from bempp.api.operators.potential import laplace as lp

# En base a los puntos donde se encuentran las cargas, calculemos el potencial u_r y u_h
# Esto luego de que podemos escribir la energía de solvatación como
# G_solv = Sum_i q_i *u_reac = Sum_i q_i * (u_h+u_r)           evaluado en cada carga.

# Se definen los operadores
slp_in_O = lp.single_layer(neumann_space, x_q.transpose()) 
dlp_in_O = lp.double_layer(dirichl_space, x_q.transpose())

# Y con la solución de las fronteras se fabrica el potencial evaluada en la posición de cada carga
u_r_O = slp_in_O * du_r  -  dlp_in_O * u_r
u_h_O = slp_in_O * du_h  +  dlp_in_O * u_s

terms =  u_r_O + u_h_O

# Donde agregando algunas constantes podemos calcular la energía de solvatacion S
K     = 0.5 * 4. * np.pi * 332.064 
S     = K * np.sum(q * terms).real
print("Three Term Splitting Solvation Energy : {:7.8f} [kCal/mol] ".format(S) )


# 6) Post-Procesamiento
# La idea de este método, es escribir la variable de forma tal que nos pueda entregar 
# el residuo en la QoI para cada elemento. Esta variable tiene la misma ecuación diferencial que
# el potencial total, con las mismas condiciones de frontera. El sistema está dado por

# Aquí se define parte del lado derecho (Que no es 0).
def q_times_G_L(x, n, domain_index, result):
    global q,x_q,ep_m,C , k
    result[:] = 1. / (4.*np.pi*ep_m)  * np.sum( q  / np.linalg.norm( x - x_q, axis=1 ) )

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

sol_r, info,it_count = bempp.api.linalg.gmres( blocked, rs_r , return_iteration_count=True, tol=1e-4)
print("The linear system for phi was solved in {0} iterations".format(it_count))
phi_r , dphi_r = sol_r

#phi_r.plot()

# Se asumirá que phi_h es 0, por lo que no se pondrá en este doc para ahorrar espacio y tiempo de cálculo.

# Se definen las integrales de superficie, válido para espacios DP0 unicamente.
def Integral_superficie( A , B ):
    A = A.real.coefficients
    B = B.real.coefficients
    Area_F = open(path+'triangleAreas_'+str(mesh_density)+'.txt').readlines()
    Areas  = np.zeros(len(Area_F))
    for i in range(len(Area_F)):
        Areas[i] = float(Area_F[i][:-1]) 
    return A*B*Areas

# Se abre el archivo de texto que guarda las áreas por elemento
area_txt = open('Molecule/{0:s}/triangleAreas_{1:s}.txt'.format(mol_name , str(mesh_density) ) ).read().split('\n')
areas    = [float(i) for i in area_txt[:-1]]

# Este código está creado para verificar que tan cierto es que <ep_m du_s , phi > corresponde a 
# la enegía de sovatación, comparado con la aproximación que entrega C. Cooper

# Energía de solvatación calculada tradicionalmente
S_phi     = K * np.sum(q * terms).real

# Aporte de cada elemento a la energía de solvatación según la expresión de Cooper
Solv_i_Cooper = K* (Integral_superficie( ep_m * (du_h+du_r) , phi_r ) - Integral_superficie( ep_m * (u_r+u_h) , dphi_r ) )

# Aporte de cada elemento a la energía de solvatación según la expresión de Zeb modificada
Solv_i_Zeb    = K * Integral_superficie( -ep_m * du_s , phi_r) + K * Integral_superficie( u_s , ep_m*dphi_r)
# Ojo que el último término fue agregado en la demostración I_1 que debe ser formalizada todavía

print('Forma Tradicional: {0:f} | Forma Cooper: {1:f} | Forma Zeb: {2:f}'.format(S_phi, Solv_i_Cooper.sum(), Solv_i_Zeb.sum()))

# Definamos una función para poder graficar esto, sobre un espacio DP-0

Solv_i_Cooper_Func = bempp.api.GridFunction(dirichl_space, fun=None, coefficients=Solv_i_Cooper)
Solv_i_Zeb_Func    = bempp.api.GridFunction(dirichl_space, fun=None, coefficients=Solv_i_Zeb)

# Lo bueno de esto, es que lo podemos graficar:
Plot_Comparacion = False
if Plot_Comparacion:
    (Solv_i_Cooper_Func - Solv_i_Zeb_Func).plot()
    
