
# All the .py files resumed in 1 python shell
# August 8th 2020

import bempp.api, numpy as np, os, time
from math import pi

from constants import values
from constants import mesh_info

from bempp.api.operators.potential import laplace as lp
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

from Grid_Maker_R2, Mesh_Ref_V2, quadrature, Potential_Solver, sphere_geom, File_converter_python2, analytical import *

txt_name = 'Results.txt'
percentaje , dens = 0.1 , 1.0
suffixes = ['-0','-1','-2','-3','-4','-5','-6','-7']  # This will refine up to '-6'.

# For the sphere case uncomment the following lines
#r = 1.
#x_q = np.array( [[  1.E-12 ,  1.E-12 ,   r/2. ]]  )
#q = np.array( [1.] )

for i in range(len(suffixes[:-1])):
    ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , tot_solving_time , S_trad_time , S_Ap_time , S_Ex_time ,
                op_time_U , assem_time_U , solv_time_U , it_count_U ) = main( 
             'methanol' , dens , suffix[i] , suffix[i+1] , percentaje ,  N=25 ,  N_ref=0 ,
             smooth = False , refine = True, Use_GAMer = False , sphere=False , Estimator = 'E_u',
             x_q = None, q=None , r = np.nan)
