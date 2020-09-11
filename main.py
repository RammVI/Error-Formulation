
# All the .py files resumed in 1 python shell
# August 8th 2020

import bempp.api, numpy as np, os, time, argparse , getpass
from math import pi

from datetime import datetime
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

username = getpass.getuser()
os.system('export LD_LIBRARY_PATH="/home/{0}/lib"'.format(username))

args   = parse_aux_func()
fine_grid_maker(args.Molecule , dens_f=40.0)
suffix = suffixes(args.N_it)

create_txt_file(args.Results_name)

for i in range(args.N_it):
    ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , tot_solving_time , S_trad_time , S_Ap_time , S_Ex_time ,
                op_time_U , assem_time_U , solv_time_U , it_count_U ) = main( 
             args.Molecule , args.Density , suffix[i] , suffix[i+1] , args.Percentaje ,  N=25 ,  N_ref=args.N_ref ,
             smooth = True , refine = True, Use_GAMer = True        , sphere=args.sphere      , Estimator = args.E,
             x_q = None, q=None , r = np.nan)
    
    text = '{0} & {1:f} & {2:d}'.format( args.Molecule, args.Density , args.N_ref  ) + ' & {0:d} & {1:d} & {2:.10f} & {3:.10f} & {4:.10f} & {5:.10f} & {6:.10f} & {7:.10f} & {8:.10f} & {9:.10f} & {10:d} \n'.format( 
                     N_elements , N_El_adj  , S_trad , S_Ap , S_Ex , tot_solving_time , S_trad_time , S_Ap_time , S_Ex_time , op_time_U , it_count_U)
    
    append_txt_file(args.Results_name , text)
