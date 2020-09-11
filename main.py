
# All the .py files resumed in 1 python shell
# August 8th 2020

import bempp.api, numpy as np, os, time, argparse
from math import pi

from datetime import datetime
from constants import values
from constants import mesh_info

from bempp.api.operators.potential import laplace as lp
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz

from Grid_Maker_R2, Mesh_Ref_V2, quadrature, Potential_Solver, sphere_geom, File_converter_python2, analytical import *

parser = argparse.ArgumentParser()

parser.add_argument("-Molecule" , "-M"   , help="Molecule name in Molecule/mol_name/.")
parser.add_argument("-Density"  , "-D"   , help="Starting mesh density", type=float)
parser.add_argument("-N_it"     , "-n"   , help="Number of iterations/consecutive mesh refinements",type=int)
parser.add_argument("-E"        , help="Error estimator. Can be E_u or E_phi",type=str)

#Optional input here
parser.add_argument("-N_ref"      , help="[Optional] Number of elements for the phi mesh. Default=0" ,
                    default=0,type=int)
parser.add_argument("-Percentaje" , "-P" , help="[Optional] Percentaje of the elements to be refined, acording to their error contribution", type=float)
parser.add_argument("-sphere"            , help="[Optional] True if the geometry is a sphere, and false if not" ,
                                    action="store_true", default=False)
parser.add_argument("-Results_name"      , help="[Optional] Result text file name in Results/." ,
                                    default="Results_{0}_{1}_{2}.txt".format(datetime.now().year ,
                                    datetime.now().month, datetime.now().day), type=str)

args = parser.parse_args()
suffix = suffixes(args.N_it)

create_txt_file(args.Results_name)

for i in range(N_it):
    ( S_trad , S_Ap , S_Ex , N_elements , N_El_adj , tot_solving_time , S_trad_time , S_Ap_time , S_Ex_time ,
                op_time_U , assem_time_U , solv_time_U , it_count_U ) = main( 
             args.Molecule , args.Density , suffix[i] , suffix[i+1] , args.Percentaje ,  N=25 ,  N_ref=args.N_ref ,
             smooth = True , refine = True, Use_GAMer = True        , sphere=args.sphere      , Estimator = args.E,
             x_q = None, q=None , r = np.nan)
    
    #sphere case should be defined with a .pqr distribution? Should I add a couple of example tests?
    
    text = '{0} & {1:f} & {2:d}'.format( molecule, dens , args.N_ref  ) +
           ' & {0:d} & {1:d} & {2:.10f} & {3:.10f} & {4:.10f} & {5:.10f} & {6:.10f} & {7:.10f} & {8:.10f} & {9:.10f} & {10:d} \n'.format( 
                     N_elements , N_El_adj  , S_trad , S_Ap , S_Ex , tot_solving_time , S_trad_time , S_Ap_time , S_Ex_time , op_time_U , it_count_U)
    
    append_txt_file(args.Results_name , text)
