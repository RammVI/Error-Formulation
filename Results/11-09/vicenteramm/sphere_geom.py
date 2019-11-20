
# 11-07-19

'''
Builds spherical boundary from decahedron
'''

import numpy as np
import os

def decahedron_points(r):
    '''
    Returns vert_array and face_array for a given radii (r) in armstrong.
    Note: Do not move the center of the sphere! This may cause troubles on smoothing function.
    '''
    
    # from http://www.sacred-geometry.es/?q=es/content/phi-en-los-s%C3%B3lidos-sagrados
    
    rho = (1. + np.sqrt(5.)) * 0.5
    
    if r==0:
        return None
    
    scale = (r) /np.sqrt(1.+rho**2.)
    
    vert_array = np.zeros((12,3))
    
    vert_array[0]  = np.array((   0. , -rho  ,   1. ))
    vert_array[1]  = np.array((   0. ,  rho  ,   1. ))
    vert_array[2]  = np.array((   0. ,  rho  ,  -1. ))
    vert_array[3]  = np.array((   0. , -rho  ,  -1. ))
    vert_array[4]  = np.array((   1. ,  0.   ,  rho ))
    vert_array[5]  = np.array((  -1. ,  0.   ,  rho ))
    vert_array[6]  = np.array((  -1. ,  0.   , -rho ))
    vert_array[7]  = np.array((   1. ,  0.   , -rho ))
    vert_array[8]  = np.array((  rho ,  1.   ,   0. ))
    vert_array[9]  = np.array(( -rho ,  1.   ,   0. ))
    vert_array[10] = np.array(( -rho , -1.   ,   0. ))
    vert_array[11] = np.array((  rho , -1.   ,   0. ))
    
    face_array = np.zeros((20,3))
    
    face_array[0]  = np.array(( 11 ,  7 ,  8 )) 
    face_array[1]  = np.array((  8 ,  7 ,  2 )) 
    face_array[2]  = np.array((  8 ,  2 ,  1 )) 
    face_array[3]  = np.array((  1 ,  2 ,  9 )) 
    face_array[4]  = np.array((  1 ,  9 ,  5 )) 
    face_array[5]  = np.array((  1 ,  5 ,  4 )) 
    face_array[6]  = np.array((  5 ,  0 ,  4 )) 
    face_array[7]  = np.array((  4 ,  0 , 11 )) 
    face_array[8]  = np.array((  7 ,  3 ,  6 )) 
    face_array[9]  = np.array(( 11 ,  8 ,  4 )) 
    face_array[10] = np.array((  8 ,  1 ,  4 ))
    face_array[11] = np.array((  6 , 10 ,  9 ))
    face_array[12] = np.array((  2 ,  6 ,  9 ))
    face_array[13] = np.array(( 10 ,  3 ,  0 ))
    face_array[14] = np.array(( 10 ,  0 ,  5 ))
    face_array[15] = np.array((  5 ,  9 , 10 ))
    face_array[16] = np.array((  6 ,  3 , 10 ))
    face_array[17] = np.array((  3 , 11 ,  0 ))
    face_array[18] = np.array((  7 , 11 ,  3 ))
    face_array[19] = np.array((  7 ,  6 ,  2 ))
        
    return vert_array*scale , face_array.astype(int)

def smothing_func( v_i , r ):
    '''
    Softens the boundary to achieve spherical behaivour.
    v_i : vertex to move to the boundary
    r   : sphere radii
    '''
    t = r /np.linalg.norm(v_i)
    return v_i * t

def pqr_assembly( x_q , q , mol_name ):
    '''
    Writtes {mol_name}.pqr file for a given position (x_q) and charges (q).
    '''
    
    pqr_text = open( os.path.join( 'Molecule' , mol_name , mol_name +'.pqr' ) , 'w+' )
    
    # .pqr format is:
    # recordName serial atomName residueName chainID residueNumber X Y Z charge radius
    
    if np.shape(x_q) == (3,):
        position = ' '.join( x_q.astype(str) )
        
        radii = 1. # Is implicitly defined on sphere radii
        
        text = 'ION      {0} None {1}   1      {2}  {3}  {4:5f} \n'.format(
        1 , mol_name , position , q , radii)
        
        pqr_text.write(text)
        pqr_text.close()
    
    else:
        count = 0
        for x in x_q:
            position = ' '.join( x.astype(str) )

            radii = 1. # Is implicitly defined on sphere radii

            text = 'ION      {0} None {1}   1      {2}  {3}  {4:5f} \n'.format(
            count+1 , mol_name , position , q[count] , radii)

            pqr_text.write(text)

            count+=1

        pqr_text.close()
    return None

def suffix_names(N_s):
    instances = np.empty( (N_s+1 , )).astype(str)
    for i in range(N_s +1 ):
        instances[i] = '-s{0:d}'.format(i)
    return instances

def moving_files(name , percentaje , s):


    path      = os.path.join( 'Molecule' , name)
    dest_path = os.path.join( path , str( int(percentaje*100)) , str(s) )
    os.system('mv {0}*.vtk {0}*.txt {0}*.face {0}*.vert {0}*.msh {0}*.off {1}'.format(path+'/' , dest_path) )
    
    return