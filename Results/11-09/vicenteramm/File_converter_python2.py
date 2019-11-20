
import numpy as np
import os

def face_and_vert_to_off(face_array , vert_array , path , file_name ):
    '''
    Creates off file from face and vert arrays.
    '''
    off_file = open( os.path.join( path , file_name + '.off' ) , 'w+')
    
    off_file.write( 'OFF\n' )
    off_file.write( '{0} {1} 0 \n'.format(len(vert_array) , len(face_array)) )
    for vert in vert_array:
        off_file.write( str(vert)[1:-1] +'\n' )
    
    for face in face_array:
        off_file.write( '3 ' + str(face - 1)[1:-1] +'\n' )
        
    off_file.close()
    
    return None

def Improve_Mesh(face_array , vert_array , path , file_name ):
    '''
    Executes ImproveSurfMesh and substitutes files.
    '''
    
    #os.system('export LD_LIBRARY_PATH=/vicenteramm/lib/')
    face_and_vert_to_off(face_array , vert_array , path , file_name)
    
    Improve_surf_Path = '/home/vicenteramm/Software/fetk/gamer/tools/ImproveSurfMesh/ImproveSurfMesh'
    os.system( Improve_surf_Path + ' --smooth --correct-normals ' + os.path.join(path , file_name +'.off')  )
    
    os.system('mv  {0}/{1}'.format(path, file_name + '_improved_0.off ') + 
                                 '{0}/{1}'.format(path, file_name + '.off '))
    
    new_off_file =  open( os.path.join( path , file_name + '.off' ) , 'r').read().split('\n')
    #print(new_off_file)
    
    num_verts = int(new_off_file[1].split()[0])
    num_faces = int(new_off_file[1].split()[1])

    new_vert_array = np.empty((0,3))
    for line in new_off_file[2:num_verts+2]:
        new_vert_array = np.vstack((new_vert_array , line.split() ))


    new_face_array = np.empty((0,3))
    for line in new_off_file[num_verts+2:-1]:
        new_face_array = np.vstack((new_face_array , line.split()[1:] ))


    new_vert_array = new_vert_array.astype(float)
    new_face_array = new_face_array.astype(int  ) + 1
    
    return new_face_array , new_vert_array