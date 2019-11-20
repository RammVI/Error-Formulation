
# Built in python3.6
# 25.07.19

import meshio
import os

path = os.path.join('Molecule' , 'sphere_cent' , '30' , 'True')

file_name = 'sphere_cent_0-s5'

def create_off_file(file_name , path):

    mesh = meshio.read( os.path.join(path , file_name + '.msh') )

    meshio.write( os.path.join( file_name + '.off' ) , mesh , file_format='off')

    #os.rename( file_name + '.off' , os.path.join( path , file_name + '.off' ) )
    
    off_text = open( file_name + '.off' , 'r')
    
    aux_text = off_text.readlines()[3:] 
    
    new_text = open( os.path.join( file_name + '_aux.off') , 'w+')
    
    new_text.write('OFF\n')
    for line in aux_text:
        
        if line == '\n': continue
            
        new_text.write(line)
    
    new_text.close()
        
    #file_name = '1CID.pdb.output.surf'
        
    improve_mesh_loc = os.path.join( 'Software','gamer','tools','ImproveSurfMesh', 'ImproveSurfMesh' )
    
    exe = str.join( ' ', ('./'+ improve_mesh_loc , '--smooth', file_name + '_aux.off' )) 
    
    print('Executing: ' + exe)
    print('./Software/gamer/tools/ImproveSurfMesh/ImproveSurfMesh --smooth  --correct-normals {0}_aux.off'.format(file_name))
    #print('./Software/gamer/tools/ImproveSurfMesh/ImproveSurfMesh --smooth sphere_cent_0-s5_aux.off ')
    print(exe)
    
    os.system( exe )
    
    #os.system('rm {0}_aux.off'.format(file_name) )
    #os.system('rm {0}.off'.format(file_name) )
    
    #print('mv {0} {1}/.'.format(file_name+'_aux_improved_0.off' , path))
    
    
    #os.system('sudo mv {0} {1}/.'.format(file_name+'_aux_improved_0.off' , path))
    #mesh = meshio.read( os.path.join(path , file_name + '_aux_improved_0.off') )    
    
    vert_array = mesh.points
    face_array = mesh.cells['triangle']

    vert_txt = open( os.path.join( path , file_name + '.vert' ) , 'w' )

    for vert in vert_array:
        vert_txt.write(str(vert)[1:-1] + '\n')

    vert_txt.close()

    face_txt = open( os.path.join( path , file_name + '.face' ) , 'w' )

    for face in face_array:
        face_txt.write(str(face)[1:-1] + '\n')

    face_txt.close()
    
    return None

create_off_file(file_name , path)