
import meshio
import os

path = os.path.join('Molecule' , 'sphere_cent' , '30' , 'True')

file_name = 'sphere_cent_0-s5'

def create_off_file(file_name , path):
    
    total_path = 'home/vicente/Documentos/BEM'

    mesh = meshio.read( os.path.join( total_path , path , file_name + '.msh') )

    meshio.write( os.path.join( total_path , file_name + '.off' ) , mesh , file_format='off')

    #os.rename( file_name + '.off' , os.path.join( path , file_name + '.off' ) )
    
    off_text = open( total_path +'/' + file_name + '.off' , 'r')
    
    aux_text = off_text.readlines()[3:] 
    
    new_text = open( os.path.join( total_path , file_name + '_aux.off') , 'w+')
    
    new_text.write('OFF\n')
    for line in aux_text:
        
        if line == '\n': continue
            
        new_text.write(line)
    
    new_text.close()
    
    improve_mesh_loc = os.path.join( total_path , 'Software','gamer','tools','ImproveSurfMesh', 'ImproveSurfMesh' )
    
    exe = str.join( ' ', ('./' +  improve_mesh_loc , '--smooth', total_path + '/' + file_name + '_aux.off' )) 
    
    print('./Software/gamer/tools/ImproveSurfMesh/ImproveSurfMesh --stats sphere_cent_0-s5_aux.off')
    print(exe)
        
    os.system( exe )
    
    return None

create_off_file(file_name , path)