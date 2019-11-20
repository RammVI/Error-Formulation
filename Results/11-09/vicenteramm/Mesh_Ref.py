
import bempp.api, numpy as np
from math import pi
import os
from Grid_Maker_R2 import *


def text_to_list(mol_name , total_suffix , txt_format , info_type=float):
    '''
    Rutine which builds lists from text files that contain the vert and face info
    mol_name : Abreviated name of the molecule
    total_suffix   : text added after molecule name, for an easier handling of files 
                     may be taken as a total and used like _{density}-{it_count}
    txt_format : file format
    info_type  : float or int
    '''
    path = os.path.join('Molecule',mol_name)
    
    list_txt = open( os.path.join(path , mol_name +total_suffix + txt_format) ).read().split('\n')

    listing = np.empty((0,3))
    for line in list_txt[:-1]:
        info    = np.array(line.split()[:3])
        listing = np.vstack( ( listing , info ) ).astype(info_type)
    return listing
     

    
def add_vertex_in_face_center( face , face_array , vert_array ):
    '''
    This function adds a point at the certer of a given face, the added vertex is at the end of the vert_array
    listing.
    face : np.array((0,3)) containing the position of 3 verts in vert_array
    face_array : np.array() which has the verts position of all the faces
    vert_array : np.array() containing sort vertex.
    '''
    # Let's load the 3 vertex and create the centered one

    v1 , v2 , v3 = vert_array[face[0]-1] , vert_array[face[1]-1] , vert_array[face[2]-1]
    v_centered   = (v1+v2+v3)/3.
    
    # And adding the vert to the vert_array list
    new_vert_array = np.vstack((vert_array , v_centered))
    
    #print(v1,v2,v3)
    
    # Let's save the position of the last added vert
    Vert_pos = len(new_vert_array)
    
    # First, the face data must be deleted from face_array
    new_face_array = np.empty((0,3))
    for data in face_array:
        #print(data,face,(data!=face).all())
        if not np.array_equal( data , face):
            new_face_array = np.vstack((new_face_array , data))
        
    # Then, we create the 3 possible combinations
    f1 = np.array((face[0],face[1],Vert_pos))
    f2 = np.array((face[1],face[2],Vert_pos))
    f3 = np.array((face[0],face[2],Vert_pos))
    
       
    # Finally, the 3 new faces are added to the face_array list
    for f in [f1,f2,f3]:
        new_face_array = np.vstack((new_face_array,f))
    
    return v_centered , new_face_array.astype(int) , new_vert_array

def mesh_flat_refinement(mol_name,faces,face_array,vert_array , suffix , dens):
    '''
    This function builds the msh file for the specified faces.
    faces : Faces arrays to be remeshed
    face_array : array containing the face vertex
    vert_array : array containing the vertex points
    '''
    for f in faces:
        new_stuff  = add_vertex_in_face_center( f , face_array , vert_array )
        vert_array , face_array  = new_stuff[2],new_stuff[1]


    vert_and_face_arrays_to_text_and_mesh(mol_name , vert_array , face_array , suffix , dens)
    print('Remesh Ready')
    return None

def random_face_list(Number_of_faces , face_array):
    '''
    Shall be used as a debug function !
    '''
    if Number_of_faces>len(face_array):
        print('Error: Numer of faces is greater than length of face_array!')
    random_faces = np.random.randint(1,len(face_array) , size =Number_of_faces)

    refined_faces = np.empty((0,3))
    for i in random_faces:
        refined_faces = np.vstack( (refined_faces , face_array[i]))
        
    return refined_faces

def mesh_ref(mol_name , refined_faces , input_file_suffix , output_file_suffixx , starting_density):
    '''
    Global function that refines the mesh.
    mol_name : Abreviated name of the molecule
    refined_faces : np.array((N,3)) which contains faces with respective verts positions to be refined.
    input_file_suffix : Suffix of the starting file
    output_file_suffix : Suffix of the outputfile
    starting_density   : Starting density (For an easy handling of files)
    '''
    mol_name = 'methanol'
    path = os.path.join('Molecule',mol_name)

    starting_density = 2.0

    suffix   = ''

    # Information is imported
    aux_path = '_'+str(starting_density)+input_file_suffix
    face_array = text_to_list(mol_name , aux_path , '.face' , info_type=int  )
    vert_array = text_to_list(mol_name , aux_path , '.vert' , info_type=float)


    #refined_faces = random_face_list(50 , face_array) 

    mesh_flat_refinement(mol_name , refined_faces,face_array,vert_array, output_file_suffixx , dens = starting_density)
    
    return None

def splitting_by_half_triangle(face , vert_array , face_array):
    
    
    v1 , v2 , v3 = vert_array[face[0]-1] , vert_array[face[1]-1] , vert_array[face[2]-1]

    v1v2 = v2 - v1
    v3v2 = v3 - v2
    v3v1 = v3 - v1
    
    mayor = 1
    
    return mayor
    
    