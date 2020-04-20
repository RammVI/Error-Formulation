
# Some rules about the mesh refinement:
# 1. A percentaje of the faces must be refined and 4 new triangles are born
#    where 3 new points are added in each edge center, and the input data is the residual
#        .                  ___ .____
#       /\                 |\  /\  /|
#      /  \      Results   | \/__\/ |
#     /    \       in      | /\  /\ |
#    /______\              |/__\/__\|

#     in the solvation energy.
# 2. Adjacent triangles are split but half UNLESS:
#    2.1 They are adjacent to 2 triangles to be refinated into 4 new triangles
# Do this until there is no more triangles in 2.1 .

# Also, the possibility to extrapolate the point to the real boundary will be
# set in a function named new_point()

# 11-09 Adjoint mesh usefull

#import bempp.api
import numpy as np
from math import pi
import os

def search_unique_position_in_array(array , main_array):
    '''
    -
    '''
    
    c=0
    array_is_contained = False
    
    for sub_array in main_array:
        
        sub_array_counter = 0
        
        equality = True
        for i in sub_array:
            
            if not isclose( array[sub_array_counter] , sub_array[sub_array_counter]  ):
                equality = False
                continue
                
            sub_array_counter+=1
            
        if equality:
            array_is_contained = True
            break
        c+=1
        
    if array_is_contained:
        return c
    
    else: 
        return -1
    
    
    
def isclose(a, b, rel_tol=1e-5, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
    
def search_multiple_positions_in_array( arrays , main_array ):
    '''
    Returns the position of each array (contained in arrays) in main_array.
    Arrays must save the information in rows, for example
    arrays = np.array((f1x , f1y , f1z),
                       f2x , f2y , f2z)....)
    '''    
    
    positions = (-1)*np.ones( (len(arrays),1) )
    
    i=0
    for array in arrays:
        
        c = 0
        
        for sub_array in main_array:
            
            if(array == sub_array).all():
                positions[i] = c                
            
            c+=1
        i+=1
    return positions

def common_verts_between_2_triangles( face_1 , face_2 ):
    '''
    Returns the 2 common positions in face_1 and face_2
    '''
    common_verts = np.zeros( (2, ) )
    c = 0    
    for i in face_1:
        for j in face_2:
            if i==j and i not in common_verts:
                common_verts[c] = i
                c+=1
                
    return common_verts.astype(int)
        
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


def value_assignor_starter(face_array , soln , percentaje):
    '''
    Assigns face's value to the number of new triangles born (Desirable is 0, 2 or 4).
    May use first!
    '''
    
    refined_faces = np.zeros( (len(soln), ) )
    
    if   type(percentaje) == int:
        
        listing    = np.zeros( (2 , len(soln)) )
        listing[0] = soln
        listing[1] = np.arange(len(soln))
        
        listing = zip(*listing)
        listing.sort(key=lambda x: x[0])
        
        listing = np.array(listing)
        
        counter = 0
        for face in listing[:,1][::-1]:

            refined_faces[int(face)] = 4
            
            if counter+1 == percentaje:
                break
            counter+=1
        
        return refined_faces
    
    elif type(percentaje) == float:
        
        percentaje_sum = 0.0
        
        listing    = np.zeros( (2 , len(soln)) )
        listing[0] = soln / np.sum( soln )  
        listing[1] = np.arange(len(soln))
        
        listing = zip(*listing)
        listing.sort(key=lambda x: x[0])
        
        listing = np.array(listing)[::-1]
        
        counter = 0
        for face in listing[:,1]:
            
            refined_faces[int(face)] = 4
            percentaje_sum += listing[ counter , 0 ]
            
            if percentaje_sum > percentaje:
                break
            
            counter +=1
        
        return refined_faces
    
    
    
    
    
        # This was used for past versions and refines a percentaje of the faces,
        # this was changed to the faces contributing to the total percentaje
        
        # This extract the class value
        #separator = np.sort(soln)[int(percentaje * len(soln) )]
        #refined_faces = np.zeros( (len(soln), ) )

        #c = 0
        #for v in soln:

        #    if v>separator:
        #        refined_faces[c] = 4
        #    c+=1

        #return refined_faces

    else:
        print('Fatal error')

        return None
    
    
def coincidence_between_1D_arrays(array1 , array2 , coinc_num ):
    '''
    Search for coincidences between 2 arrays, returns True if arrays have coinc_num of coincidences
    and also returns the coincidences.
    '''
    coincidences = np.zeros( (coinc_num , ) ) 
    
    c=0
    for a1 in array1:
        
        for a2 in array2:
            
            if a1 == a2:
                
                coincidences[c] = a1
                c+=1
                
    if 0 in coincidences:
        return False , np.zeros( (coinc_num , ) )
                
    return True , coincidences

def adjacent_faces(face_num , face_array , return_face_number):
    '''
    Searchs for adjacent faces to be refined and sets value 2 if the face is adjacent to only 1
    face which has 4 new triangles to be refined or 3 if this points has a 2 adj. faces to be refined, and 
    4 if the face has three adj triangles to be refined.
    face_num    : Position of the face in face_array
    face_status : Face position in face_array to be refined into 4 triangles
    return_face_number : Boolean True  if it is required to return the position of the face in face_array.
    '''
    
    adj_faces = np.zeros( (len(face_array), ) )
    
    pointed_face = face_array[face_num]
        
    adj = True
        
    T1 , T2 , T3 = -1 , -1 , -1
        
    cj = 0
        
    for face in face_array:
                        
        if (face == pointed_face).all():
            cj+=1
            continue
            
        Boolean , coincidences = coincidence_between_1D_arrays(pointed_face , face , coinc_num=2 )
        
        if Boolean and T1 == -1:
            T1 = cj
            cj+=1
            continue
        if Boolean and T2 == -1:
            T2 = cj
            cj+=1
            continue
        if Boolean and T3 == -1:
            T3 = cj
            cj+=1
            continue
        
        #if f1 in pointed_face:
                
        #    if (f2 in pointed_face or f3 in pointed_face) and T1 == -1:
        #        T1 = cj
        #        cj+= 1
        #        continue
    
        #if f2 in pointed_face:
                
        #    if (f1 in pointed_face or f3 in pointed_face) and T2 == -1:
        #        T2 = cj
        #        cj+= 1
        #        continue
                
        #if f3 in pointed_face:
        #    print('third cond ', face, pointed_face)
        #    if f1 in pointed_face or f2 in pointed_face:
        #        print('T3: ',face)
        #        T3 = cj
        #        continue
        cj+=1
    
    adj_faces[T1] , adj_faces[T2] , adj_faces[T3] = 2 , 2 , 2 
    
    
    if return_face_number:
        return T1,T2,T3
    
    else:
        return adj_faces

def adj_assigner_value( face_array , soln , percentaje):
    '''
    Gives each face a value depending of its solution value
    face_array : Array containing vertex position on a array of vertices
    soln       : Quantity/Face atribute to be the criteria to refinate
    percetaje  : Percentaje of the mainly faces contributing to error (if float) or
                 number of faces to be refined
    '''
    
    first_status = value_assignor_starter(face_array , soln , percentaje)
    
    adj_status   = np.zeros( (len(face_array), ) )
    
    face_num = 0    # face counter
    
    for value in first_status:
        
        if value == 4:
            
            adj1 , adj2 , adj3 = adjacent_faces( face_num , face_array , return_face_number = True)
            
            if adj_status[adj1] != 4:
                adj_status[adj1] +=1
                
            if adj_status[adj2] != 4:
                adj_status[adj2] +=1
                
            if adj_status[adj3] != 4:
                adj_status[adj3] +=1
            
            
        face_num+=1
        
    return adj_status 
        
def final_status(face_array , soln , percentaje ):
    '''
    Runs into status until there is no more triangles splitted in 3 (rule 2.1)
    '''
    
    status = value_assignor_starter(face_array , soln , percentaje )+  \
                adj_assigner_value(face_array, soln , percentaje)
    
    face_num = 0
    aux_status = status.copy()
    for s in status:
        if s >4:
            aux_status[face_num] = 4
        face_num+=1
        
    status = aux_status.copy().astype(int)        
    
    iteration_restrictor = 0
    
    Calculating = True
    
    while 2 in status and iteration_restrictor<20 or Calculating:
        
        if iteration_restrictor == 15:
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            print('Adapting the mesh is taking too long! - Breaking ')
            
        # Changing the 2 values for 4
        
        face_num = 0
        aux_status = status.copy()
        for s in status:
            if s == 2 or s == 3 :
                aux_status[face_num] = 4
            face_num += 1
        
        Calculating = True
        
        status = aux_status.copy()
        
        adj_status   = np.zeros( (len(face_array), ) )
        status_4     = np.zeros( (len(face_array), ) )
        
        face_num = 0
        for value in status:
        
            if value == 4:

                adj1 , adj2 , adj3 = adjacent_faces( face_num , face_array , return_face_number = True)

                if status[adj1] != 4:
                    adj_status[adj1] +=1

                if status[adj2] != 4:
                    adj_status[adj2] +=1

                if status[adj3] != 4:
                    adj_status[adj3] +=1
                    
                status_4[face_num] = 4
            face_num+=1

        status = adj_status + status_4  
            
        face_num = 0
        aux_status = status.copy()
        for s in status:
            if s>4:
                aux_status[face_num] = 4
                
            face_num+=1
        
        status = aux_status.copy()
        
        iteration_restrictor += 1
        
        Calculating = False
        
    return aux_status

def mesh_refiner(face_array , vert_array , soln , percentaje ):
    '''
    Refines the mesh
    '''
    
    status = final_status(face_array , soln , percentaje )
    
    new_faces_array = np.empty((0,3))
    new_vert_array = vert_array.copy()
    
    aux_face_array = face_array.copy()
    
    face_counter = 0
    
    for face in face_array:
        
        
        if status[face_counter] == 1:
             
            adj_1_pos , adj_2_pos , adj_3_pos = adjacent_faces( face_counter , face_array , True )
            
            #Search for the adjacent face that has status 4
            
            if status[adj_1_pos] == 4:
                ady_face = face_array[adj_1_pos]
                
            elif status[adj_2_pos] == 4:
                ady_face = face_array[adj_2_pos]
                
            elif status[adj_3_pos] == 4:
                ady_face = face_array[adj_3_pos]
                
            # Now search for the 2 common vertex
            v1_pos , v2_pos = common_verts_between_2_triangles( ady_face , face ) - 1 
            #v1_pos and v2_pos have an absolute value - starts from 0
            
            
            
            v1 , v2 = vert_array[v1_pos] , vert_array[v2_pos]
            
            for v3_pos in face-1:
                if v3_pos!=v1_pos and v3_pos!=v2_pos:
                    break      
            
            new_vert = newvert( v1 , v2 )
            
            # Now let's check if this new vert is already in new_vert_array
            test_position = search_unique_position_in_array(new_vert , new_vert_array )
            
            if test_position == -1:
                new_vert_array = np.vstack( (new_vert_array , new_vert ) )
                test_position = len(new_vert_array)-1  #test_position is also absolute
                
                
            #Because v1_pos and v2_pos are random-like positions, let's sort it 
            # using a case by case test
            f1 , f2 , f3 = face - 1    
            
            if   (f1 == v1_pos and f2 == v2_pos) or (f1 == v2_pos and f2 == v1_pos):
                
                new_face_1 = np.array( ( f1            , test_position , f3  ) ) +1
                new_face_2 = np.array( ( test_position , f2            , f3  ) ) +1
                
            elif (f2 == v1_pos and f3 == v2_pos) or (f2 == v2_pos and f3 == v1_pos):
                
                new_face_1 = np.array( ( f2            , test_position , f1 ) ) +1
                new_face_2 = np.array( ( f1            , test_position , f3 ) ) +1 
                
            elif (f3 == v1_pos and f1 == v2_pos) or (f3 == v2_pos and f1 == v1_pos):
                 
                new_face_1 = np.array( ( f3            , test_position , f2 ) ) +1
                new_face_2 = np.array( ( f2            , test_position , f1 ) ) +1
            
            new_faces_array = np.vstack((new_faces_array , new_face_1 ))
            new_faces_array = np.vstack((new_faces_array , new_face_2 ))
            
            # Also we have to delete the face, so let's assign the value (0,0,0) the deleted face
            aux_face_array[face_counter] = np.array((0,0,0))
                
        if status[face_counter] == 4:
            
            f1 , f2 , f3 = face_array[face_counter] - 1  # Absolute position in vert_array
            v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]
            
            new_vert_12 = newvert( v1 , v2 )
            new_vert_13 = newvert( v1 , v3 )
            new_vert_23 = newvert( v2 , v3 )
            
            # Let's check like in status == 1 if the new_vert_ij is in new_vert_array:
            
            pos = np.array((-1,-1,-1))
            c = 0
            for vert in (new_vert_12 , new_vert_13 , new_vert_23):
                test_position = search_unique_position_in_array( vert , new_vert_array )
                if test_position == -1:
                    new_vert_array = np.vstack( (new_vert_array , vert ) )
                    test_position = len(new_vert_array) - 1 # Absolute value!
                    
                pos[c] = test_position
                c+=1
            
            v12_pos , v13_pos , v23_pos = pos
            
            new_face_1 = np.array( ( f1     , v12_pos , v13_pos ) ) + 1 
            new_face_2 = np.array( ( v12_pos, v23_pos , v13_pos ) ) + 1
            new_face_3 = np.array( ( v12_pos, f2      , v23_pos ) ) + 1
            new_face_4 = np.array( ( v13_pos, v23_pos , f3      ) ) + 1
            
            new_faces_array = np.vstack( ( new_faces_array , new_face_1 ) )
            new_faces_array = np.vstack( ( new_faces_array , new_face_2 ) )
            new_faces_array = np.vstack( ( new_faces_array , new_face_3 ) )
            new_faces_array = np.vstack( ( new_faces_array , new_face_4 ) )
            
            # Also we have to delete the face, so let's assign the value (0,0,0) the deleted face
            aux_face_array[face_counter] = np.array((0,0,0))
            
        face_counter+=1
    
    final_face_array = np.empty((0,3))
    
    for face in aux_face_array:
        if (face.astype(int) == np.array( (0 , 0 , 0 ) ) ).all():
            continue
        final_face_array = np.vstack( (final_face_array , face ) )
        
    for new_face in new_faces_array:
        final_face_array = np.vstack( (final_face_array , new_face ) )
        
    return final_face_array.astype(int) , new_vert_array


def newvert(vA,vB):
    return 0.5*(vA+vB)

def smoothing_vertex( vert_array , fine_vert_array ):
    '''
    Smooths verts from a finer mesh
    '''
    
    smooted_vert_array = np.empty((0,3))
    
    for vert in vert_array:
        
        dist = 1000.
        
        for f_vert in fine_vert_array:
            
            #Searchs for closest vertex.
            radii = np.linalg.norm( vert - f_vert )
            
            if radii<dist:
                
                dist = radii
                added_vert = f_vert
                
        smooted_vert_array = np.vstack( (smooted_vert_array , added_vert ) )
    
    return smooted_vert_array 
 
    
def is_interior_triangle(adj_vertices , local_vert_array):
    '''
    Returns true if the designed triangle from the adjoint mesh is totally inside the normal triangle.
    adj_vertices : Three vertices from the same elements
    local_vert_array   : Array of vertices for a given triangle.
    '''
    
    v1a , v2a, v3a = adj_vertices
    v1 , v2 , v3 = local_vert_array
    
    v12 , v13 , v23 = 0.5*(v1+v2) , 0.5*(v1+v3) , 0.5*(v2+v3)
    
    c=0
    for adj_vertex in adj_vertices:
        for vert in (v1, v2 , v3 , v12 , v13 , v23 ):
            if isclose(adj_vertex , vert):
                print(adj_vertex , vert)
                c+=1
    if c == 3:
        return True
    else:
        return False
    
def elements_position_in_normal_grid(adj_face_array , adj_vert_array , face_array , vert_array ):
    '''
    Returns the element for which the adjoint element is contained.
    inputs:
    adj_face_array , adj_vert_array , face_array , vert_array
    '''
    c=0

    adj_relation = np.zeros((len(adj_face_array) , ))

    for element in face_array:
        f1 , f2 , f3 = element - 1
        v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]
        v12 , v13 , v23 = 0.5*(v1+v2) , 0.5*(v1+v3) , 0.5*(v2+v3)

        count_face = 0
        for adj_element in adj_face_array:
            f1a , f2a , f3a = adj_element - 1
            v1a , v2a , v3a = adj_vert_array[f1a] , adj_vert_array[f2a] , adj_vert_array[f3a]


            if (search_unique_position_in_array(v1a ,(v1 , v2 , v3 , v12 , v13 , v23) )!= -1 and
                search_unique_position_in_array(v2a ,(v1 , v2 , v3 , v12 , v13 , v23) )!= -1 and
                search_unique_position_in_array(v3a ,(v1 , v2 , v3 , v12 , v13 , v23) )!= -1):
                adj_relation[count_face] = c

            count_face +=1

        c+=1

    return adj_relation


def normals(face_array , vert_array):
    '''
    Returns an array (N,3) containing the normal vectors to each face.
    '''
    
    normals_array = np.empty((0,3))
    for f in face_array:
        f1 , f2 , f3 = f

        v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]
        
        normal = np.cross(v2-v1 , v3 - v1 )
    
        unit_n = normal/np.linalg.norm(normal)
        
        normals_array = np.vstack((normals_array , unit_n))
    
    return normals_array.astype(float)

def check_coplanar(triangle_1 , triangle_2 , n_1):
    
    v1_1 , v1_2 , v1_3 = triangle_1
    
    coplanar = True
    
    for point in triangle_2:
        
        condition = np.dot(n_1 , point - v1_1)
        
        if condition > 10.**-6:
            coplanar = False
            break
            
    return coplanar

def Area_of_a_Triangle(A,B,C):
    '''
    Returns the Area of a triangle for a trio of vertices.
    Params: 
    Vertex_1 , Vertex_2 , Vertex_3
    '''
    AB = B-A
    AC = C-A
    return 0.5 * np.linalg.norm( np.cross(AB,AC))

def check_contained_triangles_alternative_2(mesh_1 , mesh_2 , N_ref , assume_uniform = True , Tolerance=10**-10
                                           ,check_for_unnasigned_faces=True):
    '''
    Condicion de area y normales para revisar si el triangulo se encuentra dentro de la cara determinada.
    Params:
    mesh_1 : Coarse mesh
    mesh_2 : Finner mesh
    N_ref  : Number of refinements
    Assume_uniform : True o False
    Tolerance : Maximum observable difference between the area triangles and the normals.
    '''
    face_1 , vert_1 = np.transpose(mesh_1.leaf_view.elements) , np.transpose(mesh_1.leaf_view.vertices)
    face_2 , vert_2 = np.transpose(mesh_2.leaf_view.elements) , np.transpose(mesh_2.leaf_view.vertices)
    
    normal_1 , normal_2 = normals(face_1 , vert_1 ) , normals(face_2 , vert_2)
    
    relationship = np.zeros((len(face_2),)) - 1
    
    c1=0
    for f1 in face_1:
        
        v1_1 , v1_2 , v1_3 = vert_1[f1[0]] , vert_1[f1[1]] , vert_1[f1[2]]
        
        Area_1 = Area_of_a_Triangle(v1_1 , v1_2 , v1_3)
        
        n_1 = normal_1[c1]
        
        c2 = 0
        for f2 in face_2:
            if relationship[c2]<0:
                
                v2_1 , v2_2 , v2_3 = vert_2[f2[0]] , vert_2[f2[1]] , vert_2[f2[2]]
                
                n_2 = normal_2[c2]
                
                condition_counter = 0 
                for point in [v2_1 , v2_2 , v2_3]:
                    A1 , A2 , A3 = [Area_of_a_Triangle(v1_1,v1_2,point) , 
                                    Area_of_a_Triangle(v1_1,v1_3,point) ,
                                    Area_of_a_Triangle(v1_2,v1_3,point)   ]
                    
                    if np.abs((A1+A2+A3)-Area_1)<=Tolerance and np.linalg.norm(n_1 - n_2)<=Tolerance:
                        condition_counter +=1
                
                if condition_counter == 3:
                    
                    relationship[c2] = c1

                
            c2+=1
        c1+=1
        
    relationship = relationship.astype(int)
        
    if -1 in relationship:
        print('ERROR - DID NOT ASSIGN A FACE TO COARSE MESH')
        print(relationship)
        
    if check_for_unnasigned_faces:
        counter_check(relationship , mesh_1 , mesh_2 , N_ref )

    if assume_uniform:
        counter_check(relationship , mesh_1 , mesh_2 , N_ref)
    
    return relationship

def check_contained_triangles(mesh_1 , mesh_2 , N_ref , assume_uniform = True):
    
    face_1 , vert_1 = np.transpose(mesh_1.leaf_view.elements) , np.transpose(mesh_1.leaf_view.vertices)
    face_2 , vert_2 = np.transpose(mesh_2.leaf_view.elements) , np.transpose(mesh_2.leaf_view.vertices)
    
    normal_1 , normal_2 = normals(face_1 , vert_1 ) , normals(face_2 , vert_2)
    
    #print(normal_1 , normal_2)
    
    relationship = np.zeros((len(face_2),)) - 1
    
    c1=0
    for f1 in face_1:
        
        v1_1 , v1_2 , v1_3 = vert_1[f1[0]] , vert_1[f1[1]] , vert_1[f1[2]]
        
        n_1 = normal_1[c1]
        
        c2 = 0
        for f2 in face_2:
            if relationship[c2]<0:
                n_2 = normal_2[c2]
                v2_1 , v2_2 , v2_3 = vert_2[f2[0]] , vert_2[f2[1]] , vert_2[f2[2]]
                
                if np.linalg.norm(n_1-n_2)<10.**-6 or np.linalg.norm(n_1+n_2)<10.**-6 :
                    if check_coplanar([ v1_1 , v1_2 , v1_3 ] , [v2_1 , v2_2 , v2_3] , n_1 ):
                        relationship[c2] = c1
            c2+=1
        c1+=1
        
    relationship = relationship.astype(int)
        
    if -1 in relationship:
        print('ERROR - DID NOT ASSIGN A FACE TO COARSE MESH')
        print(relationship)
    
    
    
    if assume_uniform:
        counter_check(relationship , mesh_1 , mesh_2 , N_ref)
    
    #print(relationship)
    
    return relationship
                    
def counter_check(relationship , mesh_1 , mesh_2 , N_ref ):
    '''
    Counts the number of faces assigned to the coarse mesh, for a given number of unifrom refinement cycles.
    '''
    
    check_list = np.zeros(( len(np.transpose(mesh_1.leaf_view.elements)) , 2 )) 
    
    check_list[:,0] = np.arange(0, len(np.transpose(mesh_1.leaf_view.elements)) )
    
    for i in relationship:
        check_list[i,1]+=1
    
    if np.any(check_list[:,1].astype(int)!=4**N_ref):
        
        print('ERROR - TRIANGLE ASSIGNED TO OTHER TRIANGLE OF THE COARSE MESH')
        print(check_list)
        
    return None

def uniform_refinement_points(v1 , v2 , v3 ):
    '''
    Creates an array that contains ALL the points generated when doing 1 uniform refinement.
    '''
    Points = np.array([v1 , v2 , v3 , (v1+v2)/2. , (v2+v3)/2. , (v1+v3)/2.])
    N_face = np.array( [[0 , 3 , 5 ] , [3 , 1 , 4] , [4 , 2 , 5] , [3 , 4 , 5] ] )
    
    sorted_points = np.zeros([4 * 3 , 3])
    face_count = 0
    for face in N_face:
        for f in face:
            sorted_points[face_count] = Points[f]
            face_count +=1
            
    return sorted_points

def recursive_refinement_points(v1 , v2 , v3 , N_Ref , delete_repeated_points = False ):
    '''
    Creates an array that contains ALL the points generated when doing N_Ref uniform refinements.
    Params:
    v1 , v2 , v3 : Vertices
    N_Ref : Number of uniform refinements
    '''
    Ref = 0
    
    points = np.zeros((4**N_Ref*3 , 3))
    points[0] , points[1] , points[2] = v1 , v2 , v3
    
    while Ref<=N_Ref:
        
        aux_points = np.zeros((4**N_Ref*3 , 3))
        
        for c in range(4**(N_Ref-1)):
            v1 , v2 , v3 = points[3*c] , points[3*c+1] , points[3*c+2]
            
            if np.all(np.array([v1,v2,v3])<10**-8):
                break
        
            new_points = uniform_refinement_points( v1 , v2 , v3 )
            
            for aux_c in range(12):
                aux_points[12*c+aux_c] = new_points[aux_c]
            
        points = aux_points.copy()
        
        Ref+=1
    
    if delete_repeated_points:
        clean_set = np.empty((0,3))
        clean_set = np.vstack((clean_set , points[0] , points[1]))
           
        for p in points:
            
            is_contained = False
            for clean_p in clean_set:
                if np.linalg.norm(clean_p - p ) < 10**-8:
                    is_contained = True
            if not is_contained:
                clean_set = np.vstack([clean_set , p])
                
            
        points=clean_set.copy()
        
    return points*2

def similarities_counter(Small_Array , Big_Array , Tol = 10**-4):
    '''
    Returns True if all the Small_Array values are contained into the Big_Array for a given
    Tolerance
    '''
    c=0
    for small_point in Small_Array:
        for big_point in Big_Array:
            if np.linalg.norm(big_point-small_point)<Tol:
                c+=1
    
    if c==len(Small_Array)-1:
        return True
    return False

def unitary_refinement_generator(N_ref):
    '''
    Generates an array containing all the points inside the unitary triangle
    (0,0) , (1,0) , (0,1)
    Params:
    N_ref: Number of refinements
    '''
    Unit = np.array([ [0.,0.,0.] , [1. , 0. ,0.] , [0.,1., 0.] ])
    points = recursive_refinement_points(Unit[0],Unit[1],Unit[2], N_ref , delete_repeated_points=True)
    return points

def not_lineal_transformation(v1 , v2 , v3 , x_k):
    '''
    Moves the point x_k from a triangle with vertices [[0,0],[1,0],[0,1]]
    to a point in the triangle built by v1,v2,v3
    '''
    X_k = v1 + (v2-v1)*x_k[0] + (v3-v1)*x_k[1]
    return X_k
        
def check_contained_triangles_alternative(mesh_1 , mesh_2 , N_ref , assume_uniform = True):
    
    face_1 , vert_1 = np.transpose(mesh_1.leaf_view.elements) , np.transpose(mesh_1.leaf_view.vertices)
    face_2 , vert_2 = np.transpose(mesh_2.leaf_view.elements) , np.transpose(mesh_2.leaf_view.vertices)
    
    normal_1 , normal_2 = normals(face_1 , vert_1 ) , normals(face_2 , vert_2)
        
    relationship = np.zeros((len(face_2),)) - 1
    
    unitary_ref_points  = unitary_refinement_generator(N_ref)
    
    c1=0
    for f1 in face_1:
        
        v1_1 , v1_2 , v1_3 = vert_1[f1[0]] , vert_1[f1[1]] , vert_1[f1[2]]
        
        n_1 = normal_1[c1]
        
        aux_points = np.empty([0,3])
        for x_k in unitary_ref_points:
            aux_points = np.vstack([aux_points ,
                                    not_lineal_transformation(v1_1 , v1_2 , v1_3 , x_k) ])      
        
        
        c2 = 0    
        for f2 in face_2:
            
            v2_1 , v2_2 , v2_3 = vert_2[f2[0]] , vert_2[f2[1]] , vert_2[f2[2]]
            
            n_2 = normal_2[c2]
            
            if relationship[c2]<0:
                
                if np.linalg.norm(n_1-n_2)<10**-8:
                    
                    if similarities_counter(np.array([v2_1 , v2_2 , v2_3]) , aux_points ):                  
                        relationship[c2] = c1
                            
            c2+=1
        c1+=1
        
    relationship = relationship.astype(int)
    return relationship
