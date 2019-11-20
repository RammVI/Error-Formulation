
import meshio
import os

path = os.path.join('Molecule' , 'sphere_cent' , '30' , 'True')

file_name = 'sphere_cent_0-s5'

mesh = meshio.read( os.path.join(path , file_name + '_aux_improved_0.off') )

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