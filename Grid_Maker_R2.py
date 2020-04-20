
# Last revision (2) 26-May-2019

import bempp.api, numpy as np, time, os, matplotlib.pyplot as plt
from math import pi
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz
from matplotlib import pylab as plt
from File_converter_python2 import *

# This python must be saved in a directory where you have a folder named
# /Molecule/Molecule_Name, as obviusly Molecule_Name holds a .pdb or .pqr file
# ; otherwise, this function won't do anything.

# Data is saved in format {mol_name}_{mesh_density}-{it_count}
# Where mol_name is the abreviated name of the molecule
# mesh_density is the density of the mesh in elements per square amstrong
# it_count is for the mesh ref pluggin and will be treaten as -0 when is the first grid made

# IMPORTANT BUGS - 1. NANOSHAPER MUST BE REPAIRED - ONLY MSMS ALLOWED
# 2. print('juan') not printing in xyzr_to_msh !!!!!!!!!!!! This may be a reason
#    why .msh file is not being created.

#With the .pdb file we can build .pqr & .xyzr files, and they don't change when the mesh density is changed.
def pdb_to_pqr(mol_name , stern_thickness , method = 'amber' ):
    '''
    Function that makes .pqr file from .pdb using Software/apbs-pdb2pqr-master/pdb2pqr/main.py
    Be careful of the version and the save directory of the pdb2pqr python shell.
    mol_name : Abreviated name of the molecule
    stern_thicness : Length of the stern layer
    method         : This parameter is an 
    '''
    path = os.getcwd()
        
    pdb_file , pdb_directory = mol_name+'.pdb' , os.path.join('Molecule',mol_name)
    pqr_file , xyzr_file     = mol_name+'.pqr' , mol_name+'.xyzr'
    
    if os.path.isfile(os.path.join('Molecule',mol_name,pqr_file)):
        print('File already exists in directory.')
        return None
    
    # The apbs-pdb2pqr rutine, allow us to generate a .pqr file
    pdb2pqr_dir = os.path.join('Software','apbs-pdb2pqr-master','pdb2pqr','main.py')
    exe=('python2.7  ' + pdb2pqr_dir + ' '+ os.path.join(pdb_directory,pdb_file) +
         ' --ff='+method+' ' + os.path.join(pdb_directory,pqr_file)   )
    
    os.system(exe)
    
    # Now, .pqr file contains unneeded text inputs, we will save the rows starting with 'ATOM'.
    
    pqr_Text = open( os.path.join(pdb_directory , pqr_file) ).read()
    pqr_Text_xyzr = open(os.path.join(pdb_directory , xyzr_file )  ,'w+')

    
    for i in pqr_Text.split('\n'):
        row=i.split()
        if row[0]=='ATOM':
            aux=row[5]+' '+row[6]+' '+row[7]+' '+row[-1]
            pqr_Text_xyzr.write(aux + '\n')   
    pqr_Text_xyzr.close()
    
    print('Global .pqr & .xyzr ready.')
    
    # The exterior interface is easy to add, by increasing each atom radii
    if stern_thickness>0: 
        xyzr_file_stern = os.path.join(pdb_directory , mol_name +'_stern.xyzr')
        pqr_Text_xyzr_s = open(xyzr_file_stern ,'w')
        
        for i in pqr_Text.split('\n'):
            row=i.split()
            if row[0]=='ATOM':
                R_vv=float(row[-1])+stern_thickness
                pqr_Text_xyzr_s.write(row[5]+' '+row[6]+' '+row[7]+' '+str(R_vv)+'\n' )      
        pqr_Text_xyzr_s.close()
        print('Global _stern.pqr & _stern.xyzr ready.')
    
    return 

def pqr_to_xyzr(mol_name , stern_thickness , method = 'amber' ):
    '''
    Extracts .xyzr information from .pqr
    mol_name : Abreviated name of the molecule
    stern_thickness : Length of the stern layer
    method          : amber by default , a pdb2pqr parameter to build the mesh.
    '''
    path = os.getcwd()
        
    pqr_directory = os.path.join('Molecule',mol_name)
    pqr_file , xyzr_file     = mol_name+'.pqr' , mol_name+'.xyzr'
     
    # Now, .pqr file contains unneeded text inputs, we will save the rows starting with 'ATOM'.
    
    pqr_Text = open( os.path.join(pqr_directory , pqr_file) ).read()
    pqr_Text_xyzr = open(os.path.join(pqr_directory , xyzr_file )  ,'w+')

    
    for i in pqr_Text.split('\n'):
        row=i.split()
        if len(row)==0: continue
            
        if row[0]=='ATOM':
            aux=' '.join( [row[5],row[6],row[7],row[-1]] )
            pqr_Text_xyzr.write(aux + '\n')   
    pqr_Text_xyzr.close()
    
    print('.xyzr File from .pqr ready.')
    
    # The exterior interface is easy to add, by increasing each atom radii
    if stern_thickness>0: 
        xyzr_file_stern = os.path.join(pqr_directory , mol_name +'_stern.xyzr')
        pqr_Text_xyzr_s = open(xyzr_file_stern ,'w')
        
        for i in pqr_Text.split('\n'):
            row=i.split()
            if row[0]=='ATOM':
                R_vv=float(row[-1])+stern_thickness
                pqr_Text_xyzr_s.write(row[5]+' '+row[6]+' '+row[7]+' '+str(R_vv)+'\n' )      
        pqr_Text_xyzr_s.close()
        print('Global _stern.pqr & _stern.xyzr ready.')
        
    return None

def NanoShaper_config(xyzr_file , dens , probe_radius):
    '''
    Yet in beta version. Changes some data to build the mesh with NanoShaper
    xyzr_file : Directory of the xyzr_file
    dens      : mesh density
    probe_radius : might be set to 1.4
    '''
    t1 = (  'Grid_scale = {:s}'.format(str(dens)) 
                #Specify in Angstrom the inverse of the side of the grid cubes  
              , 'Grid_perfil = 80.0 '                     
                #Percentage that the surface maximum dimension occupies with
                # respect to the total grid size,
              , 'XYZR_FileName = {:s}'.format(xyzr_file)  
              ,  'Build_epsilon_maps = false'              
              , 'Build_status_map = false'                
              ,  'Save_Mesh_MSMS_Format = true'            
              ,  'Compute_Vertex_Normals = true'           
              ,  'Surface = ses  '                         
              ,  'Smooth_Mesh = true'                      
              ,  'Skin_Surface_Parameter = 0.45'           
              ,  'Cavity_Detection_Filling = false'        
              ,  'Conditional_Volume_Filling_Value = 11.4' 
              ,  'Keep_Water_Shaped_Cavities = false'      
              ,  'Probe_Radius = {:s}'.format( str(probe_radius) )                
              ,  'Accurate_Triangulation = true'           
              ,  'Triangulation = true'                    
              ,  'Check_duplicated_vertices = true'        
              ,  'Save_Status_map = false'                 
              ,  'Save_PovRay = false'                     )
    return t1

def xyzr_to_msh(mol_name , dens , probe_radius , stern_thickness , min_area , Mallador = 'MSMS',
               suffix = ''):
    '''
    Makes msh (mesh format for BEMPP) from xyzr file
    mol_name : Abreviated name of the molecule
    dens     : Mesh density
    probe_radius : might be set to 1.4[A]
    stern_thickness : Length of the stern layer
    min_area        : Discards elements with less area than this value
    Mallador        : MSMS or NanoShaper
    
    outputs : Molecule/{mol_name}/{mol_name}_{dens}-0.msh
    Where -0 was added because of probable future mesh refinement and easier handling of meshes.
    '''

    path = os.getcwd()
    mol_directory = os.path.join('Molecule',mol_name)
    path = os.path.join(path , mol_directory)
    xyzr_file     = os.path.join(mol_directory, mol_name + '.xyzr') 
    
    if stern_thickness > 0:  xyzr_s_file = os.path.join(mol_directory , mol_name + '_stern.xyzr'  )
    
    # The executable line must be:
    #  path/Software/msms/msms.x86_64Linux2.2.6.1 
    # -if path/mol_name.xyzr       (Input File)
    # -of path/mol_name -prob 1.4 -d 3.    (Output File)
   
    # The directory of msms/NS needs to be checked, it must be saved in the same folder that is this file
    if Mallador == 'MSMS':    
        M_path = os.path.join('Software','msms','msms.x86_64Linux2.2.6.1')
        mode= ' -no_header'
        prob_rad, dens_msh = ' -prob ' + str(probe_radius), ' -d ' + str(dens)
        exe= (M_path
              +' -if ' + xyzr_file
              +' -of ' + os.path.join(mol_directory , mol_name )+'_{0:s}-0'.format( str(dens) )
              + prob_rad  + dens_msh + mode )
        os.system(exe)
        print('Normal .vert & .face Done')

        grid = factory_fun_msh( mol_directory , mol_name , min_area , dens , Mallador , suffix = '-0')
        print('Normal .msh Done')
        
        # As the possibility of using a stern layer is available:
        if stern_thickness > 0:
            prob_rad, dens_msh = ' -prob ' + str(probe_radius), ' -d ' + str(dens)
            exe= (M_path+' -if '  + xyzr_s_file + 
              ' -of ' + mol_directory + mol_name +'_stern' + prob_rad  + dens_msh  + mode )
            os.system(exe)
            print('Stern .vert & .face Done')
            stern_grid= factory_fun_msh( mol_directory , mol_name+'_stern', min_area )
            print('Stern .msh Done')
        
    elif Mallador == 'NanoShaper': 
        Mallador= os.path.join('Software','nanoshaper','NanoShaper')
        config  = os.path.join('Software','nanoshaper','config')
        
        # NanoShaper options can be changed from the config file
        Config_text = open(config,'w')
        
        Conf_text = NanoShaper_config(xyzr_file , dens , probe_radius)
        
        Config_text.write('\n'.join(Conf_text))
        Config_text.close()
        
        # Files are moved to the same directory
        os.system(' '.join( ('./'+Mallador,config)  ))
        
        Files = ('triangleAreas{0:s}.txt','triangulatedSurf{0:s}.face' ,'triangulatedSurf{0:s}.vert',
         'exposedIndices{0:s}.txt','exposed{0:s}.xyz  ' , 'stderror{0:s}.txt' )
        for f in Files:
            os.system(' '.join( ('mv ', f.format('') 
                                 , os.path.join(mol_directory,f.format('_'+str(dens))))))
            
        if not os.path.isfile(os.path.join(path , 'triangulatedSurf{0:s}.vert'.format('_'+str(dens)))):
            print('Algo salio mal! .vert no se ha creado')
            return 'Error'
        
        print('Normal .vert & .face Done')
        
        grid = factory_fun_msh( mol_directory , mol_name , min_area , dens , Mallador = 'NanoShaper' )
        print('msh File Done')
        
        if stern_thickness>0:
            # NanoShaper options can be changed only in the config file
            Config_text = open(config,'w')
            Conf_text = NanoShaper_config(xyzr_file_stern , dens , probe_radius)
            
            Config_text.write('\n'.join(Conf_text))
            Config_text.close()

            # Files are moved to the same directory
            os.system(' '.join( ('./'+Mallador,config)  ))

            Files = ('triangleAreas{0:s}.txt','triangulatedSurf{0:s}.face' ,'triangulatedSurf{0:s}.vert',
             'exposedIndices{0:s}.txt','exposed{0:s}.xyz  ' , 'stderror{0:s}.txt' )
            for f in Files:
                os.system(' '.join( ('mv ', f.format('') 
                                     , os.path.join(mol_directory,f.format('_stern_'+str(dens)))))) 
            print('Stern .vert & .face Done')
            
            stern_grid= factory_fun_msh( mol_directory , mol_name+'_stern', min_area , dens , Mallador = 'NanoShaper')
            print('stern_.msh File Done')
    print('Mesh Ready')
    return

def factory_fun_msh( mol_directory , mol_name , min_area , dens , Mallador , suffix = ''):
    '''
    This functions builds msh file adding faces and respective vertices.
    mol_directory : Directory of the molecule
    mol_name      : Abreviated name of the molecule
    min_area      : Min. area set to exclude small elements
    dens          : mesh density
    Mallador      : MSMS - NanoShaper or Self (if doing the GOMR)
    suffix        : Suffix of the .vert and .face file after the mesh density ({mol_name}_{d}{suffix})
                    might be used as -{it_count}
    '''
    # Factory function for creating a .msh file from .vert & .face files
    factory = bempp.api.grid.GridFactory()
    
    # .vert and .face files are readed    
    if Mallador == 'MSMS':
        vert_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}{1}.vert'.format(str(dens),suffix) ),'r' ).read().split('\n')
        face_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}{1}.face'.format(str(dens),suffix) ),'r' ).read().split('\n')
    elif Mallador == 'NanoShaper':
        vert_Text = open( os.path.join(mol_directory , 'triangulatedSurf_{0:s}.vert'.format(str(dens)) ),'r' ).read().split('\n')
        face_Text = open( os.path.join(mol_directory , 'triangulatedSurf_{0:s}.face'.format(str(dens)) ),'r' ).read().split('\n')
    elif Mallador == 'Self':
        vert_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}{1}.vert'.format(str(dens),suffix) ),'r' ).read().split('\n')
        face_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}{1}.face'.format(str(dens),suffix) ),'r' ).read().split('\n')
    
    xcount, atotal, a_excl = 0, 0., 0.
    vertex = np.empty((0,3))
    

    # Lets load the bempp plugging to add elements
    factory = bempp.api.grid.GridFactory()
    
    # Grid assamble
    
    # Let's separe by mallator criteria because of diferent handling of info
    if Mallador !='Self': 
        
        for line in vert_Text:
            line = line.split()
            if len(line) != 9: continue
            vertex = np.vstack(( vertex, np.array([line[0:3]]).astype(float) ))
            factory.insert_vertex(vertex[-1])

        
        for line in face_Text:
            line = line.split()
            if len(line)!=5 : continue
            A, B, C, _, _ = np.array(line).astype(int)
            side1, side2  = vertex[B-1]-vertex[A-1], vertex[C-1]-vertex[A-1]
            face_area = 0.5*np.linalg.norm(np.cross(side1, side2))
            atotal += face_area
            if face_area > min_area:
                factory.insert_element([A-1, B-1, C-1])
            else:
                xcount += 1.4        
                a_excl += face_area 
    
    elif Mallador == 'Self':
        for line in vert_Text[:-1]:
            line = line.split()
            vertex = np.vstack( ( vertex, np.array(line).astype(float) )  )
            factory.insert_vertex(vertex[-1])
        
        for line in face_Text[:-1]:
            line = line.split()
            A, B, C = np.array(line).astype(int)
            side1, side2  = vertex[B-1]-vertex[A-1], vertex[C-1]-vertex[A-1]
            face_area = 0.5*np.linalg.norm(np.cross(side1, side2))
            atotal += face_area
            if face_area > min_area:
                factory.insert_element([A-1, B-1, C-1])
            else:
                xcount += 1.4        
                a_excl += face_area 
                
    
    grid = factory.finalize()
    
    export_file = os.path.join(mol_directory , mol_name +'_'+str(dens)+ suffix +'.msh' )
    print(export_file)
    bempp.api.export(grid=grid, file_name=export_file) 
    
    return grid

def triangle_areas(mol_directory , mol_name , dens , return_data = False , suffix = '', Self_build = False):
    """
    This function calculates the area of each element.
    Avoid using this with NanoShaper, only MSMS recomended
    Self_build : False if using MSMS or NanoShaper - True if building with new methods
    Has a BUG! probably not opening .vert or .face or not creating .txt or both :P .
    """
    
    vert_Text = open( os.path.join(mol_directory , mol_name +'_'+str(dens)+suffix+'.vert' ) ).read().split('\n')
    face_Text = open( os.path.join(mol_directory , mol_name +'_'+str(dens)+suffix+'.face' ) ).read().split('\n')
    area_list = np.empty((0,1))
    area_Text = open( os.path.join(mol_directory , 'triangleAreas_'+str(dens)+suffix+'.txt' ) , 'w+')
    
    vertex = np.empty((0,3))
    
    if not Self_build:
        for line in vert_Text:
            line = line.split()
            if len(line) !=9: continue
            vertex = np.vstack(( vertex, np.array(line[0:3]).astype(float) ))

        atotal=0.0
        # Grid assamble
        for line in face_Text:
            line = line.split()
            if len(line)!=5: continue
            A, B, C, _, _ = np.array(line).astype(int)
            side1, side2  = vertex[B-1]-vertex[A-1], vertex[C-1]-vertex[A-1]
            face_area = 0.5*np.linalg.norm(np.cross(side1, side2))

            area_Text.write( str(face_area)+'\n' )

            area_list = np.vstack( (area_list , face_area ) )
            atotal += face_area

        area_Text.close()

        if return_data:
            return area_list
        
    elif Self_build:
        
        for line in vert_Text[:-1]:
            line = line.split()
            
            vertex = np.vstack(( vertex, np.array(line[0:3]).astype(float) ))

        atotal=0.0
        # Grid assamble
        for line in face_Text[:-1]:
            line = line.split()
            A, B, C = np.array(line[0:3]).astype(int)
            side1, side2  = vertex[B-1]-vertex[A-1], vertex[C-1]-vertex[A-1]
            face_area = 0.5*np.linalg.norm(np.cross(side1, side2))
            area_Text.write( str(face_area)+'\n' )

            area_list = np.vstack( (area_list , face_area ) )
            atotal += face_area

        area_Text.close()

        if return_data:
            return area_list
    
    return None

def normals_to_element( face_array , vert_array , check_dir = False ):
    '''
    Calculates normals to a given element, pointint outwards.
    face_array : Array of vertex position for each triangle
    vert_array : Array of vertices
    check_dir  : checks direction of normals. WORKS ONLY FOR A SPHERE WITH RADII 1!!!!!!!!!!!!
    '''

    normals = np.empty((0,3))
    element_cent = np.empty((0,3))
    
    check_list = np.empty((0,1))
    
    for face in face_array:
        
        f1,f2,f3 = face-1
        v1 , v2 , v3 = vert_array[f1] , vert_array[f2] , vert_array[f3]
        n = np.cross( v2-v1 , v3-v1 ) 
        normals = np.vstack((normals , n/np.linalg.norm(n) )) 
        element_cent = np.vstack((element_cent, (v1+v2+v3)/3. ))
        
        if check_dir:
            v_c = v1 + v2 + v3
            pdot= np.dot( v_c , n )
            if pdot>0:
                check = True
            else:
                check = False
            check_list = np.vstack( (check_list , check ) )
            

    return normals , check_list[:,0]


def vert_and_face_arrays_to_text_and_mesh(mol_name , vert_array , face_array , suffix 
                                          , dens=2.0 , Self_build=True):
    '''
    This rutine saves the info from vert_array and face_array and creates .msh and areas.txt files
    mol_name : Abreviated name for the molecule
    dens     : Mesh density, anyway is not a parameter, just a name for the file
    vert_array: array containing verts
    face_array: array containing verts positions for each face
    suffix    : text added to diference the meshes.
    
    Returns None but creates Molecule/{mol_name}/{mol_name}_{mesh_density}{suffix}.msh file.
    '''
    normalized_path = os.path.join('Molecule',mol_name,mol_name+'_'+str(dens)+suffix)
    
    vert_txt = open( normalized_path+'.vert' , 'w+' )
    for vert in vert_array:
        txt = ' '.join( vert.astype(str) )
        vert_txt.write( txt + '\n')
    vert_txt.close()
    
    face_txt = open( normalized_path+'.face' , 'w+' )
    for face in face_array:
        txt = ' '.join( face.astype(int).astype(str) )
        face_txt.write( txt + '\n')
    face_txt.close()
    
    mol_directory = os.path.join('Molecule',mol_name)
    min_area = 0

    factory_fun_msh( mol_directory , mol_name , min_area , dens , Mallador='Self', suffix=suffix)
    #triangle_areas(mol_directory , mol_name , str(dens) , suffix = suffix , Self_build = Self_build)
    
    return None

def Grid_loader(mol_name , mesh_density , suffix , GAMer=False):
    
    path = os.path.join('Molecule',mol_name)
    grid_name_File =  os.path.join(path,mol_name + '_'+str(mesh_density)+suffix+'.msh')
    
    print(grid_name_File)
    
    if os.path.isfile(grid_name_File) and suffix == '-0':
        
        pqr_directory = os.path.join('Molecule',mol_name, mol_name+'.pqr' )
        
        if not os.path.isfile(pqr_directory):
            pdb_to_pqr(mol_name , stern_thickness , method = 'amber' )
       
    if suffix == '-0':
        pqr_to_xyzr(mol_name , stern_thickness=0 , method = 'amber' )
        xyzr_to_msh(mol_name , mesh_density , 1.4 , 0 , 0
                    , Mallador = 'MSMS', suffix = suffix )
        
    print('Working on '+grid_name_File )
    grid = bempp.api.import_grid(grid_name_File)
    if GAMer:
        face_array = np.transpose(grid.leaf_view.elements)+1
        vert_array = np.transpose(grid.leaf_view.vertices)
        
        new_face_array , new_vert_array = Improve_Mesh(face_array , vert_array , path , mol_name + '_'+str(mesh_density)+suffix )
        
        vert_and_face_arrays_to_text_and_mesh(mol_name , new_vert_array , new_face_array , suffix 
                                          , dens=mesh_density , Self_build=True)
        
        grid = bempp.api.import_grid(grid_name_File)
    
    return grid





    
