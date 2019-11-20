
# Last revision May 26 15:12

import bempp.api, numpy as np, time, os, matplotlib.pyplot as plt
from math import pi
from bempp.api.operators.boundary import sparse, laplace, modified_helmholtz
from matplotlib import pylab as plt

# This python must be saved in a directory where you have a folder named
# /Molecule/Molecule_Name, as obviusly Molecule_Name holds a .pdb or .pqr file
# ; otherwise, this function won't do anything.

# With the .pdb file we can build .pqr & .xyzr files, and they don't change when the mesh density is changed.
def pdb_to_pqr(mol_name , stern_thickness , method = 'amber' ):
    
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
        
    return

def NanoShaper_config(xyzr_file , dens , probe_radius):
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
    

def xyzr_to_msh(mol_name , dens , probe_radius , stern_thickness , min_area , Mallador = 'NanoShaper'):

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

        grid = factory_fun_msh( mol_directory , mol_name , min_area , dens , Mallador)
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
    # Factory function for creating a .msh file from .vert & .face files
    
    factory = bempp.api.grid.GridFactory()
    
    if Mallador == 'MSMS':
        vert_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}.vert'.format(str(dens)) ) ).read().split('\n')
        face_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}.face'.format(str(dens)) ) ).read().split('\n')
    elif Mallador == 'NanoShaper':
        vert_Text = open( os.path.join(mol_directory , 'triangulatedSurf_{0:s}.vert'.format(str(dens)) ) ).read().split('\n')
        face_Text = open( os.path.join(mol_directory , 'triangulatedSurf_{0:s}.face'.format(str(dens)) ) ).read().split('\n')
    elif Mallador == 'Self':
        vert_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}{1}.vert'.format(str(dens),suffix) ) ).read().split('\n')
        face_Text = open( os.path.join(mol_directory , mol_name +'_{0:s}{1}.face'.format(str(dens),suffix) ) ).read().split('\n')
    
    xcount, atotal, a_excl = 0, 0., 0.
    vertex = np.empty((0,3))

    factory = bempp.api.grid.GridFactory()
    # Create the grid with the factory method
    
    # Grid assamble
    
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
    bempp.api.export(grid=grid, file_name=export_file) 
    
    return grid

def triangle_areas(mol_directory , mol_name , dens , save = False,suffix = '', Self_build = False):
    """
    This function calculates the area of each element.
    Avoid using this with NanoShaper, only MSMS recomended
    Self_build : False if using MSMS or NanoShaper - True if building with new methods
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

        if save:
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
            print(A,B,C)
            area_Text.write( str(face_area)+'\n' )

            area_list = np.vstack( (area_list , face_area ) )
            atotal += face_area

        area_Text.close()

        if save:
            return area_list
    
    return None
    
def vert_and_face_arrays_to_text_and_mesh(mol_name , vert_array , face_array , suffix , dens=2.0 , Self_build=True):
    '''
    This rutine saves the info from vert_array and face_array and creates .msh and areas.txt files
    mol_name : Abreviated name for the molecule
    dens     : Mesh density, anyway is not a parameter, just a name for the file
    vert_array: array containing verts
    face_array: array containing verts positions for each face
    suffix    : text added to diference the meshes
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
    triangle_areas(mol_directory , mol_name , str(dens) , suffix = suffix , save = True , Self_build = Self_build)
    
    return None
    
def centroide( vert_list ):
    #  Function that returns the mid point 
    #   and the volume of each element from
    #   a list of vertex
    
    A , B , C , D = vert_list
    
    centro = (A + B + C + D)*0.25
    u , v , w = (A - B) , (B - C) , (D - C)
    Vol =np.absolute( np.dot( np.cross(u,v) , w ) ) / 6.
    return centro , Vol
    
    
def text_XML_format_To_centroid_and_Vol(mol_path):
    
    path = os.getcwd()
    
    text_XML = open( os.path.join(path , mol_path) , 'r' ).read().split('\n')
    
    vert_l = np.empty((0,3))
    vol_l  = np.empty((0,4))
    
    for line in text_XML:
        line_aux = line.split()
        if len(line_aux)>5:
            
            if line_aux[2][0] == 'x':
                x = float( line_aux[2][4:-1] )
                y = float( line_aux[3][4:-1] )
                z = float( line_aux[4][4:-1] )
                
                vert = np.array( (x,y,z) )
                
                vert_l = np.vstack( (vert_l , vert) )
        
            elif line_aux[2][0] == 'v':
            
                v0 = int(line_aux[2][4:-1])
                v1 = int(line_aux[3][4:-1])
                v2 = int(line_aux[4][4:-1])
                v3 = int(line_aux[5][4:-1])
                
                vol   = np.array( (v0,v1,v2,v3) )
                vol_l = np.vstack( (vol_l , vol) )

    centr  = np.empty((0,3))
    vol    = 0.
    for el in vol_l:
        verts = [ vert_l[index] for index in el]
        cent_i , vol_i = centroide( verts )        
        
        centr = np.vstack( (centr , cent_i) )
        vol = vol + vol_i
    
    print(vol)
        
    return centr,vol

def centr_list_to_txt(mol_name , centr , Vol=0 , Centroide = True ):
    """
    This function transform the centroid list to a txt file
    be careful, because the FIRST line is the Volume of the elements
    """
    
    if Centroide:
        
        path = os.getcwd()

        centr_text = open(os.path.join(path,'Molecule' , mol_name , mol_name)+'_centr.txt'  ,'w+' )

        centr_text.write(str(Vol)+'\n')

        for c in centr:
            centr_text.write( ' '.join(map(str , c ) )+'\n' )
        centr_text.close()
        
        return None
    
    else:
        print('Error')
        return None    