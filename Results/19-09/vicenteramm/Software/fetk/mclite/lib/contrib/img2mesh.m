function img2mesh(name,key1,key2);
%
% Usage ==> img2mesh('image_file')
% Input ==> name of 256 color gray scale file (with path and filename as string)
%           examples: img2mcsf('bridge.jpg') or img2mcsf('.\images\bowl.gif')
% Output ==>file named mcin.m containing variables vert[?,5] and simp[?,9])
%           placed into current matlab directory.
%
%This function turns an 8-bit or 24-bit, 256 color grayscale image file into an 
%mcsf format file.  It will accept any image file that Matlab recognizes 
%(bmp,jpg,gif,tiff,etc). Throughout this description, color=0 is black, RGB=(0,0,0) and
%color=255 is white, RGB=(255,255,255), while color=x is off white or black, RGB=(x,x,x).
%
%Each pixel of color=0 is turned into a square with four vertices and 2 simplices
%Each pixel of color=1 is turned into a triangle with 3 vertices and 1 simplex
%  (unless that pixel has only one neighbor then it will be left as a square.)
%*Note: keep the number of pixels small. Let MClite improve your grid's resolution.
%
%Boundary types are identified by the contrasting color of neighboring pixels.
%A neighbor of color=255 yields a Neuman (type=2) boundary on your black / off black pixel.
%Neighbor colors of 225<=c<255 give you a Dirichlet (type=1) boundary.
%Neighbor colors of 200<=c<225 give you a Neuman (type=4) boundary.
%Neighbor colors of 175<=c<200 give you a Neuman (type=6) boundary.
%150<=c<175:(type=8), 125<=c<150:(type=10), 100<=c<125:(type=12).
%Use these special Neuman types for applying mixed boundary conditions.
%
%If a pixel of color 100<=c<255 touches multiple interior pixels and you only want the 
%special boundary to be present on some of the interior pixels, make those interior pixels
%off-black, 10<c<100. (or if it was a triangle=1 pixel use offblack=6)
%*Note: The neighbor to a triangle=1 pixel is always the pixel to its left or right.
%
%To view the mcin.m file in the current directory, type img2mesh('view').  
%Optional features are activated by keys: img2mesh(name,key1,key2).  Set key1==1, to display 
%mesh right away and key2==1, to supress opening dialog and mcin.m existence warning.
%
%For additional help, see the associated pdf help file.  
%Enjoy, and email any questions and/or comments to:  CDeotte@math.ucsd.edu


%BEGIN: NAME PARAMETER MANAGEMENT.
  try %Set default filename if none provided.
      name;
  catch
      name='.\images\default.bmp';
  end;

  if strcmpi(name,'view')
      %mcin;
      try
      [vert,simp]=read(0);
      draw(vert,simp);
      catch
          fprintf('\nYou need to install mclite and run setpath first to use this function.\n');
          fprintf('If you don''t have mclite, then use img2mesh(image_file,1) That will create\n');
          fprintf('an mcin.m file and display the mesh and the same time.\n\n');
      end;
      return;
  end;

  if ~exist(name)
      fprintf('\nMatlab cannot locate ''%s'' .  Try including the path and/or extension.\n\n',name);
      return;
  end;
%END: NAME PARAMETER MANAGEMENT.

%BEGIN: KEY2 PARAMETER MANGEMENT.
  msg=1;
  try
  if key2==1
      msg=0;
  end;
  end;

  if msg %display info message if key2==1 is not present

  if exist('mcin.m')
  fprintf('\nWARNING: mcin.m already exists in the current directory.\n');  
  user=input('Do you wish to overwrite it?(Y or N)(Just hitting enter is a yes)','s');
  if isempty(user) user='y'; end;
  if ~strcmpi(user,'y')
      fprintf('\nProgram aborted.\n\n');
      return;
  end;
  end;

  fprintf('\nThe mcsf file ''mcin.m'' was created from ''%s''\n',name)
  fprintf('Type ''img2mesh(''view'')'' to view your mesh file.\n\n');
  end;
%END: KEY2 PARAMETER MANAGEMENT.

%Read in image and initialize variables.
pixels=imread(name);
[m,n,k]=size(pixels);
if k>1 pixels=pixels(:,:,1);end;
vertices=zeros(m,n,4);
vert=zeros(1,5);
simp=zeros(1,9);
vid=0;
sid=0;
usesoffblack=0;

%Image file error checking and correction.
for i=1:m
    for j=1:n
        if pixels(i,j)<100 && pixels(i,j)>5 usesoffblack=1;
        elseif pixels(i,j)>255 
            fprintf('Your image file has color values greater than 255.\n');
            fprintf('Color values should be between 0 and 255 inclusive.\n');
            return;
        elseif pixels(i,j)<0
            fprintf('Your image file has color values less than 0.\n');
            fprintf('Color values should be between 0 and 255 inclusive.\n');
            return;  
        end;
    end;
end;

%Standardize offblack & non-offback usage.
%color=0 black pixel
%1<=color<=5 black special pixel
%6<=color<=10 off black special pixel
%11<=color<=99 off black pixel
%100<=color<=175 available (but considered border pixel)
%176<=color<=254 off white pixel
%color=255 white pixel
if usesoffblack
for i=1:m
    for j=1:n
        if pixels(i,j)>=100 && pixels(i,j)~=255 %then search for an off black neighbor
            yes=0;
            if gpixels(i-1,j)>5 && gpixels(i-1,j)<100 yes=1;end;
            if gpixels(i,j+1)>5 && gpixels(i,j+1)<100 yes=1;end;
            if gpixels(i,j-1)>5 && gpixels(i,j-1)<100 yes=1;end;
            if gpixels(i+1,j)>5 && gpixels(i+1,j)<100 yes=1;end;
            if ~yes %make offblack neighbors
                if gpixels(i-1,j)==0 pixels(i-1,j)=11;end;
                if gpixels(i,j+1)==0 pixels(i,j+1)=11;end;
                if gpixels(i,j-1)==0 pixels(i,j-1)=11;end;
                if gpixels(i+1,j)==0 pixels(i+1,j)=11;end;
                if gpixels(i-1,j)<=5 pixels(i-1,j)=pixels(i-1,j)+5;end;
                if gpixels(i,j+1)<=5 pixels(i,j+1)=pixels(i,j+1)+5;end;
                if gpixels(i,j-1)<=5 pixels(i,j-1)=pixels(i,j-1)+5;end;
                if gpixels(i+1,j)<=5 pixels(i+1,j)=pixels(i+1,j)+5;end;
            end;
        end;
    end;
end;
end;
            
% neighbor identification throughout file.
% cell6 cell3 cell7
% cell2 (r,c) cell4
% cell9 cell5 cell8

%Record vertices.  Interpret each pixel's corner as a vertex.
%(Unless pixel color is 1, whereas interpret as possible triangle only)
for row=1:m
    for col=1:n
        cell6=1;cell3=1;cell7=1;cell4=1;  
        bl=1;br=1;tl=1;tr=1;cell2=1;cell5=1;
        %determine which neighbors are present.
        %Everyone writes in their bottom left vertex.  
        %Then if unwritten vertices exist, a pixel writes in its
        %bottom right vertex, next top left, next top right.
        if gpixels(row-1,col-1)>100 cell6=0; end;
        if gpixels(row-1,col)>100 cell3=0; end;
        if gpixels(row-1,col+1)>100 cell7=0; end;
        if gpixels(row,col+1)>100 cell4=0; end;
        
        if pixels(row,col)==1 || pixels(row,col)==6 %triangle pixel
            if gpixels(row,col-1)>100 cell2=0; end;
            if gpixels(row+1,col)>100 cell5=0; end;
            %Don't turn pixel into triangle unless >=2 neighbors
            if cell2+cell3+cell4+cell5>1 
                 if cell4==0 && cell5==0 br=0; 
                 elseif cell2==0 && cell3==0 tl=0; 
                 elseif cell5==0 && cell2==0 bl=0; 
                 elseif cell3==0 && cell4==0 tr=0; end;
            end;
        end;
        
        if pixels(row,col)<100 %pixel is present
          if bl wvert(vid,0,col-1,m-row,0);end; %write bottom left vertex
          if (br && ~cell4) wvert(vid,0,col,m-row,0); end; %bottom right (if color=1 no)
          if (tl && ~cell3 && ~cell6) wvert(vid,0,col-1,m-row+1,0); end; %top left (if color=2 no)
          if (tr && ~cell7 && ~cell3 && ~cell4) wvert(vid,0,col,m-row+1,0); end; %top right
        end;
    end;
end;

%Record each simplex. Making note of boundary type.
for row=1:m
    for col=1:n
        if pixels(row,col)<100
          %determine neighbors to cell  
          tri1=1;tri2=1;tri3=0;tri4=0;
          cell2=0;cell3=0;cell5=0;cell4=0;mid=0;mat=0;
          if gpixels(row+1,col)==255 cell5=2; end;
          if gpixels(row-1,col)==255 cell3=2; end;
          if gpixels(row,col-1)==255 cell2=2; end;
          if gpixels(row,col+1)==255 cell4=2; end;
          %decide if a special border is present
          if usesoffblack && pixels(row,col)<6
              if gpixels(row+1,col)>=100 && cell5==0 cell5=2; end;
              if gpixels(row-1,col)>=100 && cell3==0 cell3=2; end;
              if gpixels(row,col-1)>=100 && cell2==0 cell2=2; end;
              if gpixels(row,col+1)>=100 && cell4==0 cell4=2; end;
          end;
              if gpixels(row+1,col)>=225 && cell5==0 cell5=1; end;
              if gpixels(row-1,col)>=225 && cell3==0 cell3=1; end;
              if gpixels(row,col-1)>=225 && cell2==0 cell2=1; end;
              if gpixels(row,col+1)>=225 && cell4==0 cell4=1; end;
              
              if gpixels(row+1,col)>=200 && cell5==0 cell5=4; end;
              if gpixels(row-1,col)>=200 && cell3==0 cell3=4; end;
              if gpixels(row,col-1)>=200 && cell2==0 cell2=4; end;
              if gpixels(row,col+1)>=200 && cell4==0 cell4=4; end;
              
              if gpixels(row+1,col)>=175 && cell5==0 cell5=6; end;
              if gpixels(row-1,col)>=175 && cell3==0 cell3=6; end;
              if gpixels(row,col-1)>=175 && cell2==0 cell2=6; end;
              if gpixels(row,col+1)>=175 && cell4==0 cell4=6; end;
              
              if gpixels(row+1,col)>=150 && cell5==0 cell5=8; end;
              if gpixels(row-1,col)>=150 && cell3==0 cell3=8; end;
              if gpixels(row,col-1)>=150 && cell2==0 cell2=8; end;
              if gpixels(row,col+1)>=150 && cell4==0 cell4=8; end;
              
              if gpixels(row+1,col)>=125 && cell5==0 cell5=10; end;
              if gpixels(row-1,col)>=125 && cell3==0 cell3=10; end;
              if gpixels(row,col-1)>=125 && cell2==0 cell2=10; end;
              if gpixels(row,col+1)>=125 && cell4==0 cell4=10; end;
              
              if gpixels(row+1,col)>=100 && cell5==0 cell5=12; end;
              if gpixels(row-1,col)>=100 && cell3==0 cell3=12; end;
              if gpixels(row,col-1)>=100 && cell2==0 cell2=12; end;
              if gpixels(row,col+1)>=100 && cell4==0 cell4=12; end;

          count=0;
          %Don't turn a triangle pixel with three open borders into a triangle.
          if cell2==0 count=count+1;end;
          if cell3==0 count=count+1;end;
          if cell4==0 count=count+1;end;
          if cell5==0 count=count+1;end;
          %If a triangle pixel (c=1 or 6) is detected draw the correct triangle.
          if  count>1 && (gpixels(row,col)==1 || gpixels(row,col)==6)
              if cell5~=0 && cell2~=0 tri1=0;tri2=0;tri3=0;tri4=1;mid=cell2;end;
              if cell3~=0 && cell4~=0 tri1=0;tri2=0;tri3=1;tri4=0;mid=cell4;end;
              if cell4~=0 && cell5~=0 tri1=1;tri2=0;tri3=0;tri4=0;mid=cell4;end;
              if cell2~=0 && cell3~=0 tri1=0;tri2=1;tri3=0;tri4=0;mid=cell2;end;
          end;
          %Cell vertices are labeled 1-4 start lower left, go clockwise.
          %tri1 is tl(top left), tri2 is br(bottom right), tri3 is bl, tri4 is tr
          if tri1              
              simp(sid+1,:)=[sid 0 mat cell3 mid cell2 vertices(row,col,1) vertices(row,col,2) vertices(row,col,3)]; sid=sid+1;
          end;
          if tri2
              simp(sid+1,:)=[sid 0 mat cell4 cell5 mid vertices(row,col,1) vertices(row,col,3) vertices(row,col,4)]; sid=sid+1;
          end;
          if tri3
              simp(sid+1,:)=[sid 0 mat mid cell5 cell2 vertices(row,col,1) vertices(row,col,2) vertices(row,col,4)]; sid=sid+1;
          end;
          if tri4
              simp(sid+1,:)=[sid 0 mat cell4 mid cell3 vertices(row,col,2) vertices(row,col,3) vertices(row,col,4)]; sid=sid+1;
          end;
          end;   
     end;
end;

%Write mcin.m file.
vertices=vid;
simplices=sid;
file=fopen('mcin.m','w+');
fprintf(file,'mcsf_begin=1;\n\n');
fprintf(file,'      dim=2;         %% intrinsic manifold dimension\n');
fprintf(file,'    dimii=3;         %% imbedding manifold dimension\n');
fprintf(file,' vertices=%d;        %% number of vertices\n',vid);
fprintf(file,'simplices=%d;        %% number of simplices\n',sid);

fprintf(file,'vert=[\n');
fprintf(file,'%%-------- ---- ----------------- ----------------- -----------------\n');
fprintf(file,'%% Vert-ID Chrt X-Coordinate      Y-Coordinate      Z-Coordinate\n');
fprintf(file,'%%-------- ---- ----------------- ----------------- -----------------\n');
for x=1:vid
    fprintf(file,'%d %d %d %d %d\n',vert(x,1),vert(x,2),vert(x,3),vert(x,4),vert(x,5));
end;
fprintf(file,'];\n');
fprintf(file,'simp=[\n');
fprintf(file,'%%-------- ---- ---- ------------------- ---------------------------------------\n');
fprintf(file,'%% Simp-ID Grp  Mat  Face-Types          Vertex-Numbers\n');
fprintf(file,'%%-------- ---- ---- -------------------\n');
fprintf(file,'%%---------------------------------------\n');
for x=1:sid
    fprintf(file,'%d %d %d %d %d %d %d %d %d\n',simp(x,1),simp(x,2),simp(x,3),simp(x,4),simp(x,5),simp(x,6),simp(x,7),simp(x,8),simp(x,9));
end;
fprintf(file,'];\n');
fprintf(file,'mcsf_end=1;');
fclose(file);

%Draw the new mcsf mesh file immediately if parameter key1==1
try
if key1==1
   dim=2; dimii=3;
   [vert,simp]=read2(0);
   draw(vert,simp);
end;
end;

%Adds new vertex to variable vert and increments vid.
function wvert(a,b,c,d,e);  
    vert(a+1,:)=[a,b,c,d,e];
    vid=vid+1;
    if gpixels(m-d,c)~=255 vertices(m-d,c,4)=a;end;
    if gpixels(m-d+1,c)~=255 vertices(m-d+1,c,3)=a;end;
    if gpixels(m-d+1,c+1)~=255 vertices(m-d+1,c+1,2)=a;end;
    if gpixels(m-d,c+1)~=255 vertices(m-d,c+1,1)=a;end;
end;    

%Returns color of pixel in bitmap while handling errors.
function y=gpixels(r,c);  
    if r<1 y=255;
    elseif r>m y=255;
    elseif c<1 y=255;
    elseif c>n y=255;
    else y=pixels(r,c);
    end;
end;

%THE FOLLOWING 5 FUNCTIONS ARE FROM MICHAEL HOLST'S MCLITE PACKAGE.
%THEY ARE INCLUDED IN THIS FUNCTION TO ALLOW DISPLAY OF A MESH
%FILE WITHOUT MCLITE INSTALLED. (function draw has been modified slighty)

function draw(VERT,SIMP);
%DRAW  Draw the finite element mesh
%
% Usage: draw(VERT,SIMP);
%
% Author:   Michael Holst
% rcsid="$Id: img2mesh.m,v 1.1 2008/04/01 15:59:24 fetk Exp $"

%%% recover various problem dimensions

   [N,eight]     = size(VERT);
   [L,seventeen] = size(SIMP);
   T = 3;

%%% setup for plot

   clf;
   hold off;

%%% plot all vertices with dots

   plot(VERT(1:N,1),VERT(1:N,2), '.b') %inter
   axis equal;
   axis off;
   hold on;

%%% cycle through the simplices and draw the edges

   for ell = 1:L
      Xil(1:T,1:3) = VERT(SIMP(ell,1:T),1:3);
      plot([Xil(1:T,1)' Xil(1,1)], [Xil(1:T,2)' Xil(1,2)], 'k'); %inter

      %%% face info for the T face making up this simplex
      ftype = SIMP(ell,5:7);

      %%% plot the essential boundary vertices with circles
      for i = 1:T
         j=i+1;
         if (j>T) j=1; end;
         k=j+1;
         if (k>T) k=1; end;

         vtype = VERT(SIMP(ell,i),6);
         if ( vdiri(vtype) )
            plot(VERT(SIMP(ell,i),1), VERT(SIMP(ell,i),2), 'ok'); %Deut circ
         end
         if ( vbnd(ftype(i)) )
            Xil(1,1:3) = VERT(SIMP(ell,j),1:3);
            Xil(2,1:3) = VERT(SIMP(ell,k),1:3);
            if ( ftype(i)==1 )
               plot(Xil(1:2,1), Xil(1:2,2), 'k','linewidth',2); %Deut
            elseif ( ftype(i)==2 )
               plot(Xil(1:2,1), Xil(1:2,2), 'k'); %Neum
            elseif ( ftype(i)==4)
               plot(Xil(1:2,1), Xil(1:2,2), 's','color',[0 0 0]); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), 'w'); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), '--','color',[0 0 0],'linewidth',2); %Neum
            elseif ( ftype(i)==6)
               plot(Xil(1:2,1), Xil(1:2,2), 's','color',[0 0 0]); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), 'w'); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), '-.','color',[0 0 0],'linewidth',2); %Neum
            elseif ( ftype(i)==8)
               plot(Xil(1:2,1), Xil(1:2,2), 's','color',[0 0 0]); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), 'w'); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), ':','color',[0 0 0],'linewidth',2); %Neum
            elseif ( ftype(i)==10)
               plot(Xil(1:2,1), Xil(1:2,2), 's','color',[.35 .35 .35]); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), 'w'); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), ':','color',[.35 .35 .35],'linewidth',2); %Neum
            elseif ( ftype(i)==12)
               plot(Xil(1:2,1), Xil(1:2,2), 's','color',[.60 .60 .60]); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), 'w'); %Neum
               plot(Xil(1:2,1), Xil(1:2,2), ':','color',[.60 .60 .60],'linewidth',2); %Neum  
            end
         end
      end
   end
    
axis(axis + [-1 1 -1 1]);    
end

function [VERT_F,SIMP_F] = read2(key);
%READ  Define a 2D or 3D simplicial finite element (coordinates always in R^3)
%
% Usage: [VERT,SIMP] = read(n);
%
% Input:
%
%    An "MCSF" (MC-Simplex-Format) mesh definition file.
%
% Output:
%
%    (N = number of vertices in mesh)
%    (L = number of simplices in mesh)
%    (K = number of edges in mesh)
%
%    VERT = size(N,8) ==> vertices
%
%              VERT(:,1)  == x-coordinates of all vertices
%              VERT(:,2)  == y-coordinates of all vertices
%              VERT(:,3)  == z-coordinates of all vertices (3d case & manifold)
%              VERT(:,4)  == vertex identifier
%              VERT(:,5)  == vertex processor id or color
%              VERT(:,6)  == vertex type
%              VERT(:,7)  == blank; later used as first simplex ptr
%              VERT(:,8)  == blank; later used as first edge ptr
%
%    SIMP = size(L,17) ==> simplices
%
%              SIMP(:,1)  == 1st-vertex of the simplex
%              SIMP(:,2)  == 2nd-vertex of the simplex
%              SIMP(:,3)  == 3rd-vertex of the simplex
%              SIMP(:,4)  == 4th-vertex of the simplex (3d case)
%              SIMP(:,5)  == 1st-face type; face opposite vertex 1
%              SIMP(:,6)  == 2nd-face type; face opposite vertex 2
%              SIMP(:,7)  == 3rd-face type; face opposite vertex 3
%              SIMP(:,8)  == 4th-face type; face opposite vertex 4 (3d case)
%              SIMP(:,9)  == simplex identifier
%              SIMP(:,10) == simplex processor id or color
%              SIMP(:,11) == simplex type
%              SIMP(:,12) == blank; later used as simplex ring ptr 1
%              SIMP(:,13) == blank; later used as simplex ring ptr 2
%              SIMP(:,14) == blank; later used as simplex ring ptr 3
%              SIMP(:,15) == blank; later used as simplex ring ptr 4 (3d case)
%              SIMP(:,16) == first queue marker
%              SIMP(:,17) == second queue marker
%              SIMP(:,18) == marked marker
%              SIMP(:,19) == creation type (0=bisection,1=quad,2=quad interior)
%              SIMP(:,20) == simplex generation number
%
%    EDGE = size(K,8) ==> edges
%
%              EDGE(:,1) == 1st vertex making up edge
%              EDGE(:,2) == 2nd vertex making up edge
%              EDGE(:,3) == midpoint vertex number
%              EDGE(:,4) == edge identifier
%              EDGE(:,5) == edge processor id or color
%              EDGE(:,6) == edge type
%              EDGE(:,7) == next edge number in ring of edges about vertex 1
%              EDGE(:,8) == next edge number in ring of edges about vertex 2
%
% Author:   Michael Holst
% rcsid="$Id: img2mesh.m,v 1.1 2008/04/01 15:59:24 fetk Exp $"

%%% read in the mcf file

%    mcin

%%% this mapping allows us to grab the other two vertices given the third

    vmapOV3 = [ 2, 3;
                1, 3;
                1, 2 ];

%%% change to internal storage

    N=vertices;
    L=simplices;
    VERT_F=zeros(N,8);
    SIMP_F=zeros(L,20);
    ZERON=zeros(N,1);
    ZEROL=zeros(L,1);
    ONEN=ones(N,1);
    ONEL=ones(L,1);

    VERT_F(:,1:3) = vert(:,3:5);
    VERT_F(:,4:5) = vert(:,1:2);
    VERT_F(:,6)   = ZERON;
    VERT_F(:,4)   = VERT_F(:,4) + ONEN;

    SIMP_F(:,18)  = ONEL;
    SIMP_F(:,19)  = 2*ONEL;
    SIMP_F(:,20)  = ONEL;

    if (dim == 2)
        SIMP_F(:,1:3)  = simp(:,7:9);
        SIMP_F(:,5:7)  = simp(:,4:6);
        SIMP_F(:,9:11) = simp(:,1:3);
        SIMP_F(:,1)    = SIMP_F(:,1)  + ONEL;
        SIMP_F(:,2)    = SIMP_F(:,2)  + ONEL;
        SIMP_F(:,3)    = SIMP_F(:,3)  + ONEL;
        SIMP_F(:,9)    = SIMP_F(:,9)  + ONEL;
 
        %%% calculate the vertex types from the face types

        for element=1:L
            for vertex=1:3
                % vertex id and type (type is possibly incomplete)
                vid = SIMP_F(element,vertex);
                vtp = VERT_F(vid,6);

                % the two faces that use the vertex and their types
                f1   = vmapOV3(vertex,1);
                f2   = vmapOV3(vertex,2);
                f1tp = SIMP_F(element,4+f1);
                f2tp = SIMP_F(element,4+f2);

                % if anything touching vertex is diri, mark vertex as such
                if ( vdiri(vtp) | vdiri(f1tp) | vdiri(f2tp) )
                    tmp = max( vdiri(f1tp)*f1tp, vdiri(f2tp)*f2tp );
                    VERT_F(vid,6) = max( vdiri(vtp)*vtp, tmp );

                % else if anything touching vertex is neum, mark vertex as such
                else
                    tmp = max( vneum(f1tp)*f1tp, vneum(f2tp)*f2tp );
                    VERT_F(vid,6) = max( vneum(vtp)*vtp, tmp );
                end
            end
        end

    else
        SIMP_F(:,1:4)  = simp(:,8:11);
        SIMP_F(:,5:8)  = simp(:,4:7);
        SIMP_F(:,9:11) = simp(:,1:3);
        SIMP_F(:,1)    = SIMP_F(:,1)  + ONEL;
        SIMP_F(:,2)    = SIMP_F(:,2)  + ONEL;
        SIMP_F(:,3)    = SIMP_F(:,3)  + ONEL;
        SIMP_F(:,4)    = SIMP_F(:,4)  + ONEL;
        SIMP_F(:,9)    = SIMP_F(:,9)  + ONEL;

        %%% WARNING: we don't calculate vertex types in 3D case
        assertw( (0), 'read: refusing to calulate 3D vertex types' );

    end
end;

function [result] = vdiri(type);
%VDIRI  Yes or no answer to whether triang/edge/vertex is "dirichlet" type
%
% Usage: [result] = vdiri(type);
%
% Input:
%
%    type   = the type in question
% 
% Output:
%    
%    result = yes or no
%
% Author:   Michael Holst
% rcsid="$Id: img2mesh.m,v 1.1 2008/04/01 15:59:24 fetk Exp $"

   [m,n] = size(type);

   odd = type ./ ( 2 * ones(m,n) );
   oddint = floor(odd);
   oddtst = abs(odd - oddint);

   result = (oddtst~=0);

end;

function [result] = vneum(type);
%VNEUM  Yes or no answer to whether triang/edge/vertex is "neumann" type
%
% Usage: [result] = vneum(type);
%
% Input:
%
%    type   = the type in question
% 
% Output:
%    
%    result = yes or no
%
% Author:   Michael Holst
% rcsid="$Id: img2mesh.m,v 1.1 2008/04/01 15:59:24 fetk Exp $"

   [m,n] = size(type);

   even = type ./ ( 2 * ones(m,n) );
   evenint = floor(even);
   eventst = ( (type > 0) & (~abs(even - evenint)) );

   result = (eventst~=0);

end;

function [result] = vbnd(type);
%VBND  Yes or no answer to whether triang/edge/vertex is "boundary" type
%
% Usage: [result] = vbnd(type);
%
% Input:
%
%    type   = the type in question
% 
% Output:
%    
%    result = yes or no
%
% Author:   Michael Holst
% rcsid="$Id: img2mesh.m,v 1.1 2008/04/01 15:59:24 fetk Exp $"

   [m,n] = size(type);

   result = (type~=0);

end;

end

           
                      