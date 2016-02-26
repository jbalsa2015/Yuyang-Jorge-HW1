% This program solves a convection-diffusion problem 
% in a square domain [0,1]x[0,1] using bilinear elements.
% 
clc;
clear;

global diffusion 
diffusion = 1;

% NUMERICAL INTEGRATION
% Quadrature,Shape Function




%ttim=0; %initialize time counter

a = input('Please,input the coordinate file name:','s');
b = input('Please,input the connectivity matrix file name:','s');
c = input('Please,input the Inlet boundary condition file name:','s');
d = input('Please,input the Outlet boundary condtion file name:','s');


tic;%tic1
t1=clock;

X=load(a);
T=load(b);
% NUMERICAL INTEGRATION
% Quadrature,Shape Functions
NumberofElementNodes=size(T,2)-1;
switch NumberofElementNodes
    case 4 
        if size(X,2)<= 3
        X = X(:,2:3);
        T = T(:,2:5);
        [n,wpg,pospg,N,dNdxi] = C2D4 ;
        else
      
        X = X(:,2:4);
        T = T(:,2:5);
  
        [n,wpg,pospg,N,dNdxi] = C3D4 ;
        end 
   case 8
       if size(X,2)<=3 %using if control to swith between 2d and 3d element
       X = X(:,2:3);
       T = T(:,2:9);
       [n,wpg,pospg,N,dNdxi] = C2D8 ;
       
       else
        
        X = X(:,2:4);
        T = T(:,2:9);
        [n,wpg,pospg,N,dNdxi] = C3D8 ;
       end
   case 3
       
       X = X(:,2:3);
       T = T(:,2:4);
       [n,wpg,pospg,N,dNdxi] = C2D3 ;
       
    case 6
        
        X = X(:,2:3);
        T = T(:,2:7);
         % ngaus = 1
      
        [n,wpg,pospg,N,dNdxi] = C2D6 ;
  
    case 20
        
        X = X(:,2:4);
        T = T(:,2:21);
    
        [n,wpg,pospg,N,dNdxi] = C3D20 ;
 
    case 10
       
        X = X(:,2:4);
        T = T(:,2:11);
        [n,wpg,pospg,N,dNdxi] = C3D10 ;
        
    otherwise
            disp('No valid input. Execute again.')
end

nnode = length(X);
nelem = length(T);
% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
%[K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);
% BOUNDARY CONDITIONS 
% nodesDir1: nodes on which solution is u=1
% nodesDir0: nodes on which solution is u=0
     
switch NumberofElementNodes
    case 4
        if size(X,2)<=2;
        [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,1);
        Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
        else
        [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,7);
        Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
        end
    case 8
        if size(X,2)<=2
         [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,2);
         Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
        else
        [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,5);
         Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
        end
    case 3
         [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,3);   
         Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
    case 6
        [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,4);
          Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
   
    case 20
        [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,6);
          Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
   
    case 10
        [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi,8);
        Inlet=load(c);
        Outlet=load(d);
        nodesDir1 = Inlet';
        nodesDir0 = Outlet';
end
    % Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];

ndir = size(C,1);
neq  = size(f,1);
A = zeros(ndir,neq);
A(:,C(:,1)) = eye(ndir);
b = C(:,2);


% SOLUTION OF THE LINEAR SYSTEM
% Entire matrix
Ktot = [K A';A zeros(ndir,ndir)];
ftot = [f;b];

sol = Ktot\ftot;
Temp = sol(1:neq);
multip = sol(neq+1:end);
% printing heading to file
switch NumberofElementNodes
    case 4
        if size(X,2)<=2
        f=fopen('2D Quadrilateral_Linear.vtk','w');
        else
        f=fopen('3d_tetrahedra_lin.vtk','w');
        end
    case 8
        if size(X,2)<=2
        f=fopen('2d_Quadrilateral_Quad.vtk','w');
        else
        f=fopen('3d_hexahedra_lin.vtk','w');
        end
    case 3
        f=fopen('2d_triangle_lin.vtk','w');
    case 6
        f=fopen('2d_triangle_quad.vtk','w');
   
    case 20
        f=fopen('3d_hexahedra_quad.vtk','w');
    
    case 10
        f=fopen('3d_tetrahedra_quad.vtk','w');
end
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%

	% printing coordinates
 if size(X,2)<=2
               for i1=1:nnode
                fprintf(f,'%14.8E %14.8E %14.8E \n',X(i1,1:2),0);
               end
 else          
for i1=1:nnode
   fprintf(f,'%14.8E %14.8E %14.8E \n',X(i1,1:3));
end
 end
fprintf(f,'\n');

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%
switch NumberofElementNodes
    case 4
       if size(X,2)<=2
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
              fprintf(f,' %4i ', 9);
        end
      else 
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 10);
        end
       end
    case 8
      if size(X,2)<=2
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*9);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i %10i %10i %10i %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ',23);
        end
     else
     
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*9);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i %10i  %10i %10i %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 12);
        end
     end
    case 3
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*4);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i %10i %10i %10i %10i\n',3,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 5);
        end
    case 6
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*7);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i %10i %10i %10i %10i\n',6,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 22);
        end
    
        
    case 20
        
    fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*21);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i %10i  %10i %10i %10i %10i  %10i %10i %10i %10i  %10i %10i %10i %10i  %10i %10i %10i\n',20,T(i1,1:20)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 25);
        end    
        

    case 10
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*11);
        for i1=1:nelem     
            fprintf(f,'%4i  %10i  %10i %10i %10i %10i %10i %10i %10i %10i %10i\n',10,T(i1,1:10)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 24);
        end 
end
fprintf(f,'\n');          

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
           fprintf(f,'%14.8E\n', Temp(i1) );
end

fclose('all');

disp(['total time for the code:',num2str(etime(clock,t1))]);
  

