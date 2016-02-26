function [n,w,xi,N,dNdxi]= C3D4
%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%

n = 1;
ncoord = 3;  
nodes = 4;
%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
         xi = zeros(n,ncoord);
         xi(1,1) = 0.25;
         xi(1,2) = 0.25;
         xi(1,3) = 0.25;
         w =[1./6];
         
%
%================= SHAPE FUNCTIONS ==================================
%
%        Nij: Shape functions of the Int Point i [4x4] Ni [4x1]

N=zeros(n,nodes);
% for i1=1:n
       N(1,1) =   xi(1,1);
       N(1,2) =   xi(1,2);
       N(1,3) =   xi(1,3);
       N(1,4) = 1.-xi(1,1)-xi(1,2)-xi(1,3);
% end
%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
%        Nij,r: Dev of shape functions of the Int Point i [4x8]
%        [2*i-1 2*i] => dNi,r [4x2]
dNdxi = zeros(ncoord*n,nodes);
%for i1=1:n
       dNdxi(1,1) = 1.;
       dNdxi(2,2) = 1.;
       dNdxi(3,3) = 1.;
       dNdxi(1,4) = -1.;
       dNdxi(2,4) = -1.;
       dNdxi(3,4) = -1.;
         
%end
end
%
