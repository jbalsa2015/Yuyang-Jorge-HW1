function [n,w,xi,N,dNdxi]=C2D6
%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%

n = 3;
ncoord = 2;  
nodes = 6;
%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
         xi = zeros(n,2);
         xi(1,1) = 0.6;
         xi(1,2) = 0.2;
         xi(2,1) = 0.2;
         xi(2,2) = 0.6;
         xi(3,1) = 0.2;
         xi(3,2) = 0.2;

         w =[1./6., 1./6.,1./6.];
         
%
%================= SHAPE FUNCTIONS ==================================
%
%        Nij: Shape functions of the Int Point i [4x4] Ni [4x1]

N=zeros(n,nodes);
 for i1=1:n
       
     xi3 = 1.-xi(i1,1)-xi(i1,2);
      N(i1,1) = (2.*xi(i1,1)-1.)*xi(i1,1);
      N(i1,2) = (2.*xi(i1,2)-1.)*xi(i1,2);
      N(i1,3) = (2.*xi3-1.)*xi3;
      N(i1,4) = 4.*xi(i1,1)*xi(i1,2);
      N(i1,5) = 4.*xi(i1,2)*xi3;
      N(i1,6) = 4.*xi3*xi(i1,1);
       
 end
%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
%        Nij,r: Dev of shape functions of the Int Point i [4x8]
%        [2*i-1 2*i] => dNi,r [4x2]
dNdxi = zeros(ncoord*n,nodes);
for i1=1:n
%        dNdxi(i1*2-1,1) = 1.;
%        dNdxi(i1*2,2) = 1.;
%        dNdxi(i1*2-1,3) = -1.;
%        dNdxi(i1*2,3) = -1.;

         xi3 = 1.-xi(i1,1)-xi(i1,2);
       dNdxi(i1*2-1,1) = 4.*xi(i1,1)-1.;
       dNdxi(i1*2,2) = 4.*xi(i1,2)-1.;
       dNdxi(i1*2-1,3) = -(4.*xi3-1.);
       dNdxi(i1*2,3) = -(4.*xi3-1.);
       dNdxi(i1*2-1,4) = 4.*xi(i1,2);
       dNdxi(i1*2,4) = 4.*xi(i1,1);
       dNdxi(i1*2-1,5) = -4.*xi(i1,2);
       dNdxi(i1*2,5) = -4.*(-1.+xi(i1,1)+2.*xi(i1,2));
       dNdxi(i1*2-1,6) = -4.*(-1.+2.*xi(i1,1)+xi(i1,2));
       dNdxi(i1*2,6) = -4.*xi(i1,1);

end
end
%
