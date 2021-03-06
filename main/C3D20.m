function [n,w,xi,N,dNdxi]=C3D20
%====================== No. integration points =============================
%
%   Defines the number of integration points:be used for
%   each element type
%

n = 27;
ncoord=3;  
nodes=20;
%
%====================== INTEGRATION POINTS ==================================
%
%   Defines positions of integration points
%
         xi = zeros(n,ncoord);
          x1D = [-0.7745966692,0.,0.7745966692];
           for k = 1:3
           for j = 1:3 
            for i = 1:3
              p = 9*(k-1) + 3*(j-1) + i;
               xi(p,1) = x1D(i);
               xi(p,2) = x1D(j);
               xi(p,3) = x1D(k);
             end
           end
         end
   w1D = [0.555555555,0.888888888,0.55555555555];
   w = zeros(1,n);
         for k = 1:3
           for j = 1:3
             for i = 1:3
               n = 9*(k-1)+3*(j-1)+i;
               w(n) = w1D(i)*w1D(j)*w1D(k);
             end
           end    
         end

%
%================= SHAPE FUNCTIONS ==================================
%
%        Nij: Shape functions of the Int Point i [4x4] Ni [4x1]

N=zeros(n,nodes);
for i1=1:n
        N(i1,1) = (1.-xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3))*(-xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)/8.;
       N(i1,2) = (1.+xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3))*(xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)/8.;
       N(i1,3) = (1.+xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3))*(xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)/8.;
       N(i1,4) = (1.-xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3))*(-xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)/8.;
       N(i1,5) = (1.-xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3))*(-xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)/8.;
       N(i1,6) = (1.+xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3))*(xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)/8.;
       N(i1,7) = (1.+xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3))*(xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)/8.;
       N(i1,8) = (1.-xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3))*(-xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)/8.;
       N(i1,9)  = (1.-xi(i1,1)^2)*(1.-xi(i1,2))*(1.-xi(i1,3))/4.;
       N(i1,10) = (1.+xi(i1,1))*(1.-xi(i1,2)^2)*(1.-xi(i1,3))/4.;
       N(i1,11) = (1.-xi(i1,1)^2)*(1.+xi(i1,2))*(1.-xi(i1,3))/4.;
       N(i1,12) = (1.-xi(i1,1))*(1.-xi(i1,2)^2)*(1.-xi(i1,3))/4.;
       N(i1,13) = (1.-xi(i1,1)^2)*(1.-xi(i1,2))*(1.+xi(i1,3))/4.;
       N(i1,14) = (1.+xi(i1,1))*(1.-xi(i1,2)^2)*(1.+xi(i1,3))/4.;
       N(i1,15) = (1.-xi(i1,1)^2)*(1.+xi(i1,2))*(1.+xi(i1,3))/4.;
       N(i1,16) = (1.-xi(i1,1))*(1.-xi(i1,2)^2)*(1.+xi(i1,3))/4.;
       N(i1,17) = (1.-xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)^2)/4.;
       N(i1,18) = (1.+xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)^2)/4.;
       N(i1,19) = (1.+xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)^2)/4.;
       N(i1,20) = (1.-xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)^2)/4.;  
end
%
%================= SHAPE FUNCTION DERIVATIVES ======================
%
%        Nij,r: Dev of shape functions of the Int Point i [4x8]
%        [2*i-1 2*i] => dNi,r [4x2]
dNdxi = zeros(n*ncoord,nodes);
for i1=1:n
       dNdxi(3*i1-2,1) = (-(1.-xi(i1,2))*(1.-xi(i1,3))*(-xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)-(1.-xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1-1,1) = (-(1.-xi(i1,1))*(1.-xi(i1,3))*(-xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)-(1.-xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1,1) = (-(1.-xi(i1,1))*(1.-xi(i1,2))*(-xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)-(1.-xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)))/8.;

       dNdxi(3*i1-2,2) = ((1.-xi(i1,2))*(1.-xi(i1,3))*(xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)+(1.+xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1-1,2) = (-(1.+xi(i1,1))*(1.-xi(i1,3))*(xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)-(1.+xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1,2) = (-(1.+xi(i1,1))*(1.-xi(i1,2))*(xi(i1,1)-xi(i1,2)-xi(i1,3)-2.)-(1.+xi(i1,1))*(1.-xi(i1,2))*(1.-xi(i1,3)))/8.;

       dNdxi(3*i1-2,3) = ((1.+xi(i1,2))*(1.-xi(i1,3))*(xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)+(1.+xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1-1,3) = ((1.+xi(i1,1))*(1.-xi(i1,3))*(xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)+(1.+xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1,3) = (-(1.+xi(i1,1))*(1.+xi(i1,2))*(xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)-(1.+xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)))/8.;

       dNdxi(3*i1-2,4) = (-(1.+xi(i1,2))*(1.-xi(i1,3))*(-xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)-(1.-xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1-1,4) = ((1.-xi(i1,1))*(1.-xi(i1,3))*(-xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)+(1.-xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1,4) = (-(1.-xi(i1,1))*(1.+xi(i1,2))*(-xi(i1,1)+xi(i1,2)-xi(i1,3)-2.)-(1.-xi(i1,1))*(1.+xi(i1,2))*(1.-xi(i1,3)))/8.;
       dNdxi(3*i1-2,5) = (-(1.-xi(i1,2))*(1.+xi(i1,3))*(-xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)-(1.-xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-1,5) = (-(1.-xi(i1,1))*(1.+xi(i1,3))*(-xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)-(1.-xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1,5) = ((1.-xi(i1,1))*(1.-xi(i1,2))*(-xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)+(1.-xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-2,6) = ((1.-xi(i1,2))*(1.+xi(i1,3))*(xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)+(1.+xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-1,6) = (-(1.+xi(i1,1))*(1.+xi(i1,3))*(xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)-(1.+xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1,6) = ((1.+xi(i1,1))*(1.-xi(i1,2))*(xi(i1,1)-xi(i1,2)+xi(i1,3)-2.)+(1.+xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-2,7) = ((1.+xi(i1,2))*(1.+xi(i1,3))*(xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)+(1.+xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-1,7) = ((1.+xi(i1,1))*(1.+xi(i1,3))*(xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)+(1.+xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1,7) = ((1.+xi(i1,1))*(1.+xi(i1,2))*(xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)+(1.+xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-2,8) = (-(1.+xi(i1,2))*(1.+xi(i1,3))*(-xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)-(1.-xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-1,8) = ((1.-xi(i1,1))*(1.+xi(i1,3))*(-xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)+(1.-xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1,8) = ((1.-xi(i1,1))*(1.+xi(i1,2))*(-xi(i1,1)+xi(i1,2)+xi(i1,3)-2.)+(1.-xi(i1,1))*(1.+xi(i1,2))*(1.+xi(i1,3)))/8.;
       dNdxi(3*i1-2,9)  = -2.*xi(i1,1)*(1.-xi(i1,2))*(1.-xi(i1,3))/4.;
       dNdxi(3*i1-1,9)  = -(1.-xi(i1,1)^2)*(1.-xi(i1,3))/4.;
       dNdxi(3*i1,9)  = -(1.-xi(i1,1)^2)*(1.-xi(i1,2))/4.;
       dNdxi(3*i1-2,10)  = (1.-xi(i1,2)^2)*(1.-xi(i1,3))/4.;
       dNdxi(3*i1-1,10)  = -2.*xi(i1,2)*(1.+xi(i1,1))*(1.-xi(i1,3))/4.;
       dNdxi(3*i1,10)  = -(1.-xi(i1,2)^2)*(1.+xi(i1,1))/4.;
       dNdxi(3*i1-2,11)  = -2.*xi(i1,1)*(1.+xi(i1,2))*(1.-xi(i1,3))/4.;
       dNdxi(3*i1-1,11)  = (1.-xi(i1,1)^2)*(1.-xi(i1,3))/4.;
       dNdxi(3*i1,11)  = -(1.-xi(i1,1)^2)*(1.+xi(i1,2))/4.;
       dNdxi(3*i1-2,12)  = -(1.-xi(i1,2)^2)*(1.-xi(i1,3))/4.;
       dNdxi(3*i1-1,12)  = -2.*xi(i1,2)*(1.-xi(i1,1))*(1.-xi(i1,3))/4.;
       dNdxi(3*i1,12)  = -(1.-xi(i1,2)^2)*(1.-xi(i1,1))/4.;
       dNdxi(3*i1-2,13)  = -2.*xi(i1,1)*(1.-xi(i1,2))*(1.+xi(i1,3))/4.;
       dNdxi(3*i1-1,13)  = -(1.-xi(i1,1)^2)*(1.+xi(i1,3))/4.;
       dNdxi(3*i1,13)  = (1.-xi(i1,1)^2)*(1.-xi(i1,2))/4.;
       dNdxi(3*i1-2,14)  = (1.-xi(i1,2)^2)*(1.+xi(i1,3))/4.;
       dNdxi(3*i1-1,14)  = -2.*xi(i1,2)*(1.+xi(i1,1))*(1.+xi(i1,3))/4.;
       dNdxi(3*i1,14)  = (1.-xi(i1,2)^2)*(1.+xi(i1,1))/4.;
       dNdxi(3*i1-2,15)  = -2.*xi(i1,1)*(1.+xi(i1,2))*(1.+xi(i1,3))/4.;
       dNdxi(3*i1-1,15)  =  (1.-xi(i1,1)^2)*(1.+xi(i1,3))/4.;
       dNdxi(3*i1,15)  = (1.-xi(i1,1)^2)*(1.+xi(i1,2))/4.;
       dNdxi(3*i1-2,16)  = -(1.-xi(i1,2)^2)*(1.+xi(i1,3))/4.;
       dNdxi(3*i1-1,16)  = -2.*xi(i1,2)*(1.-xi(i1,1))*(1.+xi(i1,3))/4.;
       dNdxi(3*i1,16)  = (1.-xi(i1,2)^2)*(1.-xi(i1,1))/4.;
       dNdxi(3*i1-2,17) = -(1.-xi(i1,2))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1-1,17) = -(1.-xi(i1,1))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1,17) = -xi(i1,3)*(1.-xi(i1,1))*(1.-xi(i1,2))/2.;
       dNdxi(3*i1-2,18) = (1.-xi(i1,2))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1-1,18) = -(1.+xi(i1,1))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1,18) = -xi(i1,3)*(1.+xi(i1,1))*(1.-xi(i1,2))/2.;
       dNdxi(3*i1-2,19) = (1.+xi(i1,2))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1-1,19) = (1.+xi(i1,1))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1,19) = -xi(i1,3)*(1.+xi(i1,1))*(1.+xi(i1,2))/2.;
       dNdxi(3*i1-2,20) = -(1.+xi(i1,2))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1-1,20) = (1.-xi(i1,1))*(1.-xi(i1,3)^2)/4.;
       dNdxi(3*i1,20) = -xi(i1,3)*(1.-xi(i1,1))*(1.+xi(i1,2))/2.;
       
end
end
%
