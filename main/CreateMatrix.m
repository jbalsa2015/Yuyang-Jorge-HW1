function [K,f] = CreateMatrix(X,T,pospg,pespg,N,dNdxi,c) 
% Stiffness matrix K and r.h.s vector f 
% obtained by discretizing a heat equation
% 
% X:            nodal coordinates
% T:            connectivities (elements)
% pospg, pespg: Gauss points and weights in the reference element
% N,Nxi,Neta:   shape functions on the Gauss points 


global diffusion  h  

nu = diffusion;

% Total number of elements and number of nodes in aech one
[numel,nen] = size(T); 
% Total number of nodes
numnp = size(X,1); 
 
% Allocate storage
K = zeros(numnp,numnp); 
f = zeros(numnp,1); 

% Loop on elements
for ielem = 1:numel
    % Te: global number of the nodes in the current element
    % Xe: coordenates of the nodes in the current element
    Te = T(ielem,:); 
    % Get local information
    Xe = X(Te,:); 
    % Element matrices
    if c==1
        [Ke,fe] = MatE1(Xe,nen,pospg,pespg,N,dNdxi); 
    elseif c==2
        [Ke,fe] = MatE2(Xe,nen,pospg,pespg,N,dNdxi); 
    elseif c==3
        [Ke,fe] = MatE3(Xe,nen,pospg,pespg,N,dNdxi);
    elseif c==4
        [Ke,fe] = MatE4(Xe,nen,pospg,pespg,N,dNdxi);
    elseif c==5
        [Ke,fe] = MatE5(Xe,nen,pospg,pespg,N,dNdxi);
    elseif c==6
        [Ke,fe] = MatE6(Xe,nen,pospg,pespg,N,dNdxi);
    elseif c==7
        [Ke,fe] = MatE7(Xe,nen,pospg,pespg,N,dNdxi);
    elseif c==8
        [Ke,fe] = MatE8(Xe,nen,pospg,pespg,N,dNdxi);
    end
   
    % Assemble the element matrices
    K(Te,Te) = K(Te,Te) + Ke; 
    f(Te) = f(Te) + fe; 
    clear Ke; clear fe; 
end 
 
 
