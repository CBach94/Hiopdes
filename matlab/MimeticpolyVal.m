function [h,dhdx] = MimeticpolyVal(x,N,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [h,e] = MimeticpolyVal(x,N,k)
%
% h: Lagrange polynomials
% e: Edge polynomials
%
% x: row vector with all evaluation points
% N: number of cells in element
% k=1: GLL (Gauss-Lobatto-Legendre)
% k=2: G   (Gauss)
% k=3: EG  (Extended Gauss)
%
% For example:
% N=3 |--x--|--x--|--x--|
% has 3 cell, 4 GLL and 3 G nodes, that gives
% 3rd order GLL and 2nd order G Lagrange polynomials
% and 2nd order GLL and 1st order G Edge polynomials
% 
% Written by Jasper Kreeft - 2010
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Testscript
% clear all; close all; clc;
% N = 3;
% x = linspace(-1,1,100);%sort([-1; Gnodes(N); 1]);%
% k=1;
% nargout=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if it is a row vector, and if not transpose it
if size(x,1)>size(x,2)
    x=x';
end

% Number of values
nx = size(x,2);

% Call the interpolation nodes from the functions. Extended Gauss is just
% Gauss with end-points -1 and +1.
if k==1
    [nodes, w] = GLLnodes(N);
    n = N+1;
elseif k==2
    nodes = Gnodes(N);
    n = N;
elseif k==3
    nodes = [-1 Gnodes(N) 1];
    n = N+2;
end

% Calculate Legendre polynomial values of the evaluation points
[L,dLdx]     = LegendreVal(x,N);
% Calculate Legendre polynomial values of interpolation nodes
[L_z,dLdx_z] = LegendreVal(nodes,N);

% Initiate Lagrange polynomials
h = zeros(n,length(x));

% Write the Lagrange polynomial values at evaluation points. Using that the
% nodes are a polynomial with roots at the nodal points, so that h =
% f(x)/(f'(x_i)*(x-x_i))
if k==1
    for i=1:N+1
        h(i,:) = (x.^2-1).*dLdx./(N*(N+1)*L_z(i).*(x-nodes(i)));
    end
elseif k==2
    for i=1:N
        h(i,:) = L./(dLdx_z(i)*(x-nodes(i)));
    end
elseif k==3
    for i=1:N+2
        h(i,:) = ((x.^2-1).*L)./((2*nodes(i)*L_z(i)+(nodes(i)^2-1)*dLdx_z(i))*(x-nodes(i)));
    end
end
for i=1:n
    h(i,roundn(x,-10)==roundn(nodes(i),-10)) = 1;
end
h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
% figure
% hold on
% for i=1:n
%     plot(x,h(i,:),kleur(i));
% end
% grid
% xlabel('\xi')
% ylabel('h_i(\xi)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check is second output is assigned (edge polynomials)
if nargout>=2

% Initiate
dhdx = zeros(n,length(x));

if k==1
    % Calculate the derivative
    for i=1:N+1
        dhdx(i,:) = ( N*(N+1)*L.*(x-nodes(i))-(x.^2-1).*dLdx )./( N*(N+1)*L_z(i)*((x-nodes(i)).^2));
    end
    % Check if nodes coincide evaluation point
    for j=1:N+1
        index = roundn(x,-10)==roundn(nodes(j),-10);
        for i=1:N+1
            dhdx(i,index) = L(index)./(L_z(i)*(x(index)-nodes(i)));
        end
        dhdx(j,index) = 0;
    end
    % Check if it is an end point
    if x(1)==nodes(1)
        dhdx(1,1) = -N*(N+1)/4;
    end
    if x(nx)==nodes(N+1)
        dhdx(N+1,nx) = N*(N+1)/4;
    end

elseif k==2
    % Calculate the derivative
    for i=1:N
        dhdx(i,:) = (dLdx.*(x-nodes(i))-L)./(dLdx_z(i)*(x-nodes(i)).^2);
    end
    % Check if nodes coincide evaluation point
    for j=1:N
        index = roundn(x,-10)==roundn(nodes(j),-10);
        for i=1:N
            dhdx(i,index) = dLdx(index)./(dLdx_z(i)*(x(index)-nodes(i)));
        end
        dhdx(j,index) = nodes(j)/(1-nodes(j)^2);
    end


elseif k==3
    % Calculate the derivative
    for i=1:N+2
        dhdx(i,:) = ( (2*x.*L+(x.^2-1).*dLdx).*(x-nodes(i))-(x.^2-1).*L ) ./ ( (2*nodes(i)*L_z(i)+(nodes(i)^2-1)*dLdx_z(i))*(x-nodes(i)).^2 );
    end

    for j=1:N+2
        % Check if nodes coincide evaluation point
        index = roundn(x,-10)==roundn(nodes(j),-10);
        if max(index)==1
            for i=1:N+1
                dhdx(i,index) = ( 2*x(index)*L(index)+(x(index)^2-1)*dLdx(index) )/( (2*nodes(i)*L_z(i)+(nodes(i)^2-1)*dLdx_z(i))*(x(index)-nodes(i)) );
            end
        end
        dhdx(j,index) = -nodes(j)/(1-nodes(j)^2);
    end
    % Check if it is an end point
    if x(1)==nodes(1)
        dhdx(1,1) = -(N*(N+1)+1)/2;
    end
    if x(nx)==nodes(N+2)
        dhdx(N+2,nx) = (N*(N+1)+1)/2;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
% figure
% hold on
% for i=1:n
%     plot(x,dhdx(i,:),kleur(i));
% end
% grid
% xlabel('\xi')
% ylabel('dh_i/d\xi')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate end polynomial (sum of the differentiated Lagrange polynomial)
e = zeros(n-1,nx);
for i=1:n-1
    for j=1:nx
        e(i,j) = -sum(dhdx(1:i,j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kleur ='brgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmckybrgmcky';
% figure
% hold on
% for i=1:n-1
%     plot(x,e(i,:),kleur(i));
% end
% grid
% xlabel('\xi')
% ylabel('e_i(\xi)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end