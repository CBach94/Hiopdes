function [L,dLdx] = LegendreVal(x,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [L,dLdx] = LegendreVal(x,N)
% This function can be used to calculate the values of the Legendre
% polynomial and its derivative.
% 
% Written by Jasper Kreeft - 2009
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check if it is a row vector, and if not transpose it
if size(x,1)>size(x,2)
    x=x';
end

% Number of values
nx = size(x,2);

% Initiate
L      = zeros(N+1,nx);
% Set first row in matrix as L0=1
L(1,:) = ones(1,nx);

% If N>0 set second row in matrix as L1=x
if N>0; L(2,:) = x; end

% Using Bonnets recursion formula to calculate the next rows
for k = 2:N
    L(k+1,:) = (2*k-1)/k.*x.*L(k,:)-(k-1)/k.*L(k-1,:);
end

% Check if derivative is assigned in the output
if nargout==2

% Initiate
dLdx = zeros(N+1,nx);

% Set first row
dLdx(1,:) = zeros(1,nx);

% Check if N>0 and set second row in matrix as L1'=1
if N>0; dLdx(2,:) = ones(1,nx); end

% Using Bonnets recursion formula to calculate next rows
for k = 2:N
    dLdx(k+1,:) = dLdx(k-1,:)+(2*k-1)*L(k,:);
end

% Send last row to second output
dLdx = dLdx(N+1,:);
end

% Send Last row of Legendre polynomials to first output
L  =  L(N+1,:);