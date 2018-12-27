function retval = threeElmt2(N,n)
	close all

system(["./run " , num2str(N) , " " , num2str(n)]);

psi = load("psi2.txt");
x 	= load("x.txt");

n = length(x);
N = n^2;
[X,Y] = meshgrid(x,x);
XX = [X;X;X+2];
YY = [Y;Y+2;Y+2];

psiT 	= psi(1:N);
psiB 	= psi(N+1-n:2*N-n);
psiTR = psi(2*N+1-2*n:3*N-2*n);

psiT  = reshape(psiT,n,n);
psiTR = reshape(psiTR,n,n);
psiB  = reshape(psiB,n,n);

PSI = [psiB;psiT;psiTR];

disp("Plotting")
close all

figure
hold on
contours = 20;
contourf(X,Y,psiB,contours)
contourf(X,Y+2,psiT,contours)
contourf(X+2,Y+2,psiTR,contours)

%surf(X,Y,psiB)
%surf(X,Y+2,psiT)
%surf(X+2,Y+2,psiTR)

plot(X,Y+2,'r+')
plot(X,Y,'c+')
plot(X+2,Y+2,'m+')
title("h=3, p=30, Advection-Diffusion")
axis tight
colorbar
colormap("hot")
shading interp
xlabel "x"
ylabel "y"
endfunction