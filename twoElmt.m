clear all
close all

system("./run 30 30")

psi = load("psi2.txt");
x 	= load("x.txt");

n = length(x);
N = n^2;
[X,Y] = meshgrid(x,x);
XX = [X;X];
YY = [Y;Y+2];

psiT = psi(1:N);
psiB = psi(N+1:end);
psiT = reshape(psiT,n,n);
psiB = reshape(psiB,n,n);

PSI = [psiB;psiT];
disp("Plotting")
close all

figure
hold on

contourf(XX,YY,PSI)
plot(X,Y+2,'r.')
plot(X,Y,'c.')

axis tight
colorbar
xlabel "x"
ylabel "y"