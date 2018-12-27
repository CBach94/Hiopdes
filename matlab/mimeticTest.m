close all
clear all
clc


n = 10;
[gll,w] = GLLnodes(n);
[h,dh] = MimeticpolyVal(gll,n,1);

D = dh';

y = sin(pi*gll)';

dy = D*y;

figure
plot(gll,dy)