% This function calls a C++ program to generate solution to a PDE, then plots it.
% N: the polynomial order of the expansion basis.
% n: the amount of evaluation points on an equidistant grid for the solution output

function retval = threeElmt(n,N)
	clc
	close all
	
	%% Call C++ program to generate fields and quadPoints
	system(["./run " , num2str(n) , " " , num2str(N)]);
	
	figure
	hold on
	plotField("","grid.txt");	
endfunction