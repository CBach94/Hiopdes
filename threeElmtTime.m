% This function calls a C++ program to generate solution to a PDE, then plots it.
% N: the polynomial order of the expansion basis.
% n: the amount of evaluation points on an equidistant grid for the solution output

function retval = threeElmtTime(n,N)
	close all
	figure
	hold on
	
	system(["./cRun " , num2str(n) , " " , num2str(N)]);

	T = 10
	
	for t = 0:T-1
		field = ["_" , num2str(t)]
		plotField(field,"x.txt")
		pause(0.1)
	end
	
	
	
	
endfunction