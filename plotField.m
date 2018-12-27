function plotField(field,grid)
	x 	= load(grid);

	% Reshape fields
	n = length(x);
	N = n^2;
	[X,Y] = meshgrid(x,x);
	XX = [X;X;X+2];
	YY = [Y;Y+2;Y+2];

	% Partition fields into regions	
	psiB = load(["psiB",field,".txt"]);
	psiT = load(["psiT",field,".txt"]);
	psiTR = load(["psiTR",field,".txt"]);
	
	psiT  = reshape(psiT,n,n);
	psiTR = reshape(psiTR,n,n);
	psiB  = reshape(psiB,n,n);

	%% Plot the loaded Data	
	surf(X,Y,psiB)
	surf(X,Y+2,psiT)
	surf(X+2,Y+2,psiTR)

	title(["h=", num2str(3) ", p=",num2str(n),", Advection-Diffusion, field no.", field])
	axis tight
	colorbar
%	colormap("hot")
	shading interp
	xlabel "x"
	ylabel "y"

endfunction