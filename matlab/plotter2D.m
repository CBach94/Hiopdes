function plotter2D(n,N)
	close all
	clc
	
	% Call C++ program which generates data
	system(["./run ", num2str(n), " ", num2str(N)])

	% Relative path to location of fields+meshes
	prefix = "../";

	% Each field is associated to a mesh...
	fields = {"psi.txt"	,"Psi.txt"};
	meshes = {"x.txt"		,"X.txt"};

	% Fetch and plot 2D data.
	for i = 1:length(fields);
		fieldData = dlmread([prefix , fields{i}]);
		fieldData = fieldData(:,1);
		meshData = dlmread([prefix , meshes{i}]);
		
		[X,Y] = meshgrid(meshData,meshData);

		n = floor(sqrt(length(fieldData)));
		Phi = reshape(fieldData(:,1),n,n);

		figure
		hold on
		stem3(X,Y,Phi)
%		contourf(X,Y,Phi)
		surf(X,Y,Phi)
		shading interp;
	
		
		xlabel 'x'
		ylabel 'y'
		colorbar
	end
	

endfunction