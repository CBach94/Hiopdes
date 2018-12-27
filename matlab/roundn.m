function [y] = roundn (x,n)
	scale = 10^n;
	y = round(x*scale)/scale;
endfunction
