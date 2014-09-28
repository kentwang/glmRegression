%- TODO: - add canonical option
%-       - add statistical inference

function beta = IRWLS(X, y, n, family, epsilon)
	beta = 1;
	if strcmp(family, "binomial")
		disp(1);
	endif