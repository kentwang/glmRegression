%- This is the link function of logodds
%- TODO: catch error when mu is not in (0, 1)
function g = link_logodds(mu)
	if mu < 0 || mu > 1
		disp("Caution: Input of link_logodds is not in (0, 1)!");
	endif
	g = log(mu / (1 - mu));