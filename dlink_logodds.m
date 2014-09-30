%- This is the derivative link function of logodds
%- TODO: catch error when mu is not in (0, 1)
function dg = dlink_logodds(mu)
	if mu < 0 || mu > 1
		disp("Caution: Input of link_logodds is not in (0, 1)!");
	endif
	dg = 1 / (mu * (1 - mu));