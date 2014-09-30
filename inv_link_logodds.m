%- This is the inverse link function of logodds. For calculation of mu
%- TODO: catch error when mu is not in (0, 1)
function ginv = inv_link_logodds(eta)
	ginv = exp(eta) ./ (1 + exp(eta)); % notice the dot opertation