%- This is the variance function of binomial distribution using mu and n as the input
%- Note mu is the success count
function v = var_binomial(n, mu)
	v = n. * mu .* (n.- mu);