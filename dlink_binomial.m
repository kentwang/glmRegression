%- This is the derivative link function of binomial, input n and mu (count)
function dg = dlink_binomial(n, mu)
	dg = n ./ (mu .* (n - mu));