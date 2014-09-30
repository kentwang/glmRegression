%- This is the inverse link function of binomial count. For calculation of mu, expected count out of n
function ginv = inv_link_binomial(n, eta)
	ginv = n * exp(eta) ./ (1 + exp(eta)); % notice the dot operation