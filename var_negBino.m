function v = var_negBino(mu, alpha)
	v = mu .+ mu.^2 ./alpha;