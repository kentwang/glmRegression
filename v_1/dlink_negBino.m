function dg = dlink_negBino(mu, alpha)
	dg = alpha ./ (mu .* (alpha .+ mu));