function g = link_negBino(mu, alpha)
	g = log(mu ./ (alpha .+ mu));