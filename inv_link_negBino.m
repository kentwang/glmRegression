function ginv = inv_link_negBino(eta, alpha)
	ginv = alpha .* exp(eta) ./ (1 .- exp(eta));