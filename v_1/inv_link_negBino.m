%- This is the inverse link function for negBino
%- Notice that a non-canonical function is used
%- This is actually inverse log-link

%function ginv = inv_link_negBino(eta, alpha)
%	ginv = alpha .* exp(eta) ./ (1 .- exp(eta));

function ginv = inv_link_negBino(eta)
  ginv = exp(eta);