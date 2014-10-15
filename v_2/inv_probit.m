%- Return pi
function ginv = inv_probit(eta)
  ginv = normcdf(eta);