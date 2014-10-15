%- Return pi
function ginv = inv_logit(eta)
  ginv = 1./(1 + exp(-eta));
  