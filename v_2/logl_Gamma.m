function ll = logl_Gamma(y, r, lambda)
  ll = sum(-y./lambda - r*log(lambda) - log(gamma(r)) + (r-1)*log(y));