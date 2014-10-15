%- Choose is not included Binomial coefficients [log(NchooseK(n, y))]
function ll = logl_Binomial(mu, n, y)
  ll = sum(y.*log(mu./y) + (n-y).*log((n-mu)./n));