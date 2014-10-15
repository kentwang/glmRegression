%- input parameter mu and observed count y
%- Ignore factorials
function ll = logl_Poisson(mu, y)
  ll = sum(-mu + y.*log(mu));