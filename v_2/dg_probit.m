%- p as the input
function dg = dg_probit(p)
  dg = normpdf(p)./normcdf(p).^2;