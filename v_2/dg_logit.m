function dg = dg_logit(p)
  dg = 1./(p.*(1-p));