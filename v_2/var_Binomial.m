%- Input is n and mu 
function V = var_Binomial(n, mu)
  V = mu.*(1-mu./n);