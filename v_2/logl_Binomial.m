%- Choose is not included Binomial coefficients [log(NchooseK(n, y))]
%- Note, y is modified for logistic regression. Plus a noise
function ll = logl_Binomial(mu, n, y)
  if unique(n) == 1
    for i = 1:length(y)
      if y(i) == 0
        y(i) = 0.00001;
      endif
    endfor
  endif
  ll = sum(y.*log(mu./y) + (n-y).*log((n-mu)./n));