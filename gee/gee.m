% This function is a wrapper of the GEE for both normal and GLM case
% Todo: Convergence issue for logistic regression
function result = gee(Y, Z, s, t, workCor = "Independent", family = "Normal", n = 0, mdep = 1, epsilon = 10^-6, OLS = true)
  disp(OLS);
  if(OLS)
    b0 = inverse(Z'*Z)*Z'*Y;
  elseif(strcmp(family, "Poisson"))
    b0 = inverse(Z'*Z)*Z'*log(Y + 0.00001);
  elseif(strcmp(family, "Binomial"))
%    b0 = zeros(size(Z, 2), 1);
%    b0(1) = 1;
    p = sum(reshape(Y, t, s))' ./ n;
    for i=1:length(p)
      if(p(i) == 0)
        p(i) = 0.00001;
      elseif(p(i) == 1)
        p(i) = 0.99999;
      endif
    endfor
    eta = vec(repmat(log(p./(1-p)), 1, t)');
    b0 = inverse(Z'*Z)*Z'*eta;
  else
    b0 = ones(size(Z, 2), 1);
  endif
  bnew = b0;
  bold = ones(length(b0), 1);
  iter = 0;

  % note residual sum of squares is not a good criterion now
  while(max(abs(bnew - bold)) / sum(abs(bold)) > epsilon)
    iter++;
    printf("Iteration %d\n", iter);
    bold = bnew;
    
    if(strcmp(family, "Normal"))
      [bnew,sigma2,MVb,EVb,R] = geeNormal(Y,Z,s,t,bnew,workCor);
    else
      [bnew,sigma2,MVb,EVb,R] = geeGLM(Y,Z,s,t,bold,workCor,family,n,mdep); #only binomial/logistic GEE considered now, change b0 to bold
    endif
  endwhile
  result = struct("bnew", bnew, "sigma2", sigma2, "MVb", MVb, "EVb", EVb, "R", R);
endfunction