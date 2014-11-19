function result = gee(Y, Z, s, t, workCor = "Independent", family = "Normal", n = 0, mdep = 1, epsilon = 10^-6)
  b0 = inverse(Z'*Z)*Z'*Y;
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
      [bnew,sigma2,MVb,EVb,R] = geeGLM(Y,Z,s,t,b0,workCor,family,n,mdep); #only binomial/logistic GEE considered now
    endif
  endwhile
  result = struct("bnew", bnew, "sigma2", sigma2, "MVb", MVb, "EVb", EVb, "R", R);
endfunction