%----------------------------------------------------------------
% This is to reshape the regression under GLM framework.
% Consider combination of distribution family and link functions.
%
% Family: Bernoulli, Binomial, Poisson, NegBino, Gamma
% Link Function: identity, logit(pi), probit(pi), log, reciprocal 
%
% Previous version is v_1
%
% Todo: non-canonical link IRWLS
%----------------------------------------------------------------

%- Auxiliary input: 
%       - n vector for Binomial family
%       - alpha for Negative Binomial family
%
%- Program logit link and Binomial family for demonstration pf HW 4



function beta_new = glmReg(X, y, family, link, canonical, n)
  %- Default settings
  epsilon = 10^-6;
  use_OLS = false;
  alpha_inf = 0.05;
  iter = 0;
  
  %- Argument validation
  argFlag = glmRegArgValid(X, y, family, n);
  if argFlag > 0
    disp("Argument error");
    return;
  endif
  
  %- Check and delete .hist files
  histfile = strcat(family, ".hist"); % not need to attach beta_sbeta
  if exist(histfile, "file") % histfile is already a string
    delete(histfile);
  endif
  
  %- Initialization and conditions	
  beta_old = 10000 * ones(size(X)(2), 1);
  if use_OLS %use OLS as initial value
    beta_new = OLS(X, y); % Initial value of beta
  else
    beta_new = [1, zeros(1, size(X)(2)-1)]';  
  endif
  
  %- IRWLS algorithm
  %- Quantities to switch/case: Mu, V, Likelihood, Deviance
  f_id = fopen(histfile, "a");
  while(sum(abs(beta_new - beta_old)) / sum(abs(beta_old)) > epsilon)
    iter += 1;
    eta = X * beta_new;
    printf('Iteration %d\n', iter);
    
    if canonical % if not canonical, swtich/case on link
      switch(family)
      case "Binomial"
        mu = n.*inv_logit(eta);
        V = diag(var_Binomial(n, mu));
        deviance = 2*(logl_Binomial(y, n, y) - logl_Binomial(mu, n, y));
      endswitch
      beta_old = beta_new;
      beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V))); % pay attention to this
    endif
    
    fprintf(f_id, '%f, %f, %f\n', [deviance; beta_new]);  
  endwhile 
  fclose(f_id);  
  
