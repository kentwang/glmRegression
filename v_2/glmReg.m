%----------------------------------------------------------------
% This is to reshape the regression under GLM framework.
% Consider combination of distribution family and link functions.
%
% Family: Bernoulli, Binomial, Poisson, NegBino, Gamma
% Link Function: identity, logit(pi), probit(pi), log, reciprocal 
%
% Previous version is v_1
%
% Todo: - non-canonical link IRWLS
%       - no-offset poisson regression
%       - More output information
%       - fprintf() needs to be modified
%----------------------------------------------------------------
%----------------------------------------------------------------
%- Auxiliary input: 
%       - n vector for Binomial family
%       - alpha for Negative Binomial family
%
%- Program logit link and Binomial family for demonstration pf HW 4
%----------------------------------------------------------------
%----------------------------------------------------------------
%- Test run: hw4.m


function beta_new = glmReg(X, y, family, link, canonical, n = 0)
  %- Default settings
  epsilon = 10^-8;
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
  histfile = strcat(family, "_", link, ".hist"); % not need to attach beta_sbeta
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
  %- Quantities to switch/case: Mu, V, Likelihood, Deviance (some of them are sudo)
  f_id = fopen(histfile, "a");
  while(sum(abs(beta_new - beta_old)) / sum(abs(beta_old)) > epsilon)
    iter += 1;
    eta = X * beta_new;
    printf('Iteration %d\n', iter);
    
    switch(family) % Switch/case the family
    %%- Binomial Distribution
    case "Binomial"
      if canonical
        mu = n.*inv_logit(eta);
        V = diag(var_Binomial(n, mu));
      else 
        switch(link) % Switch/case the link
        case "probit"
          mu = n.*inv_probit(eta);
          V = diag(var_Binomial(n, mu).*n.*dg_probit(mu./n)); % ignore a(\phi). Times n?
        endswitch
      endif
      deviance = 2*(logl_Binomial(y, n, y) - logl_Binomial(mu, n, y));
      
    %%- Poisson Distribution
    case "Poisson"
      %- Check if there is an offset/exposure n. Maybe not
      if n == 0
        disp("There is no offset/exposure for Poisson. Do something else!");
        return;
      else
        y = y./n; % offset poisson
      endif
      
      if canonical
        mu = inv_log(eta);
        V = diag(var_Poisson(mu));
      else
        switch(link)
        case "logit"
          mu = inv_logit(eta); % rate of event
          V = diag(var_Poisson(mu).*dg_logit(mu));
        endswitch
      endif
      deviance = 2*(logl_Poisson(y, y) - logl_Poisson(n.*mu, y))
    endswitch
    beta_old = beta_new;
    beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V))); % pay attention to this
    
    fprintf(f_id, '%f, %f, %f\n', [deviance; beta_new]);  
  endwhile 
  fclose(f_id);  
  
  %- Convergence plot of coefficients and deviance
  M = dlmread(histfile);
  plot_convergence(M);
  
