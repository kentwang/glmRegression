%----------------------------------------------------------------
% This is to reshape the regression under GLM framework.
% Consider combination of distribution family and link functions.
% Previous version is v_1
%
% Parameters returned
%     - beta_new: Estimated coefficients of the model
%     - CI: 95% confidence intervals of all parameters
%     - seBeta: standard deviation vector matrix of b
%     - Wald: Wald statistics and pvalues
%
% Family: Bernoulli, Binomial, Poisson, NegBino, Gamma
% Link Function: identity, logit(pi), probit(pi), log, reciprocal 
%
%
% Todo: - non-canonical link IRWLS
%       - no-offset poisson regression
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


function [beta_new, CI, seBeta, Wald] = glmReg(X, y, family, link, canonical, n = 0, plotit = true,  use_OLS = false)
  %- Default settings
  epsilon = 10^-8;
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
        beta_old = beta_new;
        beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V))); % pay attention to this
        Vbeta = inverse(X' * V * X);
      else 
        switch(link) % Switch/case the link
        case "probit"
          mu = n.*inv_probit(eta);
          V = diag(var_Binomial(n, mu).*(dg_probit(mu./n)./n).^2); % V_eta, dg^2
          beta_old = beta_new;
          beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* dg_probit(mu./n)./n); %p(eta)/p(mu)
          %- Noncanonical link variance (probit link)
          Delta = diag(n.^2./(mu.*(n-mu)).*normpdf(X * beta_new));
          Vbeta = inverse(X' * Delta * diag(var_Binomial(n, mu)) * Delta * X);
        endswitch
      endif
      deviance = 2*(logl_Binomial(y, n, y) - logl_Binomial(mu, n, y));
      
    %%- Poisson Distribution
    case "Poisson"
      %- Check if there is an offset/exposure n. Maybe not
%      if n == 0
%        disp("There is no offset/exposure for Poisson. Do something else!");
%        return;
%      else
%        y = y./n; % offset poisson
%      endif
      
      if canonical
        mu = n.*inv_log(eta);
        V = diag(var_Poisson(mu));
        beta_old = beta_new;
        beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V)));
        Vbeta = inverse(X' * V * X);
      else
        switch(link)
        case "logit"
          mu = n.*inv_logit(eta); % rate of event
          V = diag(var_Poisson(mu).*(dg_logit(mu./n)./n).^2); % This is a little tricky, V_eta, dg^2
          beta_old = beta_new;
          beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* dg_logit(mu./n)./n);
          %- Noncanonical link vraince (logit link)
          Delta = diag((n-mu)./n);
          Vbeta = inverse(X' * Delta * diag(var_Poisson(n, mu)) * Delta * X);
        endswitch
      endif
      deviance = 2*(logl_Poisson(y, y) - logl_Poisson(mu, y))
    
    case "Gamma" %- Only consider log link for gamma
      r = 2; % scale parameter for process variance model
      mu = inv_log(eta);
      V = diag(var_Gamma(r, mu/r).*dg_log(mu).^2); % ignore a
      beta_old = beta_new;
      beta_new = inverse(X' * inverse(V) * X) * X' * inverse(V) * (eta .+ (y .- mu) .* dg_log(mu));
      Vbeta = inverse(X'*X)/r;
      deviance = 2*(logl_Gamma(y, r, y/r) - logl_Gamma(y, r, mu/r));
    endswitch    
    
    fdisp(f_id, [deviance; beta_new]');
%    fprintf(f_id, strcat(repeat('%f, ', length(beta_new)), ' %f'), [deviance; beta_new])
%    fprintf(f_id, '%f, %f, %f\n', [deviance; beta_new]);  
  endwhile 
  fclose(f_id);  
  
  %- Convergence plot of coefficients and deviance
  if plotit
    M = dlmread(histfile);
    plot_convergence(M);
  endif
  
  %- Confidence Interval
  seBeta = sqrt(diag(Vbeta));
  CILower = beta_new - norminv(1 - 0.5 * alpha_inf) * seBeta;
  CIUpper = beta_new + norminv(1 - 0.5 * alpha_inf) * seBeta;
  CI = [CILower, CIUpper];
  
  %- Wald test on invidual parameters
  WaldStat = beta_new ./ seBeta;
  WaldP = 1 - chi2cdf(WaldStat.^2, 1);
  Wald = [WaldStat, WaldP];
  
