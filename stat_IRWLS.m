%- TODO:
%- Add generic function names for link, dlink, var, inv_link, loglik
%- Add statistical inference functionalities
%- Change validating arguments RETURN into WARNING
%- Reduce model estimation and likihood

%- Test variables: family = "negBino"; n = 5; alpha = 0.5;
%data = csvread("HW3.csv");
%X = [ones(size(data)(1), 1), data(:, 1:3)];
%y = data(:, 4);
%- Test run: [beta, CI, varBeta, Wald] = stat_IRWLS(X, y, "negBino", n = 0, alpha = 0.5);

function [beta_new, CI, varBeta, Wald] = stat_IRWLS(X, y, family, n = 0, alpha = -1, epsilon = 10^-6, use_OLS = false, alpha_inf = 0.05)
  %- Validating arguments
  if size(X)(1) ~= length(y)
    printf("Dimension of X is %d x %d and y is %d\n", [size(X), length(y)]);
    return;
  endif

  family_set = {"Bernoulli", "Binomial", "negBino", "Poisson"};
  if ~exist("family", "var") || ~ismember(family, family_set)
    disp("Model family missing! Family can be one of the following:");
    disp(["Bernoulli, ", "Binomial, ", "negBino, ", "Poisson"]);
    return;
  elseif strcmp(family, "Binomial") && n == 0
    disp("n is not given for Binomial regression!");
    return;
  elseif strcmp(family, "negBino") && alpha == -1
    printf('alpha is %f', alpha);
    disp("alpha is not positive for Negative Binomial regression!");
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
  iter = 0;
  
  %- Generic functions involved in IRWLS for GLM
  switch(family)
  case "negBino"
    inv_link_fun = @inv_link_negBino;
    var_fun = @var_negBino;
    loglik_fun = @loglik_negBino;
    reduced_dev = @reduced_negBino_dev;
  endswitch  
  
  %- IRWLS iteration
  f_id = fopen(histfile, "a");
  while(sum(abs(beta_new - beta_old)) / sum(abs(beta_old)) > epsilon)
    iter += 1;
%    printf("Iteration %d\n", iter);
    eta = X * beta_new;
    mu = inv_link_fun(eta);
%    mu = exp(eta); % use the non-canonical link function for computation instead of the canonical
    V = diag(var_fun(mu, alpha));
    beta_old = beta_new;
    beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V))); % pay attention to this
    deviance = -2 * loglik_fun(mu, y, alpha); 
    fprintf(f_id, '%f, %f, %f, %f, %f\n', [deviance; beta_new]);
  endwhile
  fclose(f_id);
  
  %- Convergence plot of coefficients
  M = dlmread(histfile);
  plot_convergence(M);
  
  %- Confidence Interval
  varBeta = inverse(X' * V * X);
  seBeta = sqrt(diag(varBeta));
  CILower = beta_new - norminv(1 - 0.5 * alpha_inf) * seBeta;
  CIUpper = beta_new + norminv(1 - 0.5 * alpha_inf) * seBeta;
  CI = [CILower, CIUpper];
  
  %- Wald test on invidual parameters
  WaldStat = beta_new ./ seBeta;
  WaldP = 1 - chi2cdf(WaldStat.^2, 1);
  Wald = [WaldStat, WaldP];
  
  %- Likelihood ratio test. Calculate reduce deviance on each model
  %- Note simple parameter is test for each iteration and df of LRT is 1
  %%%%%%%%%%%%%% This is funky for now
  deviance_reduce = zeros(1, size(X)(2) - 1);
  for k = 1:(size(X)(2) - 1) % remove the kth feature (not column of X)
    cols = 2:size(X)(2);
    cols(k) = [];
    X_reduce = X(:, cols);
    deviance_reduce(k) = reduced_dev(X_reduce, y, alpha = alpha);
  endfor
  
  

  