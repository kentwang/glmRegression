%- Test run: reduced_negBino_dev(X, y, alpha = 0.5) Succeful Run!
function deviance = reduced_negBino_dev(X, y, alpha = -1, epsilon = 10^-6, use_OLS = false)
  %- Note X is the reduce design matrix
  %- Validating arguments
  if size(X)(1) ~= length(y)
    printf("Dimension of X is %d x %d and y is %d\n", [size(X), length(y)]);
    return;
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
  inv_link_fun = @inv_link_negBino;
  var_fun = @var_negBino;
  loglik_fun = @loglik_negBino;
  
  %- IRWLS iteration
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
  endwhile
  