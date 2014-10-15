%- Default settings
epsilon = 10^-8;
use_OLS = false;
alpha_inf = 0.05;
iter = 0;

%- Initialization and conditions	
beta_old = 10000 * ones(size(X)(2), 1);
if use_OLS %use OLS as initial value
  beta_new = OLS(X, y); % Initial value of beta
else
  beta_new = [1, zeros(1, size(X)(2)-1)]';  
endif

iter += 1;
eta = X * beta_new;
printf('Iteration %d\n', iter);
mu = n.*inv_logit(eta); % rate of event
V = diag(var_Poisson(mu).*dg_logit(mu)./n); % This is a little tricky
deviance = 2*(logl_Poisson(y, y) - logl_Poisson(mu, y))
beta_old = beta_new;
beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V))); % pay attention to this
fprintf(f_id, '%f, %f, %f\n', [deviance; beta_new]);  