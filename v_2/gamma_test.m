epsilon = 10^-8;
use_OLS = false;
alpha_inf = 0.05;
iter = 0;
r = 2;

beta_old = 10000 * ones(size(X)(2), 1);
beta_new = [1, zeros(1, size(X)(2)-1)]';  

group = [1; 1; 1; 2; 2; 2; 3; 3; 3; 4; 4; 4; 5; 5; 5; 6; 6; 6; 7; 7; 7; 8; 8; 8];
beta_OLS = inverse(X'*X)*X'* y;
s2 = (y - X*beta_OLS).^2;
scatter(group, s2);


iter += 1;
eta = X * beta_new;
printf('Iteration %d\n', iter);% scale parameter for process variance model
mu = inv_log(eta);
V = diag((1/r)*var_Gamma(r, mu/r).*dg_log(mu).^2);
beta_old = beta_new;
beta_new = inverse(X' * V * X) * X' * V * (eta .+ (s2 .- mu) .* dg_log(mu));
Vbeta = inverse(X'*X)/r;
deviance = 2*(logl_Gamma(s2, r, s2/r) - logl_Gamma(s2, r, mu/r));

disp([deviance;beta_new]');