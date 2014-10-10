%- TODO: - add canonical option
%-       - add statistical inference
%-       - add likelihood function or S(beta)
%-       - AIC and BIC maybe added

% function beta = IRWLS(X, y, n, family, epsilon)
	% beta = 1;
	% if strcmp(family, "binomial")
		% disp(1);
	% endif

data = csvread("HW3.csv");
X = [ones(size(data)(1), 1), data(:, 1:3)];
y = data(:, 4);

%- Initialization and conditions	
epsilon = 10^-6;
beta_old = 10000 * ones(size(X)(2), 1);
beta_new = OLS(X, y); % Initial value of beta
beta_new = [1; 0; 0; 0];
alpha = 0.5;
alpha_inf = 0.05
iter = 0;

%- Check existance of log files 
if exist("beta_sbeta_negBino.hist")
	delete "beta_sbeta_negBino.hist";
endif

fbeta_id = fopen("beta_sbeta_negBino.hist", "a");
while(sum(abs(beta_new - beta_old)) / sum(abs(beta_old)) > epsilon)
	iter += 1;
  printf("iteration %d\n", iter);
	eta = X * beta_new;
  %	mu = inv_link_negBino(eta, alpha);
  mu = exp(eta); % use the non-canonical link function for computation instead of the canonical
	V = diag(var_negBino(mu, alpha));
	beta_old = beta_new;
	beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V))); % pay attention to this
  deviance = -2 * loglik_negBino(mu, y, alpha); 
  fprintf(fbeta_id, '%f, %f, %f, %f, %f\n', [deviance; beta_new]);
endwhile
fclose(fbeta_id);

M = dlmread("beta_sbeta_negBino.hist");

subplot(3, 2, 1); plot(M(:, 1)); title('\beta_0'); xlabel("Iteration");
subplot(3, 2, 2); plot(M(:, 2)); title('\beta_1'); xlabel("Iteration");
subplot(3, 2, 3); plot(M(:, 3)); title('\beta_2'); xlabel("Iteration");
subplot(3, 2, 4); plot(M(:, 4)); title('\beta_3'); xlabel("Iteration");
subplot(3, 2, 5:6); plot(M(:, 5)); title('Deviance: \lambda(\beta)'); xlabel("Iteration");
suptitle("NB regression convergence using OLS as starting value for IRWLS");

%- Confidence Interval
varBeta = inverse(X' * V * X);
seBeta = sqrt(diag(varBeta));
CILower = beta_new - norminv(1 - 0.5 * alpha_inf) * seBeta;
CIUpper = beta_new + norminv(1 - 0.5 * alpha_inf) * seBeta;
CI = [CILower, CIUpper];

%- Wald test on invidual parameters
WaldStat = beta_new ./ seBeta;
WaldF = chi2inv(1 - alpha_inf, 1);
WaldP = 1 - chi2cdf(WaldStat.^2, 1);

%- Likelihood ratio test






