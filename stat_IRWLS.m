%- TODO:
%- Add automatic file saving functionality
%- Add statistical inference functionalities

%- Test variables: family = "negBino"; n = 5; alpha = 0.5;

function beta_new = stat_IRWLS(X, y, family, n = 0, alpha = -1, epsilon = 10^-6, use_OLS = false)
  %- Validating arguments
  if size(X)(1) ~= length(y)
    printf("Dimension of X is %d x %d and y is %d\n", [size(X), length(y)]);
    return;
  endif  
  if ~exist("family", "var")
    disp("Model family missing! Family can be one of the following:");
    disp(["Bernoulli, ", "Binomial", "negBino, ", "Poisson"]);
    return;
  elseif strcmp(family, "Binomial") && n == 0
    disp("n is not given for Binomial regression!");
    return;
  elseif strcmp(family, "negBino") && alpha == -1
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
    beta_new = [1; 0; 0; 0];  
  endif
  iter = 0;
  
  %- IRWLS iteration
  f_id = fopen(histfile, "a");
  while(sum(abs(beta_new - beta_old)) / sum(abs(beta_old)) > epsilon)
    iter += 1;
%    printf("Iteration %d\n", iter);
    eta = X * beta_new;
    %	mu = inv_link_negBino(eta, alpha);
    mu = exp(eta); % use the non-canonical link function for computation instead of the canonical
    V = diag(var_negBino(mu, alpha));
    beta_old = beta_new;
    beta_new = inverse(X' * V * X) * X' * V * (eta .+ (y .- mu) .* diag(inverse(V))); % pay attention to this
    deviance = -2 * loglik_negBino(mu, y, alpha); 
    fprintf(f_id, '%f, %f, %f, %f, %f\n', [deviance; beta_new]);
  endwhile
  fclose(f_id);
  
  %- Convergence plot of coefficients
  M = dlmread("beta_sbeta_negBino.hist");
  plot_convergence(M);
    