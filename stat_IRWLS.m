%- TODO:
%- Add automatic file saving functionality
%- Add statistical inference functionalities

function beta_new = stat_IRWLS(X, y, family, n = 0, alpha = 0, epsilon = 10^-6, use_OLS = true)
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
  endif
  
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
  
  