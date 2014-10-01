%- TODO:
%- Add automatic file saving functionality
%- Add statistical inference functionalities

function beta = stat_IRWLS(X, y, family, n = 0, alpha = 0, epsilon = 10^-6, use.OLS = true)
  if size(X)(1) ~= length(y)
    printf("Dimension of X is %d x %d and y is %d\n", [size(X), length(y)]);
  endif  
  
  %- Initialization and conditions	
  beta_old = 10000 * ones(size(X)(2), 1);
  if use.OLS %use OLS as initial value
    beta_new = OLS(X, y); % Initial value of beta
  endif
  beta_new = [1; 0; 0; 0];
  iter = 0;
  