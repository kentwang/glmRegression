function flag = glmRegArgValid(X, y, family, n)
  %- Validating arguments
  flag = 0; % This means valid arguments
  if size(X)(1) ~= length(y)
    printf("Dimension of X is %d x %d and y is %d\n", [size(X), length(y)]);
    flag = 1; % Dimension not matched error
  endif

  family_set = {"Bernoulli", "Binomial", "NegBino", "Poisson"};
  if ~exist("family", "var") || ~ismember(family, family_set)
    disp("Model family missing! Family can be one of the following:");
    disp(["Bernoulli, ", "Binomial, ", "NegBino, ", "Poisson"]);
    flag = 2; % family not found error
  elseif strcmp(family, "Binomial") && n == 0
    disp("n is not given for Binomial regression!");
    flag = 3; % Binomial no n error
  elseif strcmp(family, "NegBino") && alpha == -1
    printf('alpha is %f', alpha);
    disp("alpha is not positive for Negative Binomial regression!");
    flag = 4; % Negative Binomial no alpha error
  endif
  
  
