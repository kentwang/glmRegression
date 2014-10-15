data = csvread("hw4_3.csv");
X = [ones(size(data)(1), 1), data(:, 1)];
n = data(:, 2);
y = data(:, 3);


% Test run format glmReg(X, y, family, link, canonical, n = 0)

%- Logistic regression + logit
beta_logistic_logit = glmReg(X, y, "Binomial", "logit", true, n);

%- Logistic regression + probit
beta_logistic_probit = glmReg(X, y, "Binomial", "probit", false, n);

%- Poisson and logit test
glmReg(X, y, "Poisson", "logit", false, n);