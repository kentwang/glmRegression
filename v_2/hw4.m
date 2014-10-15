data = csvread("hw4_3.csv");
X = [ones(size(data)(1), 1), data(:, 1)];
n = data(:, 2);
y = data(:, 3);


% Test run format glmReg(X, y, family, link, canonical, n = 0)

%- Poisson and logit test
glmReg(X, y, "Poisson", "logit", false, n);