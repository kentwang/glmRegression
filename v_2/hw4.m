data = csvread("hw4_3.csv");
X = [ones(size(data)(1), 1), data(:, 1)];
n = data(:, 2);
y = data(:, 3);


% Test run format glmReg(X, y, family, link, canonical, n = 0)

%- Logistic regression + logit
[b1, CI1, sb1, wald1] = glmReg(X, y, "Binomial", "logit", true, n);

%- Logistic regression + probit
[b2, CI2, sb2, wald2] = glmReg(X, y, "Binomial", "probit", false, n);

%- Poisson and logit test
[b3, CI3, sb3, wald3] = glmReg(X, y, "Poisson", "logit", false, n);