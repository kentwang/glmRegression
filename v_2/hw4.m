%% Problem 3 test run

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


%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 4
%%%%%%%%%%%%%%%%%%%%%%%%

data = csvread("hw4_10.csv");
X = [ones(size(data)(1), 1), data(:, 1:(end-1))];
y = data(:, end);

groupMean = zeros(8, 1);
groupVar = zeros(8, 1);

for i = 1:8
  indice = ((i-1)*3+1):((i-1)*3+3);
  disp(y(indice)');
  groupMean(i) = mean(y(indice));
  groupVar(i) = var(y(indice));
endfor

scatter(groupMean, groupVar);

%- Try OLS
group = [1; 1; 1; 2; 2; 2; 3; 3; 3; 4; 4; 4; 5; 5; 5; 6; 6; 6; 7; 7; 7; 8; 8; 8];
beta_OLS = inverse(X'*X)*X'* y;
s2 = (y - X*beta_OLS).^2;
scatter(group, s2);

[b4, CI4, sb4, wald4] = glmReg(X, s2, "Gamma", "log", true, 0, true);