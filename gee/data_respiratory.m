data = csvread("respiratory.csv");
Y = data(:, 4);
Z = [ones(length(Y), 1), data(:, [1, 2, 3])];

s = 56;
t = 5; % the baseline measurement has been included
n = t * ones(s, 1);
b_ols = inverse(Z'*Z)*Z'*Y;


% do every case except for m-dependent first
workCorSeq = {"Unspecified", "Exchangeable", "Independent", "AR1"};

logit_unspecified = gee(Y, Z, s, t, workCor = "Unspecified", family = "Binomial", n = n);
logit_exchangeable = gee(Y, Z, s, t, workCor = "Exchangeable", family = "Binomial", n = n);
logit_independent = gee(Y, Z, s, t, workCor = "Independent", family = "Binomial", n = n);
logit_AR1= gee(Y, Z, s, t, workCor = "AR1", family = "Binomial", n = n);
logit_1dep = gee(Y, Z, s, t, workCor = "m-dependent", family = "Binomial", n = n, mdep = 1);
logit_2dep = gee(Y, Z, s, t, workCor = "m-dependent", family = "Binomial", n = n, mdep = 2);
logit_3dep = gee(Y, Z, s, t, workCor = "m-dependent", family = "Binomial", n = n, mdep = 3);