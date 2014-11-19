Yraw = [63, 89;
        54, 91;
        61, 62;
        50, 80;
        52, 72;
        59, 69;
        48, 73;
        74, 81;
        71, 69;
        54, 88;
        48, 92;
        59, 64];
Y = [Yraw([1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12], 1); Yraw([1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12], 2)];

%Z = [1, -1, 0; % this is wrong. Time is a continuous variable
%     1, -1, 0;
%     1, -1, 0;
%     1, -1, 0;
%     1, 0, 0;
%     1, 0, 0;
%     1, 0, 0;
%     1, 0, 0;
%     1, 1, 0;
%     1, 1, 0;
%     1, 1, 0;
%     1, 1, 0;
%     1, -1, 1;
%     1, -1, 1;
%     1, -1, 1;
%     1, -1, 1;
%     1, 0, 1;
%     1, 0, 1;
%     1, 0, 1;
%     1, 0, 1;
%     1, 1, 1
%     1, 1, 1;
%     1, 1, 1;
%     1, 1, 1;];
%     
     
Z = [1, 10, 0;
     1, 10, 0;
     1, 10, 0;
     1, 10, 0;
     1, 20, 0;
     1, 20, 0;
     1, 20, 0;
     1, 20, 0;
     1, 30, 0;
     1, 30, 0;
     1, 30, 0;
     1, 30, 0;
     1, 10, 1;
     1, 10, 1;
     1, 10, 1;
     1, 10, 1;
     1, 20, 1;
     1, 20, 1;
     1, 20, 1;
     1, 20, 1;
     1, 30, 1
     1, 30, 1;
     1, 30, 1;
     1, 30, 1;];
s = 6;
t = 4;

fit_exchange = gee(Y, Z, s, t, workCor = "Exchangeable");
fit_ar1 = gee(Y, Z, s, t, workCor = "AR1");

% summary
% exchangeable
btable1 = [fit_exchange.bnew, ...
sqrt(diag(fit_exchange.MVb)), ...
fit_exchange.bnew ./ sqrt(diag(fit_exchange.MVb)), ...
2*(1-normcdf(abs(fit_exchange.bnew) ./ sqrt(diag(fit_exchange.MVb)))), ...
sqrt(diag(fit_exchange.EVb)), ...
fit_exchange.bnew ./ sqrt(diag(fit_exchange.EVb)), ...
2*(1-normcdf(abs(fit_exchange.bnew) ./ sqrt(diag(fit_exchange.EVb))))];

fit_exchange.MVb
fit_exchange.EVb


% AR(1)
btable2 = [fit_ar1.bnew, ...
sqrt(diag(fit_ar1.MVb)), ...
fit_ar1.bnew ./ sqrt(diag(fit_ar1.MVb)), ...
2*(1-normcdf(abs(fit_ar1.bnew) ./ sqrt(diag(fit_ar1.MVb)))), ...
sqrt(diag(fit_ar1.EVb)), ...
fit_ar1.bnew ./ sqrt(diag(fit_ar1.EVb)), ...
2*(1-normcdf(abs(fit_ar1.bnew) ./ sqrt(diag(fit_ar1.EVb))))];

fit_ar1.MVb
fit_ar1.EVb

