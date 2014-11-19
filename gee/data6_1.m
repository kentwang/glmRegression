Yraw = [30, 34, 29, 28, 31, 31, 31, 35, 32;
     35, 41, 26, 32, 36, 30, 37, 40, 34;
     37, 38, 33, 40, 42, 32, 41, 39, 39;
     36, 42, 36, 41, 40, 40, 40, 44, 45];
     
Y = [vec(Yraw(:, [1, 4, 7])');
          vec(Yraw(:, [2, 5, 8])');
          vec(Yraw(:, [3, 6, 9])')];
          
Z = [1, 1, 200; % Z without dummies, data problem: it is 250, not 150
     1, 1, 200;
     1, 1, 200;
     1, 1, 225;
     1, 1, 225;
     1, 1, 225;
     1, 1, 250;
     1, 1, 250;
     1, 1, 250;
     1, 1, 275;
     1, 1, 275;
     1, 1, 275;
     1, 0, 200;
     1, 0, 200;
     1, 0, 200;
     1, 0, 225;
     1, 0, 225;
     1, 0, 225;
     1, 0, 250;
     1, 0, 250;
     1, 0, 250;
     1, 0, 275;
     1, 0, 275;
     1, 0, 275;
     1, -1, 200;
     1, -1, 200;
     1, -1, 200;
     1, -1, 225;
     1, -1, 225;
     1, -1, 225;
     1, -1, 250;
     1, -1, 250;
     1, -1, 250;
     1, -1, 275;
     1, -1, 275;
     1, -1, 275];

%Z = [1, 1, 0, 200; % Z with dummies
%     1, 1, 0, 200;
%     1, 1, 0, 200;
%     1, 1, 0, 225;
%     1, 1, 0, 225;
%     1, 1, 0, 225;
%     1, 1, 0, 150;
%     1, 1, 0, 150;
%     1, 1, 0, 150;
%     1, 1, 0, 275;
%     1, 1, 0, 275;
%     1, 1, 0, 275;
%     1, 0, 1, 200;
%     1, 0, 1, 200;
%     1, 0, 1, 200;
%     1, 0, 1, 225;
%     1, 0, 1, 225;
%     1, 0, 1, 225;
%     1, 0, 1, 150;
%     1, 0, 1, 150;
%     1, 0, 1, 150;
%     1, 0, 1, 275;
%     1, 0, 1, 275;
%     1, 0, 1, 275;
%     1, 0, 0, 200;
%     1, 0, 0, 200;
%     1, 0, 0, 200;
%     1, 0, 0, 225;
%     1, 0, 0, 225;
%     1, 0, 0, 225;
%     1, 0, 0, 150;
%     1, 0, 0, 150;
%     1, 0, 0, 150;
%     1, 0, 0, 275;
%     1, 0, 0, 275;
%     1, 0, 0, 275];

% problem 6.1 Data processing obmitted
s = 12;
t = 3;

hist(Y);
xlabel("Paper strength")
ylabel("Frequency");

fit_exchange = gee(Y, Z, s, t, workCor = "Exchangeable");
fit_ar1 = gee(Y, Z, s, t, workCor = "AR1");

% data summary  
% exchangeable
[fit_exchange.bnew, ...
sqrt(diag(fit_exchange.MVb)), ...
fit_exchange.bnew ./ sqrt(diag(fit_exchange.MVb)), ...
2*(1-normcdf(fit_exchange.bnew ./ sqrt(diag(fit_exchange.MVb)))), ...
sqrt(diag(fit_exchange.EVb)), ...
fit_exchange.bnew ./ sqrt(diag(fit_exchange.EVb)), ...
2*(1-normcdf(fit_exchange.bnew ./ sqrt(diag(fit_exchange.EVb))))]

fit_exchange.MVb
fit_exchange.EVb

% AR(1)
[fit_ar1.bnew, ...
sqrt(diag(fit_ar1.MVb)), ...
fit_ar1.bnew ./ sqrt(diag(fit_ar1.MVb)), ...
2*(1-normcdf(fit_ar1.bnew ./ sqrt(diag(fit_ar1.MVb)))), ...
sqrt(diag(fit_ar1.EVb)), ...
fit_ar1.bnew ./ sqrt(diag(fit_ar1.EVb)), ...
2*(1-normcdf(fit_ar1.bnew ./ sqrt(diag(fit_ar1.EVb))))]

fit_ar1.MVb
fit_ar1.EVb

