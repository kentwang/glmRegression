s = 59;
t = 5; % t should be 5 not 4, baseline should be included

tre_y_raw = dlmread("epilepsy_tre_y.txt", sep = "\t");
Y = tre_y_raw(:, 2);
tre = tre_y_raw(:, 1);

% reshape the tre, should be 59*5
tre_temp =  reshape(tre, 4, 59);
tre_temp = [tre_temp; tre_temp(1, :)];
tre = vec(tre_temp);

base_age_raw = dlmread ("epilepsy_base_age.txt", sep = " ");
base = real(base_age_raw);
age = imag(base_age_raw);

Y = vec([reshape(Y, 4, 59)', base]');

Z = [ones(length(Y), 1), tre, vec(repmat(age, 1, t)')]; % baseline is not a factor

logb_OLS = inverse(Z'*Z)*Z'*log(Y + 0.00001); # try log transformation on response. Nonconvergence to start with OLS
fit = gee(Y, Z, s, t, workCor = "AR1", family = "Poisson", OLS = false);

% sumary
btable = [fit.bnew, ...
sqrt(diag(fit.MVb)), ...
fit.bnew ./ sqrt(diag(fit.MVb)), ...
2*(1-normcdf(abs(fit.bnew) ./ sqrt(diag(fit.MVb)))), ...
sqrt(diag(fit.EVb)), ...
fit.bnew ./ sqrt(diag(fit.EVb)), ...
2*(1-normcdf(abs(fit.bnew) ./ sqrt(diag(fit.EVb))))];
