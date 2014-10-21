%-------------------------------------------------------------------
% This script runs EM algorithm on zero-inflated Poisson regression.
% For the glmReg.m, only Logistic regression and Poisson regression 
% are needed.
%
% Parameters returned
%     beta_new: model coefficients for log(\lambda) = B\beta
%     gamma_new: model coefficeints for logit(p) = G\gamma
%     z: the estimated membership vector for all obvservations
%     Note: Here B and G are equivalent to X, the design matrix
%
% Todo: Statistical inference and test on parameters
% Testrun: don't save to .hist files. Too many
%-------------------------------------------------------------------

function [beta_new, gamma_new, z] = poisZeroInflatedEM(X, y, z0, plotit = false, epsilon = 10^-6)
  %- Parameters initialization
  beta_old = gamma_old = repmat(100, size(X)(2), 1);
  z = z0;
  beta_new = glmReg(X, y, "Poisson", "log", true, repmat(1, length(z0), 1), 1 - z, false); % Octave argument is tricky. It will update the holding variables!
  gamma_new = glmReg(X, z, "Binomial", "logit", true, repmat(1, length(z0), 1), 0, false);
  
  %- Check and delete "poisZeroInflated.hist" files
  histfile = "poisZeroInflated.hist"; % not need to attach beta_sbeta
  if exist(histfile, "file") % histfile is already a string
    delete(histfile);
  endif
  
  f_id = fopen(histfile, "a");
  iter = 0;
  while(sum(abs(beta_new - beta_old)) / sum(abs(beta_old)) > epsilon || sum(abs(gamma_new - gamma_old)) / sum(abs(gamma_old)) > epsilon)
    iter += 1;
%    printf('Iteration %d\n', iter);
%    printf('Beta_new, ');
%    disp(beta_new');
%    printf('Gamma_new, ');
%    disp(gamma_new');
    
    % Updating z should be careful, LOL that simple?
    z = 1./(1 + exp(-X*gamma_new - exp(X*beta_new))).*(y == 0);
    beta_old = beta_new;
    gamma_old = gamma_new;
    beta_new = glmReg(X, y, "Poisson", "log", true, repmat(1, length(z0), 1), 1 - z, false);
    gamma_new = glmReg(X, z, "Binomial", "logit", true, repmat(1, length(z0), 1), 0, false);
    
    fdisp(f_id, [beta_new', gamma_new']);
  endwhile
  fclose(f_id); 
  
  if plotit
    M = dlmread(histfile);
    subplot(2, 2, 1); plot(M(:, 1)); title('\beta_0');
    subplot(2, 2, 2); plot(M(:, 2)); title('\beta_1');
    subplot(2, 2, 3); plot(M(:, 3)); title('\gamma_0');
    subplot(2, 2, 4); plot(M(:, 4)); title('\gamma_1');
  endif
  