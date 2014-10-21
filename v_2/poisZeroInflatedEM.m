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
% Testrun: ???
%-------------------------------------------------------------------

function [beta_new, gamm_new, z] = poisInfEM(X, y, z0, epsilon = 10^-6)
  beta_old = gamma_old = repmat(100, size(X)(2), 1);
  z = z0;
  beta_new = glmReg(X, y, "Poisson", "log", true, repmat(1, length(z0), 1), 1 - z, false);
  gamma_new = glmReg(X, z, "Binomial", "logit", true, repmat(1, length(z0), 1), 0, false);
  
  iter = 0;
  while(sum(abs(beta_new - beta_old)) / sum(abs(beta_old)) > epsilon || sum(abs(gamma_new - gamma_old)) / sum(abs(gamma_old)) > epsilon)
    iter += 1;
    printf('Iteration %d\n', iter);
    printf('Beta_new, ');
    disp(beta_new');
    printf('Gamma_new, ');
    disp(gamma_new');
    
    % Updating z should be careful
    z = 1./(1 + exp(-X*gamma_new - exp(X*beta_new)));
    beta_old = beta_new;
    gamma_old = gamma_new;
    beta_new = glmReg(X, y, "Poisson", "log", true, repmat(1, length(z0), 1), 1 - z, false);
    gamma_new = glmReg(X, z, "Binomial", "logit", true, repmat(1, length(z0), 1), 0, false);
  endwhile