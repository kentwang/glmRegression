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
  beta_new = 
  gamma_new = 
  while()