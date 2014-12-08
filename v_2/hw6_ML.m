% method from page 175, (Generalized, Linear, and Mixed Models) by McCulloch, Searle, and Neuhaus
% model: y = X*beta + Z*U. X is the design matrix, Z is the part of random effects in the design matrix
% Question: Whole plots have both fixed and random effects?

data = csvread('ST615_WindTunnel_Design.csv');
X = [ones(size(data, 1), 1), data(:, 1:10)];

% define the four response
y1 = data(:, 11);
y2 = data(:, 12);
y3 = data(:, 13);
y4 = data(:, 14);

y = y1; % pick one of the responses

% define some dimension parameters
N = 32; % total subplots
r = 8; % number of runs  of wp
n = 4; % number of sp in wp
p = 11; % dimension of fixed effect
T = 2; % there are only two variance components

%% Use the notations in Notes

% define the variance matrices for wp No need for updating
Ve = eye(N, N); %subplot error
M = [];
for i=1:r
  M = blkdiag(M, repmat(1, n, 1));
endfor
V = {Ve, M*M'};

% initialization beta, phi and thus Vy
beta_new = [1; zeros(p-1, 1)]; % fix effect coefficient
Phi_new = [1; 0.1];
Vy = Phi_new(1)*Ve + Phi_new(2)*M*M';

beta_old = [1; ones(p-1, 1)];
Phi_old = [10; 0];

epsilon = 10^-6;
iter = 0;


%%%%%%%%Start the iteration here
while(max(abs(beta_new-beta_old))/sum(abs(beta_old)) > epsilon || max(abs(Phi_new-Phi_old))/sum(abs(Phi_old)) > epsilon)
  iter = iter + 1;
  disp(['iteration', num2str(iter)]);
%  disp(beta_new');
  
  % replace old parameters
  beta_old = beta_new;
  Phi_old = Phi_new;
  
  % updating Omega and Rho
  Omega = [];
  Rho = [];

  for t = 1:T
    for m = 1:T
      Omega(t, m) = trace(inv(Vy)*V{t}*inv(Vy)*V{m});
    endfor
  endfor

  for t = 1:T
    Rho(t, 1) = (y-X*beta_old)'*inv(Vy)*V{t}*inv(Vy)*(y-X*beta_old); % Here is a problem with the iteration 
  endfor

  % updating variance components
  Phi_new = inv(Omega)*Rho;

  % updating V
  Vy = Phi_new(1)*Ve + Phi_new(2)*M*M';;

  % updating beta
  beta_new = inv(X'*inv(Vy)*X)*X'*inv(Vy)*y;
endwhile

% variance of estimates
V_beta = inv(X'*inv(Vy)*X);

%% Model summary
% WP effects
disp("--------Intercept---------");
disp(beta_new(1));
disp("--------WP A B AB---------");
disp(beta_new([2, 3, 4])');
disp("--------SP C D CD---------");
disp(beta_new([5, 6, 11])');
disp("--------WP x SP AC AD BC BD---------");
disp(beta_new([7, 8, 9, 10])');

disp('Estimate | se | t-stat | pvalue');
disp([beta_new, sqrt(diag(V_beta)), abs(beta_new) ./ sqrt(diag(V_beta)), 2*(1- tcdf(abs(beta_new) ./ sqrt(diag(V_beta)), 1))]);


disp('Variance Components');
disp(Phi_new);






