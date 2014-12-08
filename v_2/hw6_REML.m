% method from page 175, (Generalized, Linear, and Mixed Models) by McCulloch, Searle, and Neuhaus
% model: y = X*beta + Z*U. X is the design matrix, Z is the part of random effects in the design matrix
% Question: Whole plots have both fixed and random effects?

data = csvread('ST615_WindTunnel_Design.csv');
X = [ones(size(data, 1), 1), data(:, 1:10)];

% random effect matrix Z without the subplot error
Z = X(:, 2:4);

% define the four response
y1 = data(:, 11);
y2 = data(:, 12);
y3 = data(:, 13);
y4 = data(:, 14);

y = y3; % pick one of the responses

% define some dimension parameters
N = 32; % total subplots
r = 8; % number of runs  of wp
n = 4; % number of sp in wp
p = 11; % dimension of fixed effect
T = 4; % dimension of random effect epsiloin A B AB epsilon

%% Use the notations in Notes

% define the variance matrices for wp No need for updating
Ve = eye(N, N); %subplot error
Va = []; %A
Vb = []; %B
Vab = []; %AB
for i=1:r
%  Va = blkdiag(Va, ones(n));
%  Vb = blkdiag(Vb, ones(n));
%  Vab = blkdiag(Vab, ones(n));
  Va = Z(:, 1)*Z(:, 1)';
  Vb = Z(:, 2)*Z(:, 2)';
  Vab = Z(:, 3)*Z(:, 3)';
endfor
V = {Ve, Va, Vb, Vab};

% initialization beta, phi and thus Vy
beta_new = [1; zeros(p-1, 1)]; % fix effect coefficient
Phi_new = [1; 2; 2; 2];
Vy = Phi(1)*Ve + Phi(2)*Va + Phi(3)*Vb + Phi(4)*Vab;
M = eye(N) - X*inv(X'*inv(Vy)*X)*X'*inv(Vy);


beta_old = [1; zeros(p-1, 1)];
Phi_old = [1; 0; 0; 0];

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
      Omega(t, m) = trace(inv(Vy)*M*V{t}*M*inv(Vy)*V{m});
    endfor
  endfor

  for t = 1:T
    Rho(t, 1) = y'*M*inv(Vy)*V{t}*inv(Vy)*M*y;
  endfor

  % updating variance components
  Phi_new = inv(Omega)*Rho;

  % updating V
  Vy = Phi_new(1)*Ve + Phi_new(2)*Va + Phi_new(3)*Vb + Phi_new(4)*Vab;
  
  % updating M
  M = eye(N) - X*inv(X'*inv(Vy)*X)*X'*inv(Vy);
  
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

%tcdf(beta_new ./ sqrt(diag(V_beta)), 1)







