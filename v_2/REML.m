beta0 = [1; repmat(0, p-1, 1)];
theta0 = 0.5;
eta0 = 0.1;

beta_new = beta0;
theta_new = theta0;
eta_new = eta0;

% define old parameters
beta_old = beta_new;
theta_old = theta_new;
eta_old = eta_new;

ytheta = boxcox(y, theta_old);
err = ytheta - X*beta_old;



% define the M and D matrices
M = [];
for i=1:N
  M = blkdiag(M, repmat(1, n, 1));
endfor
%eta_low = -1/max(eig(M*M'));
MM = M*M';
eta_low = -0.004 % Is this a fixed number?


D = eye(S) + eta_old.*M*M';
invD = inv(D);

% define and solve the function partial(l_R)/paritial(theta) for theta
op = optimset("Display","final");
plr_ptheta = @(theta)(p-S)*(1./(err'*inv(D)*err))*(inv(D)*err)'*(theta*y.^theta.*log(y) - y.^theta + 1)/theta^2 + sum(log(y));
theta_new = fsolve(plr_ptheta, theta_old, optimset("Display", "iter"));



% I can use solver here or Newton-Raphson

% optimization test, maximize the likelihood directly
% fminbnd is finding the minimum. Put a negative sign on everything
t0 = clock();
lrtheta = @(theta) .5*(S-p)*log((boxcox(y, theta) - X*beta_old)'*invD*(boxcox(y, theta) - X*beta_old)) - (theta - 1)*sum(log(y));
fminbnd(lrtheta, -3, 3)
elapsed_time_1 = etime(clock(), t0)

t0 = clock();
lreta = @(eta).5*(S-p)*log(err'*(eye(S) + eta*MM)*err) + ...
          .5*logdet(eye(S) + eta*MM) + ...
          .5*logdet(X'*inv(eye(S) + eta*MM)*X);

fminbnd(lreta, eta_low, 10)
fminsearch(lreta, 1)
elapsed_time_2 = etime(clock(), t0)       


%fminsearch (@(x) (x(1)-5).^2+(x(2)-8).^4, [0;0])


%% Two dimensional optimization seems harder and slower
% use fminsearch with both theta and eta
% param = [theta; eta];%theta = param(1); eta = param(2);
%lr_theta_eta = @(param).5*(S-p)*log((boxcox(y, param(1)) - X*beta_old)'*inv(eye(S) + param(2)*M*M')*(boxcox(y, param(1)) - X*beta_old))...
%                +.5*logdet(eye(S) + param(2)*M*M', 'chol')...
%                +.5*logdet(X'*inv(eye(S) + param(2)*M*M')*X)...
%                +(param(1)-1)*sum(log(y));
%end;
%t0 = clock();
%[param, fval] = fminsearch(lr_theta_eta, [0.5; 0]);                
%elapsed_time = etime(clock(), t0);


%% Two dimensional optimization seems harder
%phi = @(param) .5*(S-p)*log((boxcox(y, param(1)) - X*beta_old)'*inv(eye(S) + param(2)*M*M')*(boxcox(y, param(1)) - X*beta_old))...
%                +.5*logdet(eye(S) + param(2)*M*M')...
%                +.5*logdet(X'*inv(eye(S) + param(2)*M*M')*X)...
%                -(param(1)-1)*sum(log(y));
%
%h = @(param) param(2) - eta_low;
%
%t0 = clock();
%[x, obj, info, iter, nf, lambda] = sqp([0.5; 0], phi, [], h)  
%elapsed_time = etime(clock(), t0);

