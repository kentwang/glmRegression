load("exam1_dataset.mat");
X = [ones(length(x), 1), x];
%z0 = binornd(1, 0.5, length(y), 1);

%- Test 1, intuitively define z according to the signs of y
z0 = (y == 0);
[beta_new, gamma_new, z] = poisZeroInflatedEM(X, y, z0, true);

%- Test 2, generate Bernoulli random variables for membership
%- No need They are pretty much the same
%p = 0.1:0.05:0.99;
%%nrep = 50;
%BetaGamma = repmat(length(p), 4);
%for i = 1:length(p)
%  printf('i is %d\n', i);
%  z0 = binornd(1, p(i), length(y), 1);
%  [beta_new, gamma_new, z] = poisZeroInflatedEM(X, y, z0);  
%  BetaGamma(i, :) = [beta_new', gamma_new'];  
%endfor
% 0.3879932   0.0081711  -6.4244979   0.0316391

z0 = binornd(1, 0.1, length(y), 1);
[beta_new, gamma_new, z] = poisZeroInflatedEM(X, y, z0, true);