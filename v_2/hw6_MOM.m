% model page 546 of DOE by D. Montgomery without three-way interactions

data = csvread('ST615_WindTunnel_Design.csv');

X = [ones(size(data, 1), 1), data(:, 1:10)];

% define effects individually
Const = X(:, 1);
Rep = [ones(16, 1); -ones(16, 1)];
A = X(:, 2);
B = X(:, 3);
AB = X(:, 4);
RepAB = Rep.*AB;
C = X(:, 5);
D = X(:, 6);
AC = X(:, 7);
AD = X(:, 8);
BC = X(:, 9);
BD = X(:, 10);
CD = X(:, 11);

% define the four response
y1 = data(:, 11);
y2 = data(:, 12);
y3 = data(:, 13);
y4 = data(:, 14);

y = y4; % can be replaced

[N, p] = size(data);
i = j = k = l = m = 2;


%%%% Method of Moment for variance components

%% Sum of Squares 
%grand average mu
yddd_b = Const'*y/n;

% SS total
SST = sum((y-mean(y)).^2);

%replicate effect
C_rep = Rep'*y;
SS_rep = C_rep^2/N;
MS_rep = SS_rep/1;
disp("Replicate Effect Rep")
disp(C_rep/N/2);

%wp A
C_a = A'*y;
SS_a = C_a^2/N;
MS_a = SS_a/1;

%wp B
C_b = B'*y;
SS_b = C_b^2/N;
MS_b = SS_b/1;

%wp iteraction AB
C_ab = AB'*y;
SS_ab = C_ab^2/N;
MS_ab = SS_ab/1;

%wp error RepAB
C_repab = RepAB'*y;
SS_repab = C_repab^2/N;
MS_repab = SS_repab/3; % note degrees of freedom
MSE_wp = MS_repab;

%sp C
C_c = C'*y;
SS_c = C_c^2/N;
MS_c = SS_c/1;

%sp D
C_d = D'*y;
SS_d = C_d^2/N;
MS_d = SS_d/1;

%sp CD
C_cd = CD'*y;
SS_cd = C_cd^2/N;
MS_cd = SS_cd/1;

%wp x sp, AC AD BC BD
C_ac = AC'*y;
SS_ac = C_ac^2/N;
MS_ac = SS_ac/1;

C_ad = AD'*y;
SS_ad = C_ad^2/N;
MS_ad = SS_ad/1;

C_bc = BC'*y;
SS_bc = C_bc^2/N;
MS_bc = SS_bc/1;

C_bd = BD'*y;
SS_bd = C_bd^2/N;
MS_bd = SS_bd/1;

%sp error
SS_sp = SST - SS_rep - SS_a - SS_b - SS_ab - SS_repab - SS_c - SS_d - SS_cd - SS_ac - SS_ad - SS_bc - SS_bd; # 15 cost, 17 left for sp error
MSE_sp = SS_sp/17;


%% Testing on wp effects
% check SS code
% CF = (sum(y))^2/32
% RepSS = (sum(y(1:16))^2 + (sum(y(17:32)))^2)/16-CF
disp("--------WP A B AB---------")
disp("Effect | F-stat | p-value")
F_wp = [SS_a/MSE_wp, SS_b/MSE_wp, SS_ab/MSE_wp];
pvalue_wp = [1 - fcdf(SS_a/MSE_wp, 1, 1), 1 - fcdf(SS_b/MSE_wp, 1, 1), 1 - fcdf(SS_ab/MSE_wp, 1, 1)];
effect_wp = [C_a/N/2, C_b/N/2, C_ab/N/2];
disp([effect_wp', F_wp', pvalue_wp']);

disp("--------SP C D CD---------")
disp("Effect | F-stat | p-value")
F_sp = [SS_c/MSE_sp, SS_d/MSE_sp, SS_cd/MSE_sp];
pvalue_sp = [1 - fcdf(SS_c/MSE_sp, 1, 17), 1 - fcdf(SS_d/MSE_sp, 1, 17), 1 - fcdf(SS_cd/MSE_sp, 1, 17)];
effect_sp = [C_c/N/2, C_d/N/2, C_cd/N/2];
disp([effect_sp', F_sp', pvalue_sp']);

disp("--------WP x SP AC AD BC BD---------")
disp("Effect | F-stat | p-value")
F_wpsp = [SS_ac/MSE_sp, SS_ad/MSE_sp, SS_bc/MSE_sp, SS_bd/MSE_sp];
pvalue_wpsp = [1 - fcdf(SS_ac/MSE_sp, 1, 17), 1 - fcdf(SS_ad/MSE_sp, 1, 17), 1 - fcdf(SS_bc/MSE_sp, 1, 17), 1 - fcdf(SS_bd/MSE_sp, 1, 17)];
effect_wpsp = [C_ac/N/2, C_ad/N/2, C_bc/N/2, C_bd/N/2];
disp([effect_wpsp', F_wpsp', pvalue_wpsp']);

% Model specification: A B AB C D AC AD
% Variance components
disp("--------Variance Components---------")
MSE_wp
MSE_sp

%%%%%%%%%%%%%%%%%%%%%
% Do exactly the same thing for all four responses