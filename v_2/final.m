load("250data.mat");

y = [RunHH1, RunHH2, RunHL1, RunHL2, RunLH1, RunLH2, RunLL1, RunLL2]';
S = length(y);
p = 3; % dimension of fix effects. Constant, rake angle and cutting speed
n = 250; % observations for each run, balanced design
N  = 8; % total number of runs
X = [repmat([1, 1], 2*n, 1); repmat([1, -1], 2*n, 1); repmat([-1, 1], 2*n, 1); repmat([-1, -1], 2*n, 1)];
X = [repmat(1, S, 1), X];


