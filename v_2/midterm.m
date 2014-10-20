load("exam1_dataset.mat");
X = [ones(length(x), 1), x];
z0 = (y == 0);