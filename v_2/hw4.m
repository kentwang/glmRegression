data = csvread("hw4_3.csv");
X = [ones(size(data)(1), 1), data(:, 1)];
n = data(:, 2);
y = data(:, 3);