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
