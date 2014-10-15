X1 = X(:, [1, 3, 4]);
X2 = X(:, [1, 2, 4]);
X3 = X(:, [1, 2, 3]);
reduced_negBino_dev(X1, y, alpha = 0.5)
reduced_negBino_dev(X2, y, alpha = 0.5)
reduced_negBino_dev(X3, y, alpha = 0.5)
reduced_negBino_dev(X, y, alpha = 0.5)