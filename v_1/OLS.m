function beta = OLS(X, y)
	beta = inverse(X' * X) * X' * y;