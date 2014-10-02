function f = plot_convergence(M)
  n = dim(M)(2)
  ncol_plot = ceil(sqrt(n));
  nrow_plot = ceil(n / ncol_plot)
  for i in 1:dim(M)(2)
    
  subplot(3, 2, 1); plot(M(:, 1)); title('\beta_0'); xlabel("Iteration");
  subplot(3, 2, 2); plot(M(:, 2)); title('\beta_1'); xlabel("Iteration");
  subplot(3, 2, 3); plot(M(:, 3)); title('\beta_2'); xlabel("Iteration");
  subplot(3, 2, 4); plot(M(:, 4)); title('\beta_3'); xlabel("Iteration");
  subplot(3, 2, 5:6); plot(M(:, 5)); title('Deviance: \lambda(\beta)'); xlabel("Iteration");
  suptitle("NB regression convergence using OLS as starting value for IRWLS");