function f = plot_convergence(M)
  n = size(M)(2);
  ncol_plot = 2;
  nrow_plot = ceil(n / ncol_plot);
  
  %- Reserved first graph for deviance
  figure;
  subplot(nrow_plot, ncol_plot, 1:2);
  plot(M(:, 1));
  title('Deviance: \lambda(\beta)');
  xlabel("Iteration");
  
  for i = 2:nrow_plot
    for j = 1:ncol_plot
      ind = (i - 1) * ncol_plot + j - 1;
      if ind - 1 <= n
        subplot(nrow_plot, ncol_plot, ind + 1);
        plot(M(:, ind));
        title(strcat('\beta_', int2str(ind - 2)));
        xlabel("Iteration");
      else
        break;
      endif
    endfor
  endfor
  suptitle("Convergence of IRWLS");
