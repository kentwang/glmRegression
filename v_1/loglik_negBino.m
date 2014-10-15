function ll = loglik_negBino(mu, y, alpha)
  ll = sum(y .* log(mu ./ (alpha .+ mu)) + alpha .* log(alpha ./ (alpha .+ mu)));