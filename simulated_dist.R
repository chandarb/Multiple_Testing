emp_dist = function(lambda, tao, sig){
  n = 100000
  comp = lambda * n
  nulls = rnorm(comp, 0, sig)
  signals = rnorm(n - comp, tao, sig)
  draws = c(nulls, signals)
  nume = (1 - lambda) * dnorm(draws, tao, sig)
  den = lambda * dnorm(draws, 0, sig) + nume
  hist(nume / den)
  return(0)
}

emp_dist(.5, 2, 1)