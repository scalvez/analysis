#ifndef STATS_TOOLS_H
#define STATS_TOOLS_H 1

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

int log_factorial(int n)
{
  if (n == 1 || n == 0)
    return 0;

  double s = 0;
  for(int i = 2; i <= n; ++i) {
    s += log(n);
  }
  return s;
}

double Poisson_probability(int n_obs, int n_exp)
{
  return pow(n_exp,n_obs)*exp(-1*n_exp)/factorial(n_obs);
}

#endif
