#ifndef STATS_TOOLS_H
#define STATS_TOOLS_H 1

#include "TMath.h"

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

double p0_S(double S)
{
  double S0 = 82;
  return S<S0 ? 1/S0 : 0;
}

double p0_B(double B)
{
  double B0 = 10;
  double mu_B = B0;
  double sigma = B0/2;
  return B >=0 ? TMath::Gaus(B,mu_B,sigma,true) : 0;
}

double proba_bayes(double S, double B,
                   std::vector<double> pdf_S, std::vector<double> pdf_B, std::vector<double> data) {

  double p = 1;

  for(unsigned int i = 0; i<pdf_S.size(); ++i) {
    double lambda_i = S*pdf_S[i] + B*pdf_B[i];
    p *= pow(lambda_i,data[i])*exp(-lambda_i)/factorial(int(data[i]));
  }

  return p;
}

#endif
