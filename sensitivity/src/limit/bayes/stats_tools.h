#ifndef STATS_TOOLS_H
#define STATS_TOOLS_H 1

#include "TMath.h"

// int factorial(int n)
// {
//   return std::factorial(n);
// }

int factorial(int n)
{
  if(n>=32)
    return 1e9;

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
  double B0 = 50;
  double mu_B = B0;
  double sigma = 7;
  // double sigma = B0/2;
  return B >=0 ? TMath::Gaus(B,mu_B,sigma,true) : 0;
}

double proba_bayes(double S, double B,
                   std::vector<double> pdf_S, std::vector<double> pdf_B, std::vector<double> data) {

  double p = 1;

  // std::cout << "PROBA BAYES" << std::endl;

  for(unsigned int i = 0; i<pdf_S.size(); ++i) {
    double lambda_i = S*pdf_S[i] + B*pdf_B[i];
    // std::cout << "------- " << i << std::endl;

    // std::cout << "lambda " << lambda_i << std::endl;

    // std::cout << "pow(lambda_i,data[i]) " << pow(lambda_i,data[i]) << std::endl;
    // std::cout << "exp(-lambda_i)/factorial(int(data[i])) " << exp(-lambda_i)/factorial(int(data[i])) << std::endl;

    p *= pow(lambda_i,data[i])*exp(-lambda_i)/factorial(int(data[i]));

    // std::cout << "p " << p << std::endl;

    if(std::isinf(p) || p<1e-300)
      return 0;

    // if(S<=1 && B<=1)
    //   std::cout << "p " << p << std::endl;

    // std::cout << "--  i proba bayes "<< i << "  " << pow(lambda_i,data[i])*exp(-lambda_i)/factorial(int(data[i])) << std::endl;
  }


  return p;
}

#endif
