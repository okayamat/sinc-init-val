#include <cpplapack.h>
#include <gsl/gsl_sf_expint.h>
#include <DE_trans.h>
#include <time.h>

double ex1_u(double x)
{
  return sqrt(x)*exp(-x);
}

double ex2_u(double x)
{
  return exp(-x);
}

double S(int j, double h, double x)
{
  double val = M_PI*(x - j*h)/h;
  if (val == 0)
    return 1.0;
  else
    return sin(val)/val;
}

double J(int j, double h, double x)
{
  return h * (0.5 + M_1_PI * gsl_sf_Si(M_PI * (x / h - j)));
}

double xDE(double a, double b, int k, double h)
{
  return DE_trans(a, b, k*h);
}

double del0(int i, int j)
{
  if (i == j)
    return 1.0;
  else
    return 0.0;
}

double del1(int i, int j)
{
  return 0.5 + M_1_PI * gsl_sf_Si(M_PI * (i - j));
}

int main()
{
  double a = 0.0;
  double b = 2.0;
  double d = 1.57;
  clock_t start, end;
  double time;

  for (int N = 1; N <= 70; N += 8) {
  start = clock();

  int m = (2*N+1);
  int n = (2*N+1)*2;
  double h = log(2 * d * N/(0.5)) / N;

  CPPL::dgematrix A_N(n,n);
  CPPL::dcovector f_N(n);
  CPPL::dcovector y_1(m);
  CPPL::dcovector y_2(m);

  for (int i = -N; i <= N; i++) {
    for (int j = -N; j <= N; j++) {
      A_N(i+N, j+N)     = del0(i+N, j+N)
                         + h * DE_trans_div(a, b, j*h) * del1(i, j);
      A_N(i+N  , j+N+m) = - 0.5 * h * del1(i, j) * M_PI * cosh(j*h) * sqrt(2)
              / (sqrt(1 + exp(-M_PI*sinh(j*h)))*(1 + exp(M_PI*sinh(j*h))));
      A_N(i+N+m, j+N)   =   h * del1(i, j) * M_PI * cosh(j*h) * sqrt(2)
              / (sqrt(1 + exp(-M_PI*sinh(j*h)))*(1 + exp(M_PI*sinh(j*h))));
      A_N(i+N+m, j+N+m) = del0(i+N+m, j+N+m);
    }
  }

  for (int i = 0; i < m; i++) {
    f_N(i)   = 0;
    f_N(i+m) = 1;
  }

  A_N.dgesv(f_N);

  for (int i = 0; i < m; i++) {
    y_1(i) = f_N(i);
    y_2(i) = f_N(i+m);
  }

  int SAMPLE = 1000;
  double hh = (b-a)/SAMPLE;
  double maxerr = 0;

  for (int i=1; i< SAMPLE; i++) {
    double err;
    double ans = 0;
    double x = i*hh;

    for (int j=N; j>=-N; j--) {
      ans += (0.5 * y_2( j+N)/sqrt(DE_trans(a, b, j*h)) - y_1(j+N))
            * DE_trans_div(a, b, j*h) * J( j, h, DE_trans_inv(a, b, x));
    }
    ans += 0;

    err = fabs(ex1_u(x) - ans);
    //    std::cout << x << "\t" << err << std::endl;

    if (maxerr < err)
      maxerr = err;
  }

  for (int i=1; i< SAMPLE; i++) {
    double err;
    double ans = 0;
    double x = i*hh;

    for (int j=N; j>= -N; j--) {
      ans += - y_1( j+N) * DE_trans_div(a, b, j*h)
             * J( j, h, DE_trans_inv(a, b, x))/sqrt(DE_trans(a, b, j*h));
    }
    ans += 1;

    err = fabs(ex2_u(x) - ans);
    //    std::cout << x << "\t" << err << std::endl;

    if (maxerr < err)
      maxerr = err;
  }

  end = clock();

  time = (double)(end - start) / CLOCKS_PER_SEC;

  std::cout << N << "\t" << time << "\t" << maxerr << std::endl;
  }
}
