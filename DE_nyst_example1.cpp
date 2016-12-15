#include <cpplapack.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_airy.h>
#include <DE_trans.h>
#include <time.h>

double ex1_u(double x)
{
  return gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE);
}

double ex2_u(double x)
{
  return gsl_sf_airy_Ai_deriv(x, GSL_PREC_DOUBLE);
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
  double h = log(2 * d * N) / N;

  CPPL::dgematrix A_N(n,n);
  CPPL::dcovector f_N(n);
  CPPL::dcovector y_1(m);
  CPPL::dcovector y_2(m);

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      A_N(i, j)     = del0(i, j);
      A_N(m+i, m+j) = del0(m+i, m+j);
    }
  }

  for (int i = -N; i <= N; i++) {
    for (int j = -N; j <= N; j++) {
      A_N(i+N+m, j+N)   = - h * DE_trans_div(a, b, j*h) * DE_trans(a, b, j*h) * del1(i, j);
      A_N(i+N  , j+N+m) = - h * DE_trans_div(a, b, j*h) * del1(i, j);
    }
  }

  for (int i = 0; i < m; i++) {
    f_N(i)   = 0.35502805388781724;
    f_N(i+m) =-0.25881940379280680;
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

    for (int j=N; j>0; j--) {
      ans += y_2( j+N) * DE_trans_div(a, b, j*h)
                       * J( j, h, DE_trans_inv(a, b, x));
      ans += y_2(-j+N) * DE_trans_div(a, b,-j*h)
                       * J(-j, h, DE_trans_inv(a, b, x));
    }
      ans += y_2( 0+N) * DE_trans_div(a, b, 0*h)
                       * J( 0, h, DE_trans_inv(a, b, x));
    ans += 0.35502805388781724;

    err = fabs(ex1_u(x) - ans);
    //    std::cout << x << "\t" << err << std::endl;

    if (maxerr < err)
      maxerr = err;
  }

  for (int i=1; i< SAMPLE; i++) {
    double err;
    double ans = 0;
    double x = i*hh;

    for (int j=N; j>0; j--) {
      ans += y_1( j+N) * DE_trans_div(a, b, j*h) * DE_trans(a, b, j*h)
                       * J( j, h, DE_trans_inv(a, b, x));
      ans += y_1(-j+N) * DE_trans_div(a, b,-j*h) * DE_trans(a, b,-j*h)
                       * J(-j, h, DE_trans_inv(a, b, x));
    }
      ans += y_1( 0+N) * DE_trans_div(a, b, 0*h) * DE_trans(a, b, 0*h)
                       * J( 0, h, DE_trans_inv(a, b, x));
    ans += -0.25881940379280680;

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
