#include <cpplapack.h>
#include <gsl/gsl_sf_expint.h>
#include <SE_trans.h>
#include <time.h>

double F(double x)
{
  return sqrt(cos(4 * atanh(x)) + cosh(M_PI));
}

double ex1_u(double x)
{
  return sin((1-x*x)*F(x));
}

double ex2_u(double x)
{
  return cos((1-x*x)*F(x));
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

double xSE(double a, double b, int k, double h)
{
  return SE_trans(a, b, k*h);
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
  double a =-1.0;
  double b = 1.0;
  double d = 1.57;
  double tmp;
  clock_t start, end;
  double time;

  for (int N = 1; N <= 140; N += 8) {
  start = clock();

  int m = (2*N+1);
  int n = (2*N+1)*2;
  double h = sqrt(M_PI*d/N);

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
      tmp = SE_trans(a, b, j*h);
      A_N(i+N+m, j+N)   = - 2 * h * SE_trans_div(a, b, j*h) * del1(i, j)
            * (tmp * F(tmp) * F(tmp) + sin(2*j*h))/ F(tmp);
      A_N(i+N  , j+N+m) =   2 * h * SE_trans_div(a, b, j*h) * del1(i, j)
            * (tmp * F(tmp) * F(tmp) + sin(2*j*h))/ F(tmp);
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

  for (int i=-(SAMPLE/2)+1; i< (SAMPLE/2); i++) {
    double err;
    double ans = 0;
    double x = i*hh;

    for (int j=N; j>0; j--) {
      tmp = SE_trans(a, b, j*h);
      ans += y_2( j+N) * SE_trans_div(a, b, j*h)
                       * (-2*(tmp * F(tmp) * F(tmp) + sin( 2*j*h))/ F(tmp))
                       * J( j, h, SE_trans_inv(a, b, x));
      tmp = SE_trans(a, b,-j*h);
      ans += y_2(-j+N) * SE_trans_div(a, b,-j*h)
                       * (-2*(tmp * F(tmp) * F(tmp) + sin(-2*j*h))/ F(tmp))
                       * J(-j, h, SE_trans_inv(a, b, x));
    }
      tmp = SE_trans(a, b, 0*h);
      ans += y_2( 0+N) * SE_trans_div(a, b, 0*h)
                       * (-2*(tmp * F(tmp) * F(tmp) + sin( 2*0*h))/ F(tmp))
                       * J( 0, h, SE_trans_inv(a, b, x));
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

    for (int j=N; j>0; j--) {
      tmp = SE_trans(a, b, j*h);
      ans += 2 * y_1( j+N) * SE_trans_div(a, b, j*h)
                       * ((tmp * F(tmp) * F(tmp) + sin( 2*j*h))/ F(tmp))
                       * J( j, h, SE_trans_inv(a, b, x));
      tmp = SE_trans(a, b,-j*h);
      ans += 2 * y_1(-j+N) * SE_trans_div(a, b,-j*h)
                       * ((tmp * F(tmp) * F(tmp) + sin(-2*j*h))/ F(tmp))
                       * J(-j, h, SE_trans_inv(a, b, x));
    }
      tmp = SE_trans(a, b, 0*h);
      ans += 2 * y_1( 0+N) * SE_trans_div(a, b, 0*h)
                       * ((tmp * F(tmp) * F(tmp) + sin( 2*0*h))/ F(tmp))
                       * J( 0, h, SE_trans_inv(a, b, x));
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
