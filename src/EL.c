#include "zeroin.h"
#include <math.h>

double *w1;
double *w2;
static int length1;
static int length2;

double ls1(double lambda)
{
  double sum = 0;
  for (int i = 0; i < length1; i++)
    {
      sum += w1[i] / (1 + lambda * w1[i]);
    }
  return(sum);
}

double ls2(double lambda)
{
  double sum = 0;
  for (int i = 0; i < length2; i++)
    {
      sum += w2[i] / (1 + lambda * w2[i]);
    }
  return(sum);
}

double maxw(int len, double *vals)
{
  double max = vals[0];
  for (int i = 0; i < len; i++)
    {
      if (vals[i] > max)
	max = vals[i];
    }
  return(max);
}

double minw(int len, double *vals)
{
  double min = vals[0];
  for (int i = 0; i < len; i++)
    {
      if (vals[i] < min)
	min = vals[i];
    }
  return(min);
}

void theta_equation(int *len1, double *wval1, double *aval1, 
		   int *len2, double *wval2, double *aval2,
		   double *lambda1, double *lambda2, double *res)
{
  length1 = *len1;
  length2 = *len2;
  w1 = wval1;
  w2 = wval2;

  double rmin = -1 / maxw(length1, wval1) + 1E-7;
  double rmax = -1 / minw(length1, wval1) - 1E-7;
  if (rmin < 0 && rmax > 0)
    *lambda1 = zeroin(rmin, rmax, &ls1, 1E-8);
  else
    {
      if (rmax > 0)
	*res = INFINITY;
      else
	*res = -INFINITY;
      return;
    }

  rmin = -1.0 / maxw(length2, wval2) + 1E-7;
  rmax = -1.0 / minw(length2, wval2) - 1E-7;
  if (rmin < 0 && rmax > 0)
    *lambda2 = zeroin(rmin, rmax, &ls2, 1E-8);
  else
    {
      if (rmax > 0)
	*res = INFINITY;
      else
	*res = -INFINITY;
      return;
    }

  double sum1 = 0;
  for (int i = 0; i < length1; i++)
    {
      sum1 += aval1[i] / (1 + *lambda1 * wval1[i]);
    }
  double sum2 = 0;
  for (int i = 0; i < length2; i++)
    {
      sum2 += aval2[i] / (1 + *lambda2 * wval2[i]);
    }
  *res = *lambda1 * sum1 + *lambda2 * sum2;
}
