#include <stdlib.h>
#include <math.h>

#ifndef STDDEV_H
#define STDDEV_H

double mean (double *data, int n)
{
	double mean = 0.0;
	int i;
   #ifdef _OPENMP
   #pragma omp parallel shared(mean) private(i)
   #pragma omp for reduction(+:mean)
   #endif 
	for (i=0; i<n; i++) mean += data[i];
	mean /= n;
	return mean;
}

double standard_deviation (double *data, double m, int n)
{
	double std_dev = 0.0;
	int i;
   #ifdef _OPENMP
   #pragma omp parallel shared(std_dev, data) private(i)
   #pragma omp for reduction(+:std_dev)
   #endif 
	for (i=0; i<n; i++) std_dev += (data[i]-m) * (data[i]-m);
	std_dev = sqrt(std_dev / n);
	return std_dev;
}

double drift (double *data, int n)
{
	double SUM_x=0.0, SUM_xx=0.0;		// x data counts from 0 to n-1
	double SUM_y=0.0, SUM_xy=0.0;

	int i;
   #ifdef _OPENMP
   #pragma omp parallel for reduction(+:SUM_x)
   #endif 
	for (i=0; i<n; i++) SUM_x  += i*1.0/n;
   #ifdef _OPENMP
   #pragma omp parallel for reduction(+:SUM_xx)
   #endif 
	for (i=0; i<n; i++) SUM_xx += i*i*1.0/n;
   #ifdef _OPENMP
   #pragma omp parallel for reduction(+:SUM_y)
   #endif 
	for (i=0; i<n; i++) SUM_y  += data[i]/n;
   #ifdef _OPENMP
   #pragma omp parallel for reduction(+:SUM_xy)
   #endif 
	for (i=0; i<n; i++) SUM_xy += i*data[i]/n;

	double slope = (SUM_x/n*SUM_y - SUM_xy/n)/(SUM_x/n*SUM_x - SUM_xx/n);
//	double y-intercept = SUM_y - slope*SUM_x;
	double drift = (n-1)*slope;

	return drift;
}

#endif

