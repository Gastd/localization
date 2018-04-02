#ifndef MAGNETOMETER_H
#define MAGNETOMETER_H

typedef struct{
	int FlagValidMeasure;
	double mx;
	double my;
	double mz;
	double mxvariance;
	double myvariance;
	double mzvariance;
} MagnetometerMeasure;

int magnetometer_init(MagnetometerMeasure *magnetometer_measure_ptr);
int magnetometer_close(MagnetometerMeasure *magnetometer_measure_ptr);

#endif //MAGNETOMETER_H
