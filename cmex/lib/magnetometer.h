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
} MAGNETOMETERMEASURE, *PMAGNETOMETERMEASURE;

int magnetometer_init(PMAGNETOMETERMEASURE pMagnetometerMeasure);
int magnetometer_close(PMAGNETOMETERMEASURE pMagnetometerMeasure);

#endif //MAGNETOMETER_H
