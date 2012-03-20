#ifndef IMU_H
#define IMU_H

typedef struct{
	int FlagValidAccerometerMeasure;
	int FlagValidGyrometerMeasure;
    double ax;
    double ay;
    double az;
    double wx;
    double wy;
    double wz;
    double axvariance;
    double ayvariance;
    double azvariance;
    double wxvariance;
    double wyvariance;
    double wzvariance;
} IMUMEASURE, *PIMUMEASURE;

int imu_init(PIMUMEASURE pIMUMeasure);
int imu_close(PIMUMEASURE pIMUMeasure);

#endif //IMU_H
