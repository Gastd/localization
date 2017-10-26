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
} ImuMeasure;

int imu_init(ImuMeasure *imu_measure_ptr);
int imu_close(ImuMeasure *imu_measure_ptr);

#endif //IMU_H
