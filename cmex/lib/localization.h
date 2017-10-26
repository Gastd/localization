#ifndef LOCALIZATION_H
#define LOCALIZATION_H

#include <math.h>
#include <stdio.h>

#include "gmatrix.h"
#include "gmatrix_linalg.h"
#include "imu.h"
#include "gps.h"
#include "sonar.h"
#include "magnetometer.h"
#include "rotation.h"

typedef struct{
    int AlgorithmCode;
    int FlagEstimateAccelerometerBias;
    int Nstates;
    PGMATRIX pX;
    PGMATRIX pP;
    PGMATRIX pPreset; // para o caso de reset dos filtros
    PGMATRIX pXsigma_ukf; // para o caso de UKF
    PGMATRIX pWsigma_ukf; // para o caso de UKF
    PGMATRIX pQ;
    PGMATRIX pR_convertedmeasurementtriad;
    PGMATRIX pR_pseudomeasurementnorm;
    PGMATRIX pR_magnetometer; // to use with EKF Decoupled
    PGMATRIX pR_accelerometer; // to use with EKF Decoupled
    PGMATRIX pS_previous_cekf; // para uso com o CEKF
    PGMATRIX pR_previous_cekf; // para uso com o CEKF
    PGMATRIX pA_imu_P_imu_cekf; // para uso com o CEKF
    PGMATRIX pinnovation_previous_cekf; // para uso com o CEKF
} LocalizationFilter;

#define X_q0        1
#define X_q1        2
#define X_q2        3
#define X_q3        4
#define X_px        5
#define X_py        6
#define X_pz        7
#define X_vx        8
#define X_vy        9
#define X_vz        10
#define X_bax       11
#define X_bay       12
#define X_baz       13

#define U_ax        1
#define U_ay        2
#define U_az        3
#define U_wx        4
#define U_wy        5
#define U_wz        6

#define LOCALIZATION_MAXSTATESIZE           13

#define LOCALIZATION_ALGORITHMCODE_EKF2     1
#define LOCALIZATION_ALGORITHMCODE_CEKF     2
#define LOCALIZATION_ALGORITHMCODE_UKF2     3
#define LOCALIZATION_ALGORITHMCODE_EKF_DECOUPLED 4


int localization_init(int AlgorithmCode, int FlagEstimateAccelerometerBias, LocalizationFilter *filter_struct);
int localization_close(LocalizationFilter *filter_struct);
int localization_triad(PQUATERNIONS pq, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG, PQUATERNIONS pq_previous);
int localization_filter_prediction(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T);
int localization_filter_correction(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T);

#endif //LOCALIZATION_H
