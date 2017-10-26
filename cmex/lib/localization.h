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

enum x_pos
{
    X_q0,
    X_q1,
    X_q2,
    X_q3,
    X_px,
    X_py,
    X_pz,
    X_vx,
    X_vy,
    X_vz,
    X_bax,
    X_bay,
    X_baz
};

enum u_pos
{
    U_ax,
    U_ay,
    U_az,
    U_wx,
    U_wy,
    U_wz
};

#define LOCALIZATION_MAXSTATESIZE           13

enum localization_flag
{
    LOCALIZATION_ALGORITHMCODE_EKF2,
    LOCALIZATION_ALGORITHMCODE_CEKF,
    LOCALIZATION_ALGORITHMCODE_UKF2,
    LOCALIZATION_ALGORITHMCODE_EKF_DECOUPLED
};

int localization_init(int AlgorithmCode, int FlagEstimateAccelerometerBias, LocalizationFilter *filter_struct);
int localization_close(LocalizationFilter *filter_struct);
int localization_triad(PQUATERNIONS pq, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG, PQUATERNIONS pq_previous);
int localization_filter_prediction(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T);
int localization_filter_correction(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T);

#endif //LOCALIZATION_H
