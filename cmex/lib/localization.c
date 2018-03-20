/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: localization.c 
Conteudo: Funções em código C relacionadas à localização.
Autor: G. A. Borges.
Atualizações: 
    - 19/09/2008: criação por Geovany A. Borges, modulo sem init/close
*****************************************************************************/
#include "localization.h"
#include "kalman.h"
#include "unscentedtransform.h"
#include "filter_parameters.h"

// Definicoes de uso interno

// Prototipos de funcoes internas ao modulo
int localization_filter_prediction_cekf(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T);
int localization_filter_correction_cekf(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T);
int localization_filter_prediction_ekf2(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T);
int localization_filter_correction_ekf2(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T);
int localization_filter_prediction_ekf_decoupled(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T);
int localization_filter_correction_ekf_decoupled(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T);
int localization_filter_process_model_evaluate(PGMATRIX pX, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias);
int localization_filter_process_model_df_dx(PGMATRIX pdf_dx, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias);
int localization_filter_process_model_df_du(PGMATRIX pdf_du, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias);
int localization_converted_measurement_triad(PGMATRIX pYm, PGMATRIX pPy, PGMATRIX pX_predicted, PGMATRIX pP_predicted, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias);
int localization_filter_observation_model_triad_evaluate(PGMATRIX pYm, LocalizationFilter *filter_struct, PGMATRIX pX_predicted, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG);
void localization_filter_quaternionscorrectsign(LocalizationFilter *filter_struct, PGMATRIX pX_predicted);
int localization_converted_measurement_triad_dg_du_imu(PGMATRIX pdg_du_imu, PGMATRIX pX_predicted, PGMATRIX pP_predicted, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias);

// Variaveis globais do módulo

/*****************************************************************************
******************************************************************************
** FUNCOES COM CHAMADA EXTERNA
******************************************************************************
*****************************************************************************/

/*****************************************************************************
*** int localization_init(int AlgorithmCode, int FlagEstimateAccelerometerBias, PLOCALIZATIONFILTERSTRUCT pFilterStruct)
*** Entradas: 
*** Saidas:
*****************************************************************************/
int localization_init(int AlgorithmCode, int FlagEstimateAccelerometerBias, LocalizationFilter *filter_struct)
{
    // Iniciar componentes da estrutura do filtro
    filter_struct->AlgorithmCode = AlgorithmCode;

    filter_struct->FlagEstimateAccelerometerBias = FlagEstimateAccelerometerBias;

    filter_struct->Nstates = 4 + 3 + 3; // quaternions + 3D position + 3D speed
    if(FlagEstimateAccelerometerBias)
    {
        filter_struct->Nstates += 3; // + accelerometer bias
    }

    filter_struct->pX = PGMATRIX_ALLOC(filter_struct->Nstates, 1); // %should be set by user.
    filter_struct->pP = PGMATRIX_ALLOC(filter_struct->Nstates, filter_struct->Nstates); // %should be set by user.
    filter_struct->pPreset = PGMATRIX_ALLOC(filter_struct->Nstates, filter_struct->Nstates); // %should be set by user.

    if(filter_struct->AlgorithmCode == LOCALIZATION_ALGORITHMCODE_UKF2)
    {
        filter_struct->pXsigma_ukf = PGMATRIX_ALLOC(filter_struct->Nstates, 2*filter_struct->Nstates+1);
        filter_struct->pWsigma_ukf = PGMATRIX_ALLOC(1, 2*filter_struct->Nstates+1);
    }

    filter_struct->pQ = PGMATRIX_ALLOC(filter_struct->Nstates, filter_struct->Nstates);
    PGMATRIX_ZEROES(filter_struct->pQ);
    PGMATRIX_DATA(filter_struct->pQ, 1, 1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_Q0);
    PGMATRIX_DATA(filter_struct->pQ, 2, 2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_Q1);
    PGMATRIX_DATA(filter_struct->pQ, 3, 3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_Q2);
    PGMATRIX_DATA(filter_struct->pQ, 4, 4) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_Q3);

    PGMATRIX_DATA(filter_struct->pQ, 5, 5) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_POSITION_X);
    PGMATRIX_DATA(filter_struct->pQ, 6, 6) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_POSITION_Y);
    PGMATRIX_DATA(filter_struct->pQ, 7, 7) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_POSITION_Z);

    PGMATRIX_DATA(filter_struct->pQ, 8, 8) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_VELOCITY_X);
    PGMATRIX_DATA(filter_struct->pQ, 9, 9) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_VELOCITY_Y);
    PGMATRIX_DATA(filter_struct->pQ, 10, 10) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_VELOCITY_Z);

    if(filter_struct->FlagEstimateAccelerometerBias)
    {
        PGMATRIX_DATA(filter_struct->pQ, 11, 11) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_ACCELEROMETER_BIAS_X);
        PGMATRIX_DATA(filter_struct->pQ, 12, 12) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_ACCELEROMETER_BIAS_Y);
        PGMATRIX_DATA(filter_struct->pQ, 13, 13) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_Q_MATRIX_ACCELEROMETER_BIAS_Z);
    }

    filter_struct->pR_pseudomeasurementnorm = PGMATRIX_ALLOC(1, 1);
    PGMATRIX_DATA(filter_struct->pR_pseudomeasurementnorm, 1, 1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_PSEUDOMEASUREMENTNORM);

    filter_struct->pR_convertedmeasurementtriad = PGMATRIX_ALLOC(4, 4);
    PGMATRIX_ZEROES(filter_struct->pR_convertedmeasurementtriad);

    filter_struct->pR_accelerometer = PGMATRIX_ALLOC(3, 3);
    PGMATRIX_ZEROES(filter_struct->pR_accelerometer);

    filter_struct->pR_magnetometer = PGMATRIX_ALLOC(3, 3);
    PGMATRIX_ZEROES(filter_struct->pR_magnetometer);

    if(filter_struct->AlgorithmCode == LOCALIZATION_ALGORITHMCODE_CEKF)
    {
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 1, 1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_1_1);
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 2, 2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_2_2);
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 3, 3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_3_3);
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 4, 4) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_CEKF_4_4);
    }
    else
    {
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 1, 1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_1_1);
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 2, 2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_2_2);
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 3, 3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_3_3);
        PGMATRIX_DATA(filter_struct->pR_convertedmeasurementtriad, 4, 4) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_CONVERTEDMEASUREMENTTRIAD_4_4);
    }

    if(filter_struct->AlgorithmCode == LOCALIZATION_ALGORITHMCODE_CEKF)
    {
        filter_struct->pS_previous_cekf = PGMATRIX_ALLOC(filter_struct->Nstates, 4);
        PGMATRIX_ZEROES(filter_struct->pS_previous_cekf);
        filter_struct->pR_previous_cekf = PGMATRIX_ALLOC(4, 4);
        PGMATRIX_IDENTITY(filter_struct->pR_previous_cekf);
        filter_struct->pA_imu_P_imu_cekf = PGMATRIX_ALLOC(filter_struct->Nstates, 6);
        PGMATRIX_ZEROES(filter_struct->pA_imu_P_imu_cekf);
        filter_struct->pinnovation_previous_cekf = PGMATRIX_ALLOC(4, 1);
        PGMATRIX_ZEROES(filter_struct->pinnovation_previous_cekf);
    }

    if(filter_struct->AlgorithmCode == LOCALIZATION_ALGORITHMCODE_EKF_DECOUPLED)
    {
        PGMATRIX_DATA(filter_struct->pR_accelerometer, 1, 1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_1_1);
        PGMATRIX_DATA(filter_struct->pR_accelerometer, 2, 2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_2_2);
        PGMATRIX_DATA(filter_struct->pR_accelerometer, 3, 3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_3_3);

        PGMATRIX_DATA(filter_struct->pR_magnetometer, 1, 1) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_MAGNETOMETER_EKF_DECOUPLED_1_1);
        PGMATRIX_DATA(filter_struct->pR_magnetometer, 2, 2) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_MAGNETOMETER_EKF_DECOUPLED_2_2);
        PGMATRIX_DATA(filter_struct->pR_magnetometer, 3, 3) = GMATRIXMACRO_SQR(FILTER_PARAMETERS_R_MAGNETOMETER_EKF_DECOUPLED_3_3);
    }

    return 1; 
}

/*****************************************************************************
*** int localization_close(PLOCALIZATIONFILTERSTRUCT pFilterStruct)
*** Entradas: 
*** Saidas:
*****************************************************************************/
int localization_close(LocalizationFilter *filter_struct)
{
    // Desalocar componentes da estrutura do filtro
    PGMATRIX_FREE(filter_struct->pX);
    PGMATRIX_FREE(filter_struct->pP);
    PGMATRIX_FREE(filter_struct->pPreset);

    if(filter_struct->AlgorithmCode == LOCALIZATION_ALGORITHMCODE_UKF2)
    {
        PGMATRIX_FREE(filter_struct->pXsigma_ukf);
        PGMATRIX_FREE(filter_struct->pWsigma_ukf);
    }

    PGMATRIX_FREE(filter_struct->pQ);

    PGMATRIX_FREE(filter_struct->pR_accelerometer);
    PGMATRIX_FREE(filter_struct->pR_magnetometer);
    PGMATRIX_FREE(filter_struct->pR_pseudomeasurementnorm);

    PGMATRIX_FREE(filter_struct->pR_convertedmeasurementtriad);

    if(filter_struct->AlgorithmCode == LOCALIZATION_ALGORITHMCODE_CEKF)
    {
        PGMATRIX_FREE(filter_struct->pS_previous_cekf);
        PGMATRIX_FREE(filter_struct->pR_previous_cekf);
        PGMATRIX_FREE(filter_struct->pA_imu_P_imu_cekf);
        PGMATRIX_FREE(filter_struct->pinnovation_previous_cekf);
    }

    // Retorna
    return 1; 
}

/*****************************************************************************
*** int localization_filter_correction(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGPSMEASURE pGPSMeasure, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PSONARMEASURE pSonarMeasure, PGMATRIX pM, PGMATRIX pG, double T)
*** Entradas: 
*** Saidas:
*****************************************************************************/
int localization_filter_correction(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T)
{
    // Corre��o conforme a arquitetura
    switch(filter_struct->AlgorithmCode)
    {
    case LOCALIZATION_ALGORITHMCODE_EKF2:
        if(!localization_filter_correction_ekf2(filter_struct, gps_measure_ptr, imu_measure_ptr, magnetometer_measure_ptr, sonar_measure_ptr, pM, pG, T)) return 0;
        break;
    case LOCALIZATION_ALGORITHMCODE_CEKF:
        if(!localization_filter_correction_cekf(filter_struct, gps_measure_ptr, imu_measure_ptr, magnetometer_measure_ptr, sonar_measure_ptr, pM, pG, T)) return 0;
        break;
    case LOCALIZATION_ALGORITHMCODE_EKF_DECOUPLED:
        if(!localization_filter_correction_ekf_decoupled(filter_struct, gps_measure_ptr, imu_measure_ptr, magnetometer_measure_ptr, sonar_measure_ptr, pM, pG, T)) return 0;
        break;
    default:
        return 0;
    }

    // Reset em caso de matriz P negativa:
/*  if(sum(eig(kf_structure.P)<0)>0)
        % filter reset:
        disp('*** Warning: Kalman filter reset');
        kf_structure.P = kf_structure.Preset;
    end
*/
    // Retorna
    return 1; 
}


void localization_filter_quaternionscorrectsign(LocalizationFilter *filter_struct, PGMATRIX pX_predicted)
{

    QUATERNIONS_DECLARE(qcurrent);
    QUATERNIONS_DECLARE(qpredicted);

    QUATERNIONS_Q0(qcurrent) = PGMATRIX_DATA(filter_struct->pX,X_q0, 1);
    QUATERNIONS_Q1(qcurrent) = PGMATRIX_DATA(filter_struct->pX,X_q1, 1);
    QUATERNIONS_Q2(qcurrent) = PGMATRIX_DATA(filter_struct->pX,X_q2, 1);
    QUATERNIONS_Q3(qcurrent) = PGMATRIX_DATA(filter_struct->pX,X_q3, 1);
    QUATERNIONS_Q0(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q0, 1);
    QUATERNIONS_Q1(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q1, 1);
    QUATERNIONS_Q2(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q2, 1);
    QUATERNIONS_Q3(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q3, 1);
    rotation_quaternionscorrectsign(&qcurrent, &qpredicted);
    PGMATRIX_DATA(filter_struct->pX,X_q0,1) = QUATERNIONS_Q0(qcurrent);
    PGMATRIX_DATA(filter_struct->pX,X_q1,1) = QUATERNIONS_Q1(qcurrent);
    PGMATRIX_DATA(filter_struct->pX,X_q2,1) = QUATERNIONS_Q2(qcurrent);
    PGMATRIX_DATA(filter_struct->pX,X_q3,1) = QUATERNIONS_Q3(qcurrent);
}

int localization_filter_correction_cekf(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T)
{
// %%% corrige usando medidas do algoritmo TRIAD: 
// %%% Fun�ao de medi�ao: [eye(4) zeros(4,6)]*X (a medida � dada pelo estimador TRIAD.
// %%% Estimar a medi�ao do quaternion assim como da matriz R usando UT;

    GMATRIX_DECLARE(dg_du_imu, 4, 6);
    GMATRIX_DECLARE(H,LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(V, 4, 1);
    GMATRIX_DECLARE(Ym, 4, 1);
    GMATRIX_DECLARE(Py, 4, 4);
    GMATRIX_DECLARE(Rtilde, 4, 4);
    GMATRIX_DECLARE(Stilde, LOCALIZATION_MAXSTATESIZE, 4);
    GMATRIX_DECLARE(R, 4, 4);
    GMATRIX_DECLARE(X_predicted, LOCALIZATION_MAXSTATESIZE, 1);
    GMATRIX_DECLARE(P_predicted, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy1, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy2, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy3, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy4, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    DUMMY_MATRICES MatDummy;

    double aux;

    MatDummy.pMat1 = &MatDummy1;
    MatDummy.pMat2 = &MatDummy2;
    MatDummy.pMat3 = &MatDummy3;
    MatDummy.pMat4 = &MatDummy4;

    // Defini��o dos tamanhos das matrizes deve ser feito aqui, o que evita aloca��o dim�mica
    GMATRIX_SETSIZE(X_predicted, filter_struct->Nstates, 1);
    GMATRIX_SETSIZE(P_predicted, filter_struct->Nstates, filter_struct->Nstates);
    GMATRIX_SETSIZE(Stilde, filter_struct->Nstates, 1);

    // Medi��o pelo TRIAD
    if(magnetometer_measure_ptr->FlagValidMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);

        // if kf_structure.flagestimateaccelerometerbias
        //     H = [eye(4) zeros(4,9)];
        // else
        //     H = [eye(4) zeros(4,6)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 4, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 4, 10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H, 1, 1) = 1.0;
        GMATRIX_DATA(H, 2, 2) = 1.0;
        GMATRIX_DATA(H, 3, 3) = 1.0;
        GMATRIX_DATA(H, 4, 4) = 1.0;

        // [Ym,Py] = ConvertedMeasurementTRIAD('cekf',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
        GMATRIX_SETSIZE(Ym, 4, 1);
        GMATRIX_SETSIZE(Py, 4, 4);
        localization_converted_measurement_triad(&Ym, &Py, &X_predicted, &P_predicted, imu_measure_ptr, magnetometer_measure_ptr, pM, pG, filter_struct->FlagEstimateAccelerometerBias);
//      GMATRIX_PRINT_MATLABFORM(Ym);
//      GMATRIX_PRINT_MATLABFORM(Py);

        // v = (Ym - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, &Ym, MatDummy.pMat1);

        // Rtilde = Py + kf_structure.R_convertedmeasurementtriad_cekf;
        PGMATRIX_ADD_COPY(&Rtilde, &Py, filter_struct->pR_convertedmeasurementtriad);

        // Stilde = kf_structure.A_imu_P_imu_cekf * localization_filter_model_gtriad('dg_du_imu',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias,1)';
        localization_converted_measurement_triad_dg_du_imu(&dg_du_imu, &X_predicted, &P_predicted, imu_measure_ptr, magnetometer_measure_ptr, pM, pG, filter_struct->FlagEstimateAccelerometerBias);
        PGMATRIX_MULTIPLY_COPY_EXTENDED(&Stilde, filter_struct->pA_imu_P_imu_cekf, 0, &dg_du_imu, 1);

        // K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad_ekf);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_ADD_COPY(&R, &Py, filter_struct->pR_convertedmeasurementtriad);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);

        // kf_structure.R_previous_cekf = Rtilde;
        PGMATRIX_COPY(filter_struct->pR_previous_cekf, &Rtilde);

        // kf_structure.S_previous_cekf = Stilde;
        PGMATRIX_COPY(filter_struct->pS_previous_cekf, &Stilde);

        // kf_structure.innovation_previous_cekf = (Ym - H*kf_structure.X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, filter_struct->pX);
        PGMATRIX_SUBTRACT_COPY(&V, &Ym, MatDummy.pMat1);
        PGMATRIX_COPY(filter_struct->pinnovation_previous_cekf, &V);
    }

    // %%% corrige usando medidas do GPS:
    // %%% Fun�ao de medi�ao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
    if(gps_measure_ptr->FlagValidPositionMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
        // else
        //     H = [zeros(3,4) eye(3,3) zeros(3,3)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H,3,13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H,3,10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H,1,5) = 1.0;
        GMATRIX_DATA(H,2,6) = 1.0;
        GMATRIX_DATA(H,3,7) = 1.0;
        // v = (gpsmeasure.p - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, gps_measure_ptr->pPosition, MatDummy.pMat1);

        // K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_COPY(&R,gps_measure_ptr->pPPosition);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do GPS:
    // %%% Fun�ao de medi�ao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
    if(gps_measure_ptr->FlagValidVelocityMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
        // else
        //     H = [zeros(3,4) zeros(3,3) eye(3,3)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 3, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 3, 10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H, 1, 8) = 1.0;
        GMATRIX_DATA(H, 2, 9) = 1.0;
        GMATRIX_DATA(H, 3, 10) = 1.0;
        // v = (gpsmeasure.v - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, gps_measure_ptr->pVelocity, MatDummy.pMat1);

        // K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_COPY(&R,gps_measure_ptr->pPVelocity);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do sonar:
    // %%% Fun�ao de medi�ao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
    // %%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
    if(sonar_measure_ptr->FlagValidMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // dHdX1 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(1);
        // dHdX2 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(2);
        // dHdX3 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(3);
        // dHdX4 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(4);
        // dHdX7 =  1/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0];
        // else
        //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H,1,13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H,1,10);
            GMATRIX_ZEROES(H);
        }

        aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,1,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,2,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,3,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,4,1));

        if(aux != 0.0)
        { // situa��o em que a medi��o seria diferente de infinito
            GMATRIX_DATA(H,1,1) = -GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,1,1);
            GMATRIX_DATA(H,1,2) =  GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,2,1);
            GMATRIX_DATA(H,1,3) =  GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,3,1);
            GMATRIX_DATA(H,1,4) = -GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,4,1);
            GMATRIX_DATA(H,1,7) =  1.0/(aux);

            // v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
            GMATRIX_SETSIZE(V,1,1);
            GMATRIX_DATA(V,1,1) = sonar_measure_ptr->range - GMATRIX_DATA(X_predicted,7,1) / aux;

            // K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
            // kf_structure.X = X + K*v;
            // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
            GMATRIX_SETSIZE(R,1,1);
            GMATRIX_DATA(R,1,1) = sonar_measure_ptr->rangevariance;
            kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

            // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
            localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
        }
    }

    // %%% corrige a norma do quaternion usando pseudo-medi�oes:
    // %%% Fun�ao de medi�ao: 1 = q0^2 + q1^2 + q2^2 + q3^2
    // %%% Fun�ao de medi�ao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
    // X = kf_structure.X;
    PGMATRIX_COPY(&X_predicted, filter_struct->pX);

    // P = kf_structure.P;
    PGMATRIX_COPY(&P_predicted, filter_struct->pP);
    // dHdX1 = 2*X(1);
    // dHdX2 = 2*X(2);
    // dHdX3 = 2*X(3);
    // dHdX4 = 2*X(4);
    // if kf_structure.flagestimateaccelerometerbias
    //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
    // else
    //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
    // end
    if (filter_struct->FlagEstimateAccelerometerBias)
    {
        GMATRIX_SETSIZE(H,1,13);
        GMATRIX_ZEROES(H);
    }
    else
    {
        GMATRIX_SETSIZE(H,1,10);
        GMATRIX_ZEROES(H);
    }
    GMATRIX_DATA(H, 1, 1) = 2.0 * GMATRIX_DATA(X_predicted, 1, 1);
    GMATRIX_DATA(H, 1, 2) = 2.0 * GMATRIX_DATA(X_predicted, 2, 1);
    GMATRIX_DATA(H, 1, 3) = 2.0 * GMATRIX_DATA(X_predicted, 3, 1);
    GMATRIX_DATA(H, 1, 4) = 2.0 * GMATRIX_DATA(X_predicted, 4, 1);
    // v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
    GMATRIX_SETSIZE(V, 1, 1);
    aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 1, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 2, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 3, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 4, 1));
    GMATRIX_DATA(V,1,1) = 1.0 - aux;


    // K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
    // kf_structure.X = X + K*v;
    // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
    GMATRIX_SETSIZE(R, 1, 1);
    PGMATRIX_COPY(&R, filter_struct->pR_pseudomeasurementnorm);
    kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);
    localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);

    // Retorna
    return 1; 
}

int localization_filter_correction_ekf2(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T)
{
// %%% corrige usando medidas do algoritmo TRIAD: 
// %%% Fun�ao de medi�ao: [eye(4) zeros(4,6)]*X (a medida � dada pelo estimador TRIAD.
// %%% Estimar a medi�ao do quaternion assim como da matriz R usando UT;

    GMATRIX_DECLARE(H, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(V, 4, 1);
    GMATRIX_DECLARE(Ym, 4, 1);
    GMATRIX_DECLARE(Py, 4, 4);
    GMATRIX_DECLARE(R, 4, 4);
    GMATRIX_DECLARE(X_predicted, LOCALIZATION_MAXSTATESIZE, 1);
    GMATRIX_DECLARE(P_predicted, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy1, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy2, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy3, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy4, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    DUMMY_MATRICES MatDummy;
    double aux;

    MatDummy.pMat1 = &MatDummy1;
    MatDummy.pMat2 = &MatDummy2;
    MatDummy.pMat3 = &MatDummy3;
    MatDummy.pMat4 = &MatDummy4;

    // Defini��o dos tamanhos das matrizes deve ser feito aqui, o que evita aloca��o dim�mica
    GMATRIX_SETSIZE(X_predicted, filter_struct->Nstates, 1);
    GMATRIX_SETSIZE(P_predicted, filter_struct->Nstates, filter_struct->Nstates);

    // Medi��o pelo TRIAD
    if(magnetometer_measure_ptr->FlagValidMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);

        // if kf_structure.flagestimateaccelerometerbias
        //     H = [eye(4) zeros(4,9)];
        // else
        //     H = [eye(4) zeros(4,6)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 4, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 4, 10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H, 1, 1) = 1.0;
        GMATRIX_DATA(H, 2, 2) = 1.0;
        GMATRIX_DATA(H, 3, 3) = 1.0;
        GMATRIX_DATA(H, 4, 4) = 1.0;

        // [Ym,Py] = ConvertedMeasurementTRIAD('ekf2',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
        GMATRIX_SETSIZE(Ym, 4, 1);
        GMATRIX_SETSIZE(Py, 4, 4);
        localization_converted_measurement_triad(&Ym, &Py, &X_predicted, &P_predicted, imu_measure_ptr, magnetometer_measure_ptr, pM, pG, filter_struct->FlagEstimateAccelerometerBias);
     // GMATRIX_PRINT_MATLABFORM(Ym);
     // GMATRIX_PRINT_MATLABFORM(Py);

        // v = (Ym - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, &Ym, MatDummy.pMat1);

        // K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad_ekf);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_ADD_COPY(&R, &Py, filter_struct->pR_convertedmeasurementtriad);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do GPS:
    // %%% Fun�ao de medi�ao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
    if(gps_measure_ptr->FlagValidPositionMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
        // else
        //     H = [zeros(3,4) eye(3,3) zeros(3,3)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 3, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 3, 10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H, 1, 5) = 1.0;
        GMATRIX_DATA(H, 2, 6) = 1.0;
        GMATRIX_DATA(H, 3, 7) = 1.0;
        // v = (gpsmeasure.p - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, gps_measure_ptr->pPosition, MatDummy.pMat1);

        // K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_COPY(&R,gps_measure_ptr->pPPosition);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do GPS:
    // %%% Fun�ao de medi�ao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
    if(gps_measure_ptr->FlagValidVelocityMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
        // else
        //     H = [zeros(3,4) zeros(3,3) eye(3,3)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 3, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 3, 10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H, 1, 8) = 1.0;
        GMATRIX_DATA(H, 2, 9) = 1.0;
        GMATRIX_DATA(H, 3, 10) = 1.0;
        // v = (gpsmeasure.v - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, gps_measure_ptr->pVelocity, MatDummy.pMat1);

        // K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_COPY(&R,gps_measure_ptr->pPVelocity);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do sonar:
    // %%% Fun�ao de medi�ao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
    // %%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
    if(sonar_measure_ptr->FlagValidMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // dHdX1 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(1);
        // dHdX2 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(2);
        // dHdX3 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(3);
        // dHdX4 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(4);
        // dHdX7 =  1/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0];
        // else
        //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 1, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 1, 10);
            GMATRIX_ZEROES(H);
        }

        aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 1, 1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 2, 1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 3, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 4, 1));

        if(aux != 0.0)
        { // situa��o em que a medi��o seria diferente de infinito
            GMATRIX_DATA(H, 1, 1) = -GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 1, 1);
            GMATRIX_DATA(H, 1, 2) =  GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 2, 1);
            GMATRIX_DATA(H, 1, 3) =  GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 3, 1);
            GMATRIX_DATA(H, 1, 4) = -GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 4, 1);
            GMATRIX_DATA(H,1,7) =  1.0/(aux);

            // v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
            GMATRIX_SETSIZE(V, 1, 1);
            GMATRIX_DATA(V, 1, 1) = sonar_measure_ptr->range - GMATRIX_DATA(X_predicted, 7, 1) / aux;

            // K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
            // kf_structure.X = X + K*v;
            // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
            GMATRIX_SETSIZE(R, 1, 1);
            GMATRIX_DATA(R, 1, 1) = sonar_measure_ptr->rangevariance;
            kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

            // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
            localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
        }
    }
    // %%% corrige a norma do quaternion usando pseudo-medi�oes:
    // %%% Fun�ao de medi�ao: 1 = q0^2 + q1^2 + q2^2 + q3^2
    // %%% Fun�ao de medi�ao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
    // X = kf_structure.X;
    PGMATRIX_COPY(&X_predicted, filter_struct->pX);

    // P = kf_structure.P;
    PGMATRIX_COPY(&P_predicted, filter_struct->pP);
    // dHdX1 = 2*X(1);
    // dHdX2 = 2*X(2);
    // dHdX3 = 2*X(3);
    // dHdX4 = 2*X(4);
    // if kf_structure.flagestimateaccelerometerbias
    //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
    // else
    //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
    // end
    if (filter_struct->FlagEstimateAccelerometerBias)
    {
        GMATRIX_SETSIZE(H, 1, 13);
        GMATRIX_ZEROES(H);
    }
    else
    {
        GMATRIX_SETSIZE(H, 1, 10);
        GMATRIX_ZEROES(H);
    }
    GMATRIX_DATA(H, 1, 1) = 2.0*GMATRIX_DATA(X_predicted, 1, 1);
    GMATRIX_DATA(H, 1, 2) = 2.0*GMATRIX_DATA(X_predicted, 2, 1);
    GMATRIX_DATA(H, 1, 3) = 2.0*GMATRIX_DATA(X_predicted, 3, 1);
    GMATRIX_DATA(H, 1, 4) = 2.0*GMATRIX_DATA(X_predicted, 4, 1);
    // v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
    GMATRIX_SETSIZE(V, 1, 1);
    aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 1, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 2, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 3, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 4, 1));
    GMATRIX_DATA(V, 1, 1) = 1.0 - aux;


    // K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
    // kf_structure.X = X + K*v;
    // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
    GMATRIX_SETSIZE(R, 1, 1);
    PGMATRIX_COPY(&R, filter_struct->pR_pseudomeasurementnorm);
    kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);
    localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);

    // Retorna
    return 1; 
}

void localization_converted_measurement_triad_UT_nonbiased(PGMATRIX pY, PGMATRIX pX, PGMATRIX pParameters)
{
    QUATERNIONS_DECLARE(q_predicted);
    ImuMeasure imu_measure;
    MagnetometerMeasure magnetometer_measure;
    GMATRIX_DECLARE(M, 3, 1);
    GMATRIX_DECLARE(G, 3, 1);

    QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX, X_q0, 1);
    QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX, X_q1, 1);
    QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX, X_q2, 1);
    QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX, X_q3, 1);

    imu_measure.ax = PGMATRIX_DATA(pX, X_vz + 1, 1);
    imu_measure.ay = PGMATRIX_DATA(pX, X_vz + 2, 1);
    imu_measure.az = PGMATRIX_DATA(pX, X_vz + 3, 1);

    magnetometer_measure.mx = PGMATRIX_DATA(pX, X_vz + 4, 1);
    magnetometer_measure.my = PGMATRIX_DATA(pX, X_vz + 5, 1);
    magnetometer_measure.mz = PGMATRIX_DATA(pX, X_vz + 6, 1);

    GMATRIX_DATA(M, 1, 1) = PGMATRIX_DATA(pParameters, 1, 1);
    GMATRIX_DATA(M, 2, 1) = PGMATRIX_DATA(pParameters, 2, 1);
    GMATRIX_DATA(M, 3, 1) = PGMATRIX_DATA(pParameters, 3, 1);
    GMATRIX_DATA(G, 1, 1) = PGMATRIX_DATA(pParameters, 4, 1);
    GMATRIX_DATA(G, 2, 1) = PGMATRIX_DATA(pParameters, 5, 1);
    GMATRIX_DATA(G, 3, 1) = PGMATRIX_DATA(pParameters, 6, 1);

    localization_triad(pY, &imu_measure, &magnetometer_measure, &M, &G, &q_predicted);
}

void localization_converted_measurement_triad_UT_biased(PGMATRIX pY, PGMATRIX pX, PGMATRIX pParameters)
{
    QUATERNIONS_DECLARE(q_predicted);
    ImuMeasure imu_measure;
    MagnetometerMeasure magnetometer_measure;
    GMATRIX_DECLARE(M, 3, 1);
    GMATRIX_DECLARE(G, 3, 1);

    QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX, X_q0, 1);
    QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX, X_q1, 1);
    QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX, X_q2, 1);
    QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX, X_q3, 1);

    imu_measure.ax = PGMATRIX_DATA(pX, X_baz + 1, 1) - PGMATRIX_DATA(pX, X_bax, 1);
    imu_measure.ay = PGMATRIX_DATA(pX, X_baz + 2, 1) - PGMATRIX_DATA(pX, X_bay, 1);
    imu_measure.az = PGMATRIX_DATA(pX, X_baz + 3, 1) - PGMATRIX_DATA(pX, X_baz, 1);

    magnetometer_measure.mx = PGMATRIX_DATA(pX, X_baz + 4, 1);
    magnetometer_measure.my = PGMATRIX_DATA(pX, X_baz + 5, 1);
    magnetometer_measure.mz = PGMATRIX_DATA(pX, X_baz + 6, 1);

    GMATRIX_DATA(M, 1, 1) = PGMATRIX_DATA(pParameters, 1, 1);
    GMATRIX_DATA(M, 2, 1) = PGMATRIX_DATA(pParameters, 2, 1);
    GMATRIX_DATA(M, 3, 1) = PGMATRIX_DATA(pParameters, 3, 1);
    GMATRIX_DATA(G, 1, 1) = PGMATRIX_DATA(pParameters, 4, 1);
    GMATRIX_DATA(G, 2, 1) = PGMATRIX_DATA(pParameters, 5, 1);
    GMATRIX_DATA(G, 3, 1) = PGMATRIX_DATA(pParameters, 6, 1);

    localization_triad(pY, &imu_measure, &magnetometer_measure, &M, &G, &q_predicted);
}

int localization_converted_measurement_triad(PGMATRIX pYm, PGMATRIX pPy, PGMATRIX pX_predicted, PGMATRIX pP_predicted, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias)
{
    /* Nesse caso, ser� aplicada a unscented transform para o seguinte problema:
            Y = f(Z)
       com Z = [X; u_imu_acc; u_mag] e Pz = diag[Px, Pu_imu_acc, Pu_mag];
    */

    int i;
    UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM_INIT(4, LOCALIZATION_MAXSTATESIZE+6);
    QUATERNIONS_DECLARE(q_predicted);
    GMATRIX_DECLARE(Parameters, 6, 1);

    QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q0, 1);
    QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q1, 1);
    QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q2, 1);
    QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q3, 1);

    // Defini��o dos tamanhos das matrizes deve ser feito aqui,  o que evita aloca��o dim�mica
    GMATRIX_SETSIZE(XMean, pX_predicted->Nr + 6, 1);
    GMATRIX_SETSIZE(XCov, (pX_predicted->Nr + 6), (pX_predicted->Nr + 6));
    GMATRIX_SETSIZE(_utX, (pX_predicted->Nr + 6), 1);
    GMATRIX_SETSIZE(_utXSamples, (pX_predicted->Nr + 6), (2*(pX_predicted->Nr + 6) + 1));
    GMATRIX_SETSIZE(_utYSamples, (4), (2*(pX_predicted->Nr + 6) + 1));
    GMATRIX_SETSIZE(_utWSamples, (1), (2*(pX_predicted->Nr + 6) + 1));
    GMATRIX_SETSIZE(_utQSamples, (pX_predicted->Nr + 6),  (pX_predicted->Nr + 6));

    for(i=1;i<=pX_predicted->Nr; +  + i)
    {
        GMATRIX_DATA(XMean, i, 1) = PGMATRIX_DATA(pX_predicted, i, 1);
    }
    GMATRIX_DATA(XMean, pX_predicted->Nr + 1, 1) = imu_measure_ptr->ax;
    GMATRIX_DATA(XMean, pX_predicted->Nr + 2, 1) = imu_measure_ptr->ay;
    GMATRIX_DATA(XMean, pX_predicted->Nr + 3, 1) = imu_measure_ptr->az;
    GMATRIX_DATA(XMean, pX_predicted->Nr + 4, 1) = magnetometer_measure_ptr->mx;
    GMATRIX_DATA(XMean, pX_predicted->Nr + 5, 1) = magnetometer_measure_ptr->my;
    GMATRIX_DATA(XMean, pX_predicted->Nr + 6, 1) = magnetometer_measure_ptr->mz;
    GMATRIX_ZEROES(XCov);
    PGMATRIX_SUBMATRIX_COPY(&XCov, 1, 1, pP_predicted);
    GMATRIX_DATA(XCov, pX_predicted->Nr + 1, pX_predicted->Nr + 1) = imu_measure_ptr->axvariance;
    GMATRIX_DATA(XCov, pX_predicted->Nr + 2, pX_predicted->Nr + 2) = imu_measure_ptr->ayvariance;
    GMATRIX_DATA(XCov, pX_predicted->Nr + 3, pX_predicted->Nr + 3) = imu_measure_ptr->azvariance;
    GMATRIX_DATA(XCov, pX_predicted->Nr + 4, pX_predicted->Nr + 4) = magnetometer_measure_ptr->mxvariance;
    GMATRIX_DATA(XCov, pX_predicted->Nr + 5, pX_predicted->Nr + 5) = magnetometer_measure_ptr->myvariance;
    GMATRIX_DATA(XCov, pX_predicted->Nr + 6, pX_predicted->Nr + 6) = magnetometer_measure_ptr->mzvariance;

    GMATRIX_DATA(Parameters, 1, 1) = PGMATRIX_DATA(pM, 1, 1);
    GMATRIX_DATA(Parameters, 2, 1) = PGMATRIX_DATA(pM, 2, 1);
    GMATRIX_DATA(Parameters, 3, 1) = PGMATRIX_DATA(pM, 3, 1);
    GMATRIX_DATA(Parameters, 4, 1) = PGMATRIX_DATA(pG, 1, 1);
    GMATRIX_DATA(Parameters, 5, 1) = PGMATRIX_DATA(pG, 2, 1);
    GMATRIX_DATA(Parameters, 6, 1) = PGMATRIX_DATA(pG, 3, 1);

    if (FlagEstimateAccelerometerBias)
    {
        UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM(localization_converted_measurement_triad_UT_biased);
    }
    else
    {
        UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM(localization_converted_measurement_triad_UT_nonbiased);
    }

    PGMATRIX_COPY(pYm,&YMean);
    PGMATRIX_COPY(pPy,&YCov);

    rotation_quaternionsnormalize(pYm);
    rotation_quaternionscorrectsign(pYm, &q_predicted);

    return 1; 
}

int localization_converted_measurement_triad_dg_du_imu(PGMATRIX pdg_du_imu, PGMATRIX pX_predicted, PGMATRIX pP_predicted, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias)
{
    ImuMeasure imu_measure;
    QUATERNIONS_DECLARE(q_predicted);
    QUATERNIONS_DECLARE(y);
    QUATERNIONS_DECLARE(y_ax);
    QUATERNIONS_DECLARE(y_ay);
    QUATERNIONS_DECLARE(y_az);
    double delta;

    // Define q_predicted
    QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q0, 1);
    QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q1, 1);
    QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q2, 1);
    QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX_predicted, X_q3, 1);

    // Considera polariza��o dos aceler�metros
    imu_measure.wx = imu_measure_ptr->wx;
    imu_measure.wy = imu_measure_ptr->wy;
    imu_measure.wz = imu_measure_ptr->wz;
    if (FlagEstimateAccelerometerBias)
    {
        imu_measure.ax = imu_measure_ptr->ax - PGMATRIX_DATA(pX_predicted, X_bax, 1);
        imu_measure.ay = imu_measure_ptr->ay - PGMATRIX_DATA(pX_predicted, X_bay, 1);
        imu_measure.az = imu_measure_ptr->az - PGMATRIX_DATA(pX_predicted, X_baz, 1);
    }
    else
    {
        imu_measure.ax = imu_measure_ptr->ax;
        imu_measure.ay = imu_measure_ptr->ay;
        imu_measure.az = imu_measure_ptr->az;
    }

    // M�todo aproximado de calculo de dg_du_imu para o caso da triad
    delta = 1e-5;

    // y = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    localization_triad(&y, &imu_measure, magnetometer_measure_ptr, pM, pG, &q_predicted);

    // imumeasure.ax = imumeasure.ax + delta;
    // y_ax = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    // imumeasure.ax = imumeasure.ax - delta;
    imu_measure.ax += delta;
    localization_triad(&y_ax, &imu_measure, magnetometer_measure_ptr, pM, pG, &q_predicted);
    imu_measure.ax -= delta;

    // imumeasure.ay = imumeasure.ay + delta;
    // y_ay = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    // imumeasure.ay = imumeasure.ay - delta;
    imu_measure.ay += delta;
    localization_triad(&y_ay, &imu_measure, magnetometer_measure_ptr, pM, pG, &q_predicted);
    imu_measure.ay -= delta;

    // imumeasure.az = imumeasure.az + delta;
    // y_az = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    // imumeasure.az = imumeasure.az - delta;
    imu_measure.az += delta;
    localization_triad(&y_az, &imu_measure, magnetometer_measure_ptr, pM, pG, &q_predicted);
    imu_measure.az -= delta;
    
    // dg_du_imu = [(y_ax-y)/delta (y_ay-y)/delta (y_az-y)/delta zeros(4,3)];
    PGMATRIX_SETSIZE(pdg_du_imu, 4, 6);
    PGMATRIX_ZEROES(pdg_du_imu);
    PGMATRIX_DATA(pdg_du_imu, 1, 1) = ( GMATRIX_DATA(y_ax, 1, 1) - GMATRIX_DATA(y, 1, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 2, 1) = ( GMATRIX_DATA(y_ax, 2, 1) - GMATRIX_DATA(y, 2, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 3, 1) = ( GMATRIX_DATA(y_ax, 3, 1) - GMATRIX_DATA(y, 3, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 4, 1) = ( GMATRIX_DATA(y_ax, 4, 1) - GMATRIX_DATA(y, 4, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 1, 2) = ( GMATRIX_DATA(y_ay, 1, 1) - GMATRIX_DATA(y, 1, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 2, 2) = ( GMATRIX_DATA(y_ay, 2, 1) - GMATRIX_DATA(y, 2, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 3, 2) = ( GMATRIX_DATA(y_ay, 3, 1) - GMATRIX_DATA(y, 3, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 4, 2) = ( GMATRIX_DATA(y_ay, 4, 1) - GMATRIX_DATA(y, 4, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 1, 3) = ( GMATRIX_DATA(y_az, 1, 1) - GMATRIX_DATA(y, 1, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 2, 3) = ( GMATRIX_DATA(y_az, 2, 1) - GMATRIX_DATA(y, 2, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 3, 3) = ( GMATRIX_DATA(y_az, 3, 1) - GMATRIX_DATA(y, 3, 1) ) / delta;
    PGMATRIX_DATA(pdg_du_imu, 4, 3) = ( GMATRIX_DATA(y_az, 4, 1) - GMATRIX_DATA(y, 4, 1) ) / delta;

    return 1;
}

int localization_filter_prediction(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T)
{
    // Predi��o conforme a arquitetura
    switch(filter_struct->AlgorithmCode)
    {
        case LOCALIZATION_ALGORITHMCODE_EKF2:
            if(!localization_filter_prediction_ekf2(filter_struct, imu_measure_ptr, pG, T)) return 0;
            break;
        case LOCALIZATION_ALGORITHMCODE_EKF_DECOUPLED:
            if(!localization_filter_prediction_ekf_decoupled(filter_struct, imu_measure_ptr, pG, T)) return 0;
            break;
        case LOCALIZATION_ALGORITHMCODE_CEKF:
            if(!localization_filter_prediction_cekf(filter_struct, imu_measure_ptr, pG, T)) return 0;
            break;
        default:
            return 0;
    }

    // Reset em caso de matriz P negativa:
/*  if(sum(eig(kf_structure.P)<0)>0)
        % filter reset:
        disp('*** Warning: Kalman filter reset');
        kf_structure.P = kf_structure.Preset;
    end
*/
    // Retorna
    return 1; 
}                      

int localization_filter_prediction_cekf(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T)
{
    GMATRIX_DECLARE(X_previous, LOCALIZATION_MAXSTATESIZE, 1);
    GMATRIX_DECLARE(P_previous, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(u_imu, 6, 1);
    GMATRIX_DECLARE(P_u_imu, 6, 6);
    GMATRIX_DECLARE(df_dx, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(df_du, LOCALIZATION_MAXSTATESIZE, 6);
    GMATRIX_DECLARE(K_previous, LOCALIZATION_MAXSTATESIZE, 4);
    GMATRIX_DECLARE(H, 4, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(Qtilde, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy1, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy2, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);

    /*
        filter_struct->pS_previous_cekf = PGMATRIX_ALLOC(filter_struct->Nstates,4);
        PGMATRIX_ZEROES(filter_struct->pS_previous_cekf);
        filter_struct->pR_previous_cekf = PGMATRIX_ALLOC(4,4);
        PGMATRIX_IDENTITY(filter_struct->pR_previous_cekf);
        filter_struct->pA_imu_P_imu_cekf = PGMATRIX_ALLOC(filter_struct->Nstates,6);
        PGMATRIX_ZEROES(filter_struct->pA_imu_P_imu_cekf);
        filter_struct->pinnovation_previous_cekf = PGMATRIX_ALLOC(4,1);
        PGMATRIX_ZEROES(filter_struct->pinnovation_previous_cekf);
    */

    // Defini��o dos tamanhos das matrizes deve ser feito aqui, o que evita aloca��o dim�mica
    GMATRIX_SETSIZE(X_previous, filter_struct->Nstates, 1);
    GMATRIX_SETSIZE(P_previous, filter_struct->Nstates, filter_struct->Nstates);
    GMATRIX_SETSIZE(df_dx, filter_struct->Nstates, filter_struct->Nstates);
    GMATRIX_SETSIZE(df_du, filter_struct->Nstates, 6);
    GMATRIX_SETSIZE(K_previous, filter_struct->Nstates, 4);
    GMATRIX_SETSIZE(Qtilde, filter_struct->Nstates, filter_struct->Nstates);

    // X_previous = kf_structure.X;
    PGMATRIX_COPY(&X_previous, filter_struct->pX);

    // P_previous = kf_structure.P;
    PGMATRIX_COPY(&P_previous, filter_struct->pP);

    // u_imu = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
    GMATRIX_DATA(u_imu, 1, 1) = imu_measure_ptr->ax;
    GMATRIX_DATA(u_imu, 2, 1) = imu_measure_ptr->ay;
    GMATRIX_DATA(u_imu, 3, 1) = imu_measure_ptr->az;
    GMATRIX_DATA(u_imu, 4, 1) = imu_measure_ptr->wx;
    GMATRIX_DATA(u_imu, 5, 1) = imu_measure_ptr->wy;
    GMATRIX_DATA(u_imu, 6, 1) = imu_measure_ptr->wz;

    // Pu_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
    GMATRIX_ZEROES(P_u_imu);
    GMATRIX_DATA(P_u_imu, 1, 1) = imu_measure_ptr->axvariance;
    GMATRIX_DATA(P_u_imu, 2, 2) = imu_measure_ptr->ayvariance;
    GMATRIX_DATA(P_u_imu, 3, 3) = imu_measure_ptr->azvariance;
    GMATRIX_DATA(P_u_imu, 4, 4) = imu_measure_ptr->wxvariance;
    GMATRIX_DATA(P_u_imu, 5, 5) = imu_measure_ptr->wyvariance;
    GMATRIX_DATA(P_u_imu, 6, 6) = imu_measure_ptr->wzvariance;

    // K_previous = kf_structure.S_previous_cekf * inv(kf_structure.R_previous_cekf);
    PGMATRIX_INVERSE_COPY(&MatDummy1, filter_struct->pR_previous_cekf);
    PGMATRIX_MULTIPLY_COPY(&K_previous, filter_struct->pS_previous_cekf, &MatDummy1);


    // kf_structure.X = localization_filter_model_f('evaluate',x_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias) + K_previous*kf_structure.innovation_previous_cekf;
    localization_filter_process_model_evaluate(filter_struct->pX, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);
    PGMATRIX_MULTIPLY_ADD(filter_struct->pX, &K_previous, filter_struct->pinnovation_previous_cekf);

    // df_dx = localization_filter_model_f('df_dx',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_df_dx(&df_dx, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // df_du = localization_filter_model_f('df_du',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_df_du(&df_du, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // Qtilde = kf_structure.Q_cekf + df_du_imu*P_u_imu*df_du_imu';
    PGMATRIX_TRIPLEMULTIPLY_COPY_EXTENDED(&Qtilde, &df_du, 0, &P_u_imu, 0, &df_du, 1, &MatDummy1);
    PGMATRIX_ADD(&Qtilde, filter_struct->pQ);

    // if kf_structure.flagestimateaccelerometerbias
    //     H = [eye(4) zeros(4,9)];
    // else
    //     H = [eye(4) zeros(4,6)];
    // end
    if (filter_struct->FlagEstimateAccelerometerBias)
    {
        GMATRIX_SETSIZE(H, 4, 13);
        GMATRIX_ZEROES(H);
    }
    else
    {
        GMATRIX_SETSIZE(H, 4, 10);
        GMATRIX_ZEROES(H);
    }
    GMATRIX_DATA(H, 1, 1) = 1.0;
    GMATRIX_DATA(H, 2, 2) = 1.0;
    GMATRIX_DATA(H, 3, 3) = 1.0;
    GMATRIX_DATA(H, 4, 4) = 1.0;

    // kf_structure.P = (df_dx-K_previous*H)*P_previous*(df_dx-K_previous*H)' + Qtilde - K_previous*kf_structure.R_previous_cekf*K_previous';
    PGMATRIX_TRIPLEMULTIPLY_COPY_EXTENDED(filter_struct->pP, &K_previous, 0, filter_struct->pR_previous_cekf, 0, &K_previous, 1, &MatDummy1);
    PGMATRIX_MULTIPLY_CONST(filter_struct->pP, -1.0); // kf_structure.P = - K_previous*kf_structure.R_previous_cekf*K_previous';
    PGMATRIX_ADD(filter_struct->pP, &Qtilde); // kf_structure.P = Qtilde - K_previous*kf_structure.R_previous_cekf*K_previous';
    PGMATRIX_MULTIPLY_COPY(&MatDummy1, &K_previous, &H);
    PGMATRIX_SUBTRACT_COPY(&MatDummy2, &df_dx, &MatDummy1); // MatDummy2 = df_dx-K_previous*H
    PGMATRIX_TRIPLEMULTIPLY_ADD_EXTENDED(filter_struct->pP, &MatDummy2, 0, &P_previous, 0, &MatDummy2, 1, &MatDummy1);

    // kf_structure.A_imu_P_imu_cekf = df_du_imu*P_u_imu;
    PGMATRIX_MULTIPLY_COPY(filter_struct->pA_imu_P_imu_cekf, &df_du, &P_u_imu);

    // Retorna
    return 1; 
}  

int localization_filter_prediction_ekf2(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T)
{
    GMATRIX_DECLARE(X_previous, LOCALIZATION_MAXSTATESIZE, 1);
    GMATRIX_DECLARE(P_previous, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(u_imu, 6, 1);
    GMATRIX_DECLARE(P_u_imu, 6, 6);
    GMATRIX_DECLARE(df_dx, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(df_du, LOCALIZATION_MAXSTATESIZE, 6);
    GMATRIX_DECLARE(MatDummy, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);

    // Defini��o dos tamanhos das matrizes deve ser feito aqui, o que evita aloca��o dim�mica
    GMATRIX_SETSIZE(X_previous,filter_struct->Nstates, 1);
    GMATRIX_SETSIZE(P_previous, filter_struct->Nstates, filter_struct->Nstates);
    GMATRIX_SETSIZE(df_dx, filter_struct->Nstates, filter_struct->Nstates);
    GMATRIX_SETSIZE(df_du, filter_struct->Nstates, 6);

    // X_previous = kf_structure.X;
    PGMATRIX_COPY(&X_previous, filter_struct->pX);

    // P_previous = kf_structure.P;
    PGMATRIX_COPY(&P_previous, filter_struct->pP);

    // u_imu = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
    GMATRIX_DATA(u_imu, 1, 1) = imu_measure_ptr->ax;
    GMATRIX_DATA(u_imu, 2, 1) = imu_measure_ptr->ay;
    GMATRIX_DATA(u_imu, 3, 1) = imu_measure_ptr->az;
    GMATRIX_DATA(u_imu, 4, 1) = imu_measure_ptr->wx;
    GMATRIX_DATA(u_imu, 5, 1) = imu_measure_ptr->wy;
    GMATRIX_DATA(u_imu, 6, 1) = imu_measure_ptr->wz;

    // Pu_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
    GMATRIX_ZEROES(P_u_imu);
    GMATRIX_DATA(P_u_imu, 1, 1) = imu_measure_ptr->axvariance;
    GMATRIX_DATA(P_u_imu, 2, 2) = imu_measure_ptr->ayvariance;
    GMATRIX_DATA(P_u_imu, 3, 3) = imu_measure_ptr->azvariance;
    GMATRIX_DATA(P_u_imu, 4, 4) = imu_measure_ptr->wxvariance;
    GMATRIX_DATA(P_u_imu, 5, 5) = imu_measure_ptr->wyvariance;
    GMATRIX_DATA(P_u_imu, 6, 6) = imu_measure_ptr->wzvariance;

    // kf_structure.X = localization_filter_model_f('evaluate',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_evaluate(filter_struct->pX, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // df_dx = localization_filter_model_f('df_dx',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_df_dx(&df_dx, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // df_du = localization_filter_model_f('df_du',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_df_du(&df_du, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // kf_structure.P = df_dx*P_previous*df_dx' + df_du*P_u_imu*df_du' + kf_structure.Q_ekf;
    PGMATRIX_TRIPLEMULTIPLY_COPY_EXTENDED(filter_struct->pP, &df_dx, 0, &P_previous, 0, &df_dx, 1, &MatDummy);
    PGMATRIX_TRIPLEMULTIPLY_ADD_EXTENDED(filter_struct->pP, &df_du, 0, &P_u_imu, 0, &df_du, 1, &MatDummy);
    PGMATRIX_ADD(filter_struct->pP, filter_struct->pQ);

    // Retorna
    return 1; 
}                      

int localization_filter_prediction_ekf_decoupled(LocalizationFilter *filter_struct, ImuMeasure *imu_measure_ptr, PGMATRIX pG, double T)
{
    GMATRIX_DECLARE(X_previous, LOCALIZATION_MAXSTATESIZE, 1);
    GMATRIX_DECLARE(P_previous, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(u_imu, 6, 1);
    GMATRIX_DECLARE(P_u_imu, 6, 6);
    GMATRIX_DECLARE(df_dx, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(df_du, LOCALIZATION_MAXSTATESIZE, 6);
    GMATRIX_DECLARE(MatDummy, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);

    GMATRIX_SETSIZE(X_previous, filter_struct->Nstates, 1);
    GMATRIX_SETSIZE(P_previous, filter_struct->Nstates, filter_struct->Nstates);
    GMATRIX_SETSIZE(df_dx, filter_struct->Nstates, filter_struct->Nstates);
    GMATRIX_SETSIZE(df_du, filter_struct->Nstates, 6);

    PGMATRIX_COPY(&X_previous, filter_struct->pX);

    PGMATRIX_COPY(&P_previous, filter_struct->pP);

    // u_imu = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
    GMATRIX_DATA(u_imu,1,1) = imu_measure_ptr->ax;
    GMATRIX_DATA(u_imu,2,1) = imu_measure_ptr->ay;
    GMATRIX_DATA(u_imu,3,1) = imu_measure_ptr->az;
    GMATRIX_DATA(u_imu,4,1) = imu_measure_ptr->wx;
    GMATRIX_DATA(u_imu,5,1) = imu_measure_ptr->wy;
    GMATRIX_DATA(u_imu,6,1) = imu_measure_ptr->wz;

    // Pu_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
    GMATRIX_ZEROES(P_u_imu);
    GMATRIX_DATA(P_u_imu,1,1) = imu_measure_ptr->axvariance;
    GMATRIX_DATA(P_u_imu,2,2) = imu_measure_ptr->ayvariance;
    GMATRIX_DATA(P_u_imu,3,3) = imu_measure_ptr->azvariance;
    GMATRIX_DATA(P_u_imu,4,4) = imu_measure_ptr->wxvariance;
    GMATRIX_DATA(P_u_imu,5,5) = imu_measure_ptr->wyvariance;
    GMATRIX_DATA(P_u_imu,6,6) = imu_measure_ptr->wzvariance;

    // kf_structure.X = localization_filter_model_f('evaluate',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_evaluate(filter_struct->pX, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // df_dx = localization_filter_model_f('df_dx',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_df_dx(&df_dx, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // df_du = localization_filter_model_f('df_du',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
    localization_filter_process_model_df_du(&df_du, &X_previous, &u_imu, pG, T, filter_struct->FlagEstimateAccelerometerBias);

    // kf_structure.P = df_dx*P_previous*df_dx' + df_du*P_u_imu*df_du' + kf_structure.Q_ekf;
    PGMATRIX_TRIPLEMULTIPLY_COPY_EXTENDED(filter_struct->pP, &df_dx, 0, &P_previous, 0, &df_dx, 1, &MatDummy);
    PGMATRIX_TRIPLEMULTIPLY_ADD_EXTENDED(filter_struct->pP, &df_du, 0, &P_u_imu, 0, &df_du, 1, &MatDummy);
    PGMATRIX_ADD(filter_struct->pP, filter_struct->pQ);

    // Retorna
    return 1;
}

int localization_filter_correction_ekf_decoupled(LocalizationFilter *filter_struct, GpsMeasure *gps_measure_ptr, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, SonarMeasure *sonar_measure_ptr, PGMATRIX pM, PGMATRIX pG, double T)
{
    GMATRIX_DECLARE(X_predicted, LOCALIZATION_MAXSTATESIZE, 1);
    GMATRIX_DECLARE(P_predicted, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(H, 3, LOCALIZATION_MAXSTATESIZE);
    DCM_DECLARE(R);
    QUATERNIONS_DECLARE(q);

    GMATRIX_DECLARE(tempMatrix, 3, 3);
    GMATRIX_DECLARE(tempMatrix2, 3, 1);
    GMATRIX_DECLARE(tempMatrix3, 3, 3);
    GMATRIX_DECLARE(V, 3, 3);

    GMATRIX_DECLARE(magnetometer_measurement, 3, 1);
    GMATRIX_DECLARE(accelerometer_measurement, 3, 1);

    GMATRIX_SETSIZE(X_predicted, filter_struct->Nstates, 1);
    GMATRIX_SETSIZE(P_predicted, filter_struct->Nstates, filter_struct->Nstates);

    GMATRIX_DECLARE(MatDummy1, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy2, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy3, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    GMATRIX_DECLARE(MatDummy4, LOCALIZATION_MAXSTATESIZE, LOCALIZATION_MAXSTATESIZE);
    DUMMY_MATRICES MatDummy;

    MatDummy.pMat1 = &MatDummy1;
    MatDummy.pMat2 = &MatDummy2;
    MatDummy.pMat3 = &MatDummy3;
    MatDummy.pMat4 = &MatDummy4;

    double M_norm = 0.0;
    double M_1 = 0.0;
    double M_2 = 0.0;
    double M_3 = 0.0;
    double G_norm = 0.0;
    double G_1 = 0.0;
    double G_2 = 0.0;
    double G_3 = 0.0;
    double q0 = 0.0;
    double q1 = 0.0;
    double q2 = 0.0;
    double q3 = 0.0;

    double accelerometer_R = 10000.0;

    static double accelerometer_filtered = 0.0;

    double aux = 0.0;

    // First run: set new alpha
    if(accelerometer_filtered == 0.0)
    {
        // remove the pow function
        accelerometer_filtered = sqrt(pow(imu_measure_ptr->ax, 2) + pow(imu_measure_ptr->ay, 2) + pow(imu_measure_ptr->az, 2));
    }

    //        % Magnetometer Correction - Attitude
    //        % H is the measurement function for the magnetometer - in this
    //        % case, we have:
    //        % b = A(q)r
    //        % where b is the actual measurement (data from the magnetometer),
    //        % A(q) is the DCM written in terms of the quaternions, and r is the
    //        % reference frame (the local magnetic field in NED coordinates).
    //        % Therefore: h(X) = A(q)r
    //        % Here, A(q) = [ q0^2 + q1^2 - q2^2 - q3^2, 2*q1*q2 - 2*q0*q3,                 2*q0*q2 + 2*q1*q3]
    //        %              [ 2*q0*q3 + 2*q1*q2,         q0^2 - q1^2 + q2^2 - q3^2,         2*q2*q3 - 2*q0*q1]
    //        %              [ 2*q1*q3 - 2*q0*q2,         2*q0*q1 + 2*q2*q3,                 q0^2 - q1^2 - q2^2 + q3^2]
    //        % h(X) is the actual magnetometer measurement (normalized) and r is
    //        % local magnetic field (normalized as well).
    //        % Reference: model: Unscented Filtering for Spacecraft Attitude
    //        % Estimation, Crassidis, Markley
    //        % Rotation matrix: Strapdown Intertial navigation Technology, page 48

    if(magnetometer_measure_ptr->FlagValidMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);

        M_norm = sqrt(pow(PGMATRIX_DATA(pM, 1, 1), 2) + pow(PGMATRIX_DATA(pM, 2, 1), 2) + pow(PGMATRIX_DATA(pM, 3, 1), 2));
        M_1 = PGMATRIX_DATA(pM, 1, 1);
        M_2 = PGMATRIX_DATA(pM, 2, 1);
        M_3 = PGMATRIX_DATA(pM, 3, 1);
        q0 = GMATRIX_DATA(X_predicted, X_q0, 1);
        q1 = GMATRIX_DATA(X_predicted, X_q1, 1);
        q2 = GMATRIX_DATA(X_predicted, X_q2, 1);
        q3 = GMATRIX_DATA(X_predicted, X_q3, 1);

        // % Calculating the inovation term directly:
        // R = quaternions2dcm(x(1:4));
        QUATERNIONS_Q0(q) = GMATRIX_DATA(X_predicted, X_q0, 1);
        QUATERNIONS_Q1(q) = GMATRIX_DATA(X_predicted, X_q1, 1);
        QUATERNIONS_Q2(q) = GMATRIX_DATA(X_predicted, X_q2, 1);
        QUATERNIONS_Q3(q) = GMATRIX_DATA(X_predicted, X_q3, 1);
        rotation_quaternions2dcm(&q, &R);

        // tempMatrix = transpose(R)
        GMATRIX_TRANSPOSE_COPY(tempMatrix, R);
        // tempMatrix2 = M/norm(M)
        PGMATRIX_MULTIPLY_CONST_COPY(&tempMatrix2, pM, 1/M_norm);

        // tempMatrix3 = transpose(quaternions2dcm(x(1:4, k)))*(M/norm(M));
        GMATRIX_MULTIPLY_COPY(tempMatrix3, tempMatrix, tempMatrix2);

        // magnetometer_measurement = [mx(k); my(k); mz(k);]/sqrt(mx(k).^2+my(k).^2+mz(k).^2)
        GMATRIX_DATA(magnetometer_measurement, 1, 1) = magnetometer_measure_ptr->mx / sqrt(pow(magnetometer_measure_ptr->mx, 2) + pow(magnetometer_measure_ptr->my, 2) + pow(magnetometer_measure_ptr->mz, 2));
        GMATRIX_DATA(magnetometer_measurement, 2, 1) = magnetometer_measure_ptr->my / sqrt(pow(magnetometer_measure_ptr->mx, 2) + pow(magnetometer_measure_ptr->my, 2) + pow(magnetometer_measure_ptr->mz, 2));
        GMATRIX_DATA(magnetometer_measurement, 3, 1) = magnetometer_measure_ptr->mz / sqrt(pow(magnetometer_measure_ptr->mx, 2) + pow(magnetometer_measure_ptr->my, 2) + pow(magnetometer_measure_ptr->mz, 2));

        // % Calculating the inovation term directly:
        // v = [mx(k); my(k); mz(k);]/sqrt(mx(k).^2+my(k).^2+mz(k).^2) - transpose(quaternions2dcm(x(1:4, k)))*(M/norm(M));
        GMATRIX_SUBTRACT_COPY(V, magnetometer_measurement, tempMatrix3);

        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H ,3 ,13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H ,3 ,10);
            GMATRIX_ZEROES(H);
        }

        //        % First column: derivative in relation to q1
        //        % 2*m1*q0 + 2*m2*q3 - 2*m3*q2
        //        % 2*m2*q0 - 2*m1*q3 + 2*m3*q1
        //        % 2*m1*q2 - 2*m2*q1 + 2*m3*q0
        //        dh_dX(1, 1) = 2*M(1)/norm(M)*q0 + 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q2;
        //        dh_dX(2, 1) = 2*M(2)/norm(M)*q0 - 2*M(1)/norm(M)*q3 + 2*M(3)/norm(M)*q1;
        //        dh_dX(3, 1) = 2*M(1)/norm(M)*q2 - 2*M(2)/norm(M)*q1 + 2*M(3)/norm(M)*q0;
        GMATRIX_DATA(H, 1, 1) = (2.0 * M_1 / M_norm) * q0 + (2.0 * M_2 / M_norm) * q3 - (2.0 * M_3 / M_norm) * q2;
        GMATRIX_DATA(H, 2, 1) = (2.0 * M_2 / M_norm) * q0 - (2.0 * M_1 / M_norm) * q3 + (2.0 * M_3 / M_norm) * q1;
        GMATRIX_DATA(H, 3, 1) = (2.0 * M_1 / M_norm) * q2 - (2.0 * M_2 / M_norm) * q1 + (2.0 * M_3 / M_norm) * q0;

        //        % Second column: derivative in relation to q2
        //        % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
        //        % 2*m1*q2 - 2*m2*q1 + 2*m3*q0
        //        % 2*m1*q3 - 2*m2*q0 - 2*m3*q1
        //        dh_dX(1, 2) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
        //        dh_dX(2, 2) = 2*M(1)/norm(M)*q2 - 2*M(2)/norm(M)*q1 + 2*M(3)/norm(M)*q0;
        //        dh_dX(3, 2) = 2*M(1)/norm(M)*q3 - 2*M(2)/norm(M)*q0 - 2*M(3)/norm(M)*q1;
        GMATRIX_DATA(H, 1, 2) = (2.0 * M_1 / M_norm) * q1 + (2.0 * M_2 / M_norm) * q2 + (2.0 * M_3 / M_norm) * q3;
        GMATRIX_DATA(H, 2, 2) = (2.0 * M_1 / M_norm) * q2 - (2.0 * M_2 / M_norm) * q1 + (2.0 * M_3 / M_norm) * q0;
        GMATRIX_DATA(H, 3, 2) = (2.0 * M_1 / M_norm) * q3 - (2.0 * M_2 / M_norm) * q0 - (2.0 * M_3 / M_norm) * q1;

        //        % Third column: derivative in relation to q3
        //        % 2*m2*q1 - 2*m1*q2 - 2*m3*q0
        //        % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
        //        % 2*m1*q0 + 2*m2*q3 - 2*m3*q2
        //        dh_dX(1, 3) = 2*M(2)/norm(M)*q1 - 2*M(1)/norm(M)*q2 - 2*M(3)/norm(M)*q0;
        //        dh_dX(2, 3) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
        //        dh_dX(3, 3) = 2*M(1)/norm(M)*q0 + 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q2;
        GMATRIX_DATA(H, 1, 3) = (2.0 * M_2 / M_norm) * q1 - (2.0 * M_1 / M_norm) * q2 - (2.0 * M_3 / M_norm) * q0;
        GMATRIX_DATA(H, 2, 3) = (2.0 * M_1 / M_norm) * q1 + (2.0 * M_2 / M_norm) * q2 + (2.0 * M_3 / M_norm) * q3;
        GMATRIX_DATA(H, 3, 3) = (2.0 * M_1 / M_norm) * q0 + (2.0 * M_2 / M_norm) * q3 - (2.0 * M_3 / M_norm) * q2;

        //        % Fourth column: derivative in relation to q4
        //        % 2*m2*q0 - 2*m1*q3 + 2*m3*q1
        //        % 2*m3*q2 - 2*m2*q3 - 2*m1*q0
        //        % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
        //        dh_dX(1, 4) = 2*M(2)/norm(M)*q0 - 2*M(1)/norm(M)*q3 + 2*M(2)/norm(M)*q1;
        //        dh_dX(2, 4) = 2*M(3)/norm(M)*q2 - 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q0;
        //        dh_dX(3, 4) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
        GMATRIX_DATA(H, 1, 4) = (2.0 * M_2 / M_norm) * q0 - (2.0 * M_1 / M_norm) * q3 + (2.0 * M_2 / M_norm) * q1;
        GMATRIX_DATA(H, 2, 4) = (2.0 * M_3 / M_norm) * q2 - (2.0 * M_2 / M_norm) * q3 - (2.0 * M_3 / M_norm) * q0;
        GMATRIX_DATA(H, 3, 4) = (2.0 * M_1 / M_norm) * q1 + (2.0 * M_2 / M_norm) * q2 + (2.0 * M_3 / M_norm) * q3;

       // % Calculate the innovation covariance S
       // S = H*P*transpose(H)+ magnetometer_R;
       // % Kalman gain
       // K = P*transpose(H)*inv(S);
       // % Save previous estimate so we can correct the quaternion sign
       // x_previous = x(:, k);
       // % Update the state estimate
       // x(:, k) = x_previous + K*v;
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, filter_struct->pR_magnetometer, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

       //     % Accelerometer Correction - Attitude
       // if accelerometer_validmeasure(k) == 1
       //     % H is the measurement function for the accelerometer - in this
       //     % case, we have the same as the magnetometer. See above for
       //     % details.
    if(1)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);

        G_norm = sqrt(pow(PGMATRIX_DATA(pG, 1, 1), 2) + pow(PGMATRIX_DATA(pG, 2, 1), 2) + pow(PGMATRIX_DATA(pG, 3, 1), 2));
        G_1 = -1.0 * PGMATRIX_DATA(pG, 1, 1);
        G_2 = -1.0 * PGMATRIX_DATA(pG, 2, 1);
        G_3 = -1.0 * PGMATRIX_DATA(pG, 3, 1);
        q0 = GMATRIX_DATA(X_predicted, X_q0, 1);
        q1 = GMATRIX_DATA(X_predicted, X_q1, 1);
        q2 = GMATRIX_DATA(X_predicted, X_q2, 1);
        q3 = GMATRIX_DATA(X_predicted, X_q3, 1);

        // % Calculating the inovation term directly:
        // R = quaternions2dcm(x(1:4));
        QUATERNIONS_Q0(q) = GMATRIX_DATA(X_predicted, X_q0, 1);
        QUATERNIONS_Q1(q) = GMATRIX_DATA(X_predicted, X_q1, 1);
        QUATERNIONS_Q2(q) = GMATRIX_DATA(X_predicted, X_q2, 1);
        QUATERNIONS_Q3(q) = GMATRIX_DATA(X_predicted, X_q3, 1);
        rotation_quaternions2dcm(&q, &R);

        // tempMatrix = transpose(R)
        GMATRIX_TRANSPOSE_COPY(tempMatrix, R);
        // tempMatrix2 = M/norm(M)
        PGMATRIX_MULTIPLY_CONST_COPY(&tempMatrix2, pG, -1.0/G_norm);

        // tempMatrix3 = transpose(quaternions2dcm(x(1:4, k)))*(-1.0*G/norm(G));
        GMATRIX_MULTIPLY_COPY(tempMatrix3, tempMatrix, tempMatrix2);

        // accelerometer_measurement = [mx(k); my(k); mz(k);]/sqrt(mx(k).^2+my(k).^2+mz(k).^2)
        GMATRIX_DATA(accelerometer_measurement, 1, 1) = imu_measure_ptr->ax / sqrt(pow(imu_measure_ptr->ax, 2.0) + pow(imu_measure_ptr->ay, 2.0) + pow(imu_measure_ptr->az, 2.0));
        GMATRIX_DATA(accelerometer_measurement, 2, 1) = imu_measure_ptr->ay / sqrt(pow(imu_measure_ptr->ax, 2.0) + pow(imu_measure_ptr->ay, 2.0) + pow(imu_measure_ptr->az, 2.0));
        GMATRIX_DATA(accelerometer_measurement, 3, 1) = imu_measure_ptr->az / sqrt(pow(imu_measure_ptr->ax, 2.0) + pow(imu_measure_ptr->ay, 2.0) + pow(imu_measure_ptr->az, 2.0));

        // % Calculating the inovation term directly:
        // v = [ax(k); ay(k); az(k);]/sqrt(ax(k).^2+ay(k).^2+az(k).^2) - transpose(quaternions2dcm(x(1:4, k)))*(-1.0*G/norm(G));
        GMATRIX_SUBTRACT_COPY(V, accelerometer_measurement, tempMatrix3);

        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 3, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 3, 10);
            GMATRIX_ZEROES(H);
        }

        //        % First column: derivative in relation to q1
        //        % 2*m1*q0 + 2*m2*q3 - 2*m3*q2
        //        % 2*m2*q0 - 2*m1*q3 + 2*m3*q1
        //        % 2*m1*q2 - 2*m2*q1 + 2*m3*q0
        //        dh_dX(1, 1) = 2*M(1)/norm(M)*q0 + 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q2;
        //        dh_dX(2, 1) = 2*M(2)/norm(M)*q0 - 2*M(1)/norm(M)*q3 + 2*M(3)/norm(M)*q1;
        //        dh_dX(3, 1) = 2*M(1)/norm(M)*q2 - 2*M(2)/norm(M)*q1 + 2*M(3)/norm(M)*q0;
        GMATRIX_DATA(H, 1, 1) = (2.0 * G_1 / G_norm) * q0 + (2.0 * G_2 / G_norm) * q3 - (2.0 * G_3 / G_norm) * q2;
        GMATRIX_DATA(H, 2, 1) = (2.0 * G_2 / G_norm) * q0 - (2.0 * G_1 / G_norm) * q3 + (2.0 * G_3 / G_norm) * q1;
        GMATRIX_DATA(H, 3, 1) = (2.0 * G_1 / G_norm) * q2 - (2.0 * G_2 / G_norm) * q1 + (2.0 * G_3 / G_norm) * q0;

        //        % Second column: derivative in relation to q2
        //        % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
        //        % 2*m1*q2 - 2*m2*q1 + 2*m3*q0
        //        % 2*m1*q3 - 2*m2*q0 - 2*m3*q1
        //        dh_dX(1, 2) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
        //        dh_dX(2, 2) = 2*M(1)/norm(M)*q2 - 2*M(2)/norm(M)*q1 + 2*M(3)/norm(M)*q0;
        //        dh_dX(3, 2) = 2*M(1)/norm(M)*q3 - 2*M(2)/norm(M)*q0 - 2*M(3)/norm(M)*q1;
        GMATRIX_DATA(H, 1, 2) = (2.0 * G_1 / G_norm) * q1 + (2.0 * G_2 / G_norm) * q2 + (2.0 * G_3 / G_norm) * q3;
        GMATRIX_DATA(H, 2, 2) = (2.0 * G_1 / G_norm) * q2 - (2.0 * G_2 / G_norm) * q1 + (2.0 * G_3 / G_norm) * q0;
        GMATRIX_DATA(H, 3, 2) = (2.0 * G_1 / G_norm) * q3 - (2.0 * G_2 / G_norm) * q0 - (2.0 * G_3 / G_norm) * q1;

        //        % Third column: derivative in relation to q3
        //        % 2*m2*q1 - 2*m1*q2 - 2*m3*q0
        //        % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
        //        % 2*m1*q0 + 2*m2*q3 - 2*m3*q2
        //        dh_dX(1, 3) = 2*M(2)/norm(M)*q1 - 2*M(1)/norm(M)*q2 - 2*M(3)/norm(M)*q0;
        //        dh_dX(2, 3) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
        //        dh_dX(3, 3) = 2*M(1)/norm(M)*q0 + 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q2;
        GMATRIX_DATA(H, 1, 3) = (2.0 * G_2 / G_norm) * q1 - (2.0 * G_1 / G_norm) * q2 - (2.0 * G_3 / G_norm) * q0;
        GMATRIX_DATA(H, 2, 3) = (2.0 * G_1 / G_norm) * q1 + (2.0 * G_2 / G_norm) * q2 + (2.0 * G_3 / G_norm) * q3;
        GMATRIX_DATA(H, 3, 3) = (2.0 * G_1 / G_norm) * q0 + (2.0 * G_2 / G_norm) * q3 - (2.0 * G_3 / G_norm) * q2;

        //        % Fourth column: derivative in relation to q4
        //        % 2*m2*q0 - 2*m1*q3 + 2*m3*q1
        //        % 2*m3*q2 - 2*m2*q3 - 2*m1*q0
        //        % 2*m1*q1 + 2*m2*q2 + 2*m3*q3
        //        dh_dX(1, 4) = 2*M(2)/norm(M)*q0 - 2*M(1)/norm(M)*q3 + 2*M(2)/norm(M)*q1;
        //        dh_dX(2, 4) = 2*M(3)/norm(M)*q2 - 2*M(2)/norm(M)*q3 - 2*M(3)/norm(M)*q0;
        //        dh_dX(3, 4) = 2*M(1)/norm(M)*q1 + 2*M(2)/norm(M)*q2 + 2*M(3)/norm(M)*q3;
        GMATRIX_DATA(H, 1, 4) = (2.0 * G_2 / G_norm) * q0 - (2.0 * G_1 / G_norm) * q3 + (2.0 * G_2 / G_norm) * q1;
        GMATRIX_DATA(H, 2, 4) = (2.0 * G_3 / G_norm) * q2 - (2.0 * G_2 / G_norm) * q3 - (2.0 * G_3 / G_norm) * q0;
        GMATRIX_DATA(H, 3, 4) = (2.0 * G_1 / G_norm) * q1 + (2.0 * G_2 / G_norm) * q2 + (2.0 * G_3 / G_norm) * q3;

       // % Calculate the innovation covariance S
       // S = H*P*transpose(H)+ magnetometer_R;
       // % Kalman gain
       // K = P*transpose(H)*inv(S);
       // % Save previous estimate so we can correct the quaternion sign
       // x_previous = x(:, k);
       // % Update the state estimate
       // x(:, k) = x_previous + K*v;

        // Modulate the accelerometer R matrix based on how far we are from the gravity vector
       // % OK, here comes another trick. We modulate the R matrix based on
       // % how far off it is from the expected G measurement. The idea is to
       // % trust it only in steady-state conditions.
       // % To do so, we first do a low-pass filter in the accelerometer
       // % measurements:
       // % RC = Ts*(1-alpha)/alpha
       // % For Ts = 0.02 and RC = 0.03 (5 Hz)

        accelerometer_filtered = 0.385869545095038 * sqrt(pow(imu_measure_ptr->ax, 2) + pow(imu_measure_ptr->ay, 2) + pow(imu_measure_ptr->az, 2)) + (1 - 0.385869545095038) * accelerometer_filtered;
        accelerometer_R = exp(20 * fabs(accelerometer_filtered - G_norm));
        if(accelerometer_R > 1000.0) accelerometer_R = 1000.0;
        PGMATRIX_DATA(filter_struct->pR_accelerometer, 1, 1) = FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_1_1 * accelerometer_R;
        PGMATRIX_DATA(filter_struct->pR_accelerometer, 2, 2) = FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_2_2 * accelerometer_R;
        PGMATRIX_DATA(filter_struct->pR_accelerometer, 3, 3) = FILTER_PARAMETERS_R_ACCELEROMETER_EKF_DECOUPLED_3_3 * accelerometer_R;

        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, filter_struct->pR_accelerometer, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do GPS:
    // %%% Fun�ao de medi�ao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
    if(gps_measure_ptr->FlagValidPositionMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
        // else
        //     H = [zeros(3,4) eye(3,3) zeros(3,3)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 3, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 3, 10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H, 1, 5) = 1.0;
        GMATRIX_DATA(H, 2, 6) = 1.0;
        GMATRIX_DATA(H, 3, 7) = 1.0;
        // v = (gpsmeasure.p - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, gps_measure_ptr->pPosition, MatDummy.pMat1);

        // K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_COPY(&R, gps_measure_ptr->pPPosition);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do GPS:
    // %%% Fun�ao de medi�ao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
    if(gps_measure_ptr->FlagValidVelocityMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
        // else
        //     H = [zeros(3,4) zeros(3,3) eye(3,3)];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 3, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 3, 10);
            GMATRIX_ZEROES(H);
        }
        GMATRIX_DATA(H, 1, 8) = 1.0;
        GMATRIX_DATA(H, 2, 9) = 1.0;
        GMATRIX_DATA(H, 3, 10) = 1.0;
        // v = (gpsmeasure.v - H*X);
        PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
        PGMATRIX_SUBTRACT_COPY(&V, gps_measure_ptr->pVelocity, MatDummy.pMat1);

        // K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
        // kf_structure.X = X + K*v;
        // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
        PGMATRIX_COPY(&R, gps_measure_ptr->pPVelocity);
        kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

        // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
        localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
    }

    // %%% corrige usando medidas do sonar:
    // %%% Fun�ao de medi�ao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
    // %%% Fun�ao de medi�ao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
    if(sonar_measure_ptr->FlagValidMeasure)
    {
        // X = kf_structure.X;
        PGMATRIX_COPY(&X_predicted, filter_struct->pX);

        // P = kf_structure.P;
        PGMATRIX_COPY(&P_predicted, filter_struct->pP);
        // dHdX1 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(1);
        // dHdX2 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(2);
        // dHdX3 =  X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(3);
        // dHdX4 = -X(7)*(((X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2))^(-2))*2*X(4);
        // dHdX7 =  1/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
        // if kf_structure.flagestimateaccelerometerbias
        //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0 0 0 0];
        // else
        //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 dHdX7 0 0 0];
        // end
        if (filter_struct->FlagEstimateAccelerometerBias)
        {
            GMATRIX_SETSIZE(H, 1, 13);
            GMATRIX_ZEROES(H);
        }
        else
        {
            GMATRIX_SETSIZE(H, 1, 10);
            GMATRIX_ZEROES(H);
        }

        aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,1,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,2,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,3,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,4,1));

        if(aux != 0.0)
        { // situa��o em que a medi��o seria diferente de infinito
            GMATRIX_DATA(H, 1, 1) = -GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 1, 1);
            GMATRIX_DATA(H, 1, 2) =  GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 2, 1);
            GMATRIX_DATA(H, 1, 3) =  GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 3, 1);
            GMATRIX_DATA(H, 1, 4) = -GMATRIX_DATA(X_predicted, 7, 1) * (1.0 / (GMATRIXMACRO_SQR(aux))) * 2.0 * GMATRIX_DATA(X_predicted, 4, 1);
            GMATRIX_DATA(H,1,7) =  1.0/(aux);

            // v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
            GMATRIX_SETSIZE(V, 1, 1);
            GMATRIX_DATA(V, 1, 1) = sonar_measure_ptr->range - GMATRIX_DATA(X_predicted, 7, 1)/aux;

            // K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
            // kf_structure.X = X + K*v;
            // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
            GMATRIX_SETSIZE(R, 1, 1);
            GMATRIX_DATA(R, 1, 1) = sonar_measure_ptr->rangevariance;
            kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);

            // kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
            localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);
        }
    }
    PGMATRIX_COPY(&X_predicted, filter_struct->pX);

    // P = kf_structure.P;
    PGMATRIX_COPY(&P_predicted, filter_struct->pP);
    // dHdX1 = 2*X(1);
    // dHdX2 = 2*X(2);
    // dHdX3 = 2*X(3);
    // dHdX4 = 2*X(4);
    // if kf_structure.flagestimateaccelerometerbias
    //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
    // else
    //     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
    // end
    if (filter_struct->FlagEstimateAccelerometerBias)
    {
        GMATRIX_SETSIZE(H, 1, 13);
        GMATRIX_ZEROES(H);
    }
    else
    {
        GMATRIX_SETSIZE(H, 1, 10);
        GMATRIX_ZEROES(H);
    }
    GMATRIX_DATA(H, 1, 1) = 2.0*GMATRIX_DATA(X_predicted, 1, 1);
    GMATRIX_DATA(H, 1, 2) = 2.0*GMATRIX_DATA(X_predicted, 2, 1);
    GMATRIX_DATA(H, 1, 3) = 2.0*GMATRIX_DATA(X_predicted, 3, 1);
    GMATRIX_DATA(H, 1, 4) = 2.0*GMATRIX_DATA(X_predicted, 4, 1);
    // v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
    GMATRIX_SETSIZE(V, 1, 1);
    aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 1, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 2, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 3, 1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted, 4, 1));
    GMATRIX_DATA(V,1,1) = 1.0 - aux;


    // K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
    // kf_structure.X = X + K*v;
    // kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
    GMATRIX_SETSIZE(R, 1, 1);
    PGMATRIX_COPY(&R, filter_struct->pR_pseudomeasurementnorm);
    kalman_EKF_update_innovationform(filter_struct->pX, &X_predicted, &V, filter_struct->pP, &P_predicted, &R, &H, &MatDummy, 1);
    localization_filter_quaternionscorrectsign(filter_struct, &X_predicted);

    // Retorna
    return 1;
}

int localization_filter_process_model_df_du(PGMATRIX pdf_du, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias)
{
    double sx, sy, sz, a, ca2, sa2;
    double x1, x2, x3, x4;
    double da_dwx, da_dwy, da_dwz;
    DCM_DECLARE(R);
    QUATERNIONS_DECLARE(q);

    sx = PGMATRIX_DATA(pu_imu,U_wx,1) * T;
    sy = PGMATRIX_DATA(pu_imu,U_wy,1) * T;
    sz = PGMATRIX_DATA(pu_imu,U_wz,1) * T;
    a = sqrt(sx*sx + sy*sy + sz*sz);

    ca2 = cos(a/2.0);
    sa2 = sin(a/2.0);

    // zera df_du:
    PGMATRIX_ZEROES(pdf_du);

    // %Atualiza a atitude
    if (a>0.0){
        // da_dwx = (1/2)*(1/a)*(2*sx*T);
        // da_dwy = (1/2)*(1/a)*(2*sy*T);
        // da_dwz = (1/2)*(1/a)*(2*sz*T);
        da_dwx = (0.5) * (1.0/a) * (2.0 * sx * T);
        da_dwy = (0.5) * (1.0/a) * (2.0 * sy * T);
        da_dwz = (0.5) * (1.0/a) * (2.0 * sz * T);

        x1 = PGMATRIX_DATA(pX_previous, X_q0, 1);
        x2 = PGMATRIX_DATA(pX_previous, X_q1, 1);
        x3 = PGMATRIX_DATA(pX_previous, X_q2, 1);
        x4 = PGMATRIX_DATA(pX_previous, X_q3, 1);

        // df_du(1,3) = -sa2*(1/2)*da_dwx*x1 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x2 ); % df_dwx
        // df_du(1,4) = -sa2*(1/2)*da_dwy*x1 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x3 ); % df_dwy
        // df_du(1,5) = -sa2*(1/2)*da_dwz*x1 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x4 ); % df_dwz

        PGMATRIX_DATA(pdf_du, 1, 3) = -sa2 * (0.5) * da_dwx * x1 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwx * ( sx*x2 + sy*x3 + sz*x4) - sa2 / a * ( T*x2 );
        PGMATRIX_DATA(pdf_du, 1, 4) = -sa2 * (0.5) * da_dwy * x1 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwy * ( sx*x2 + sy*x3 + sz*x4) - sa2 / a * ( T*x3 );
        PGMATRIX_DATA(pdf_du, 1, 5) = -sa2 * (0.5) * da_dwz * x1 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwz * ( sx*x2 + sy*x3 + sz*x4) - sa2 / a * ( T*x4 );

        // df_du(2, 3) = -sa2*(1/2)*da_dwx*x2 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*(-T*x1 ); % df_dwx
        // df_du(2, 4) = -sa2*(1/2)*da_dwy*x2 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*( T*x4 ); % df_dwy
        // df_du(2, 5) = -sa2*(1/2)*da_dwz*x2 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*(-T*x3 ); % df_dwz

        PGMATRIX_DATA(pdf_du, 2, 3) = -sa2 * (0.5) * da_dwx * x2 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwx * (-sx*x1 - sz*x3 + sy*x4) - sa2 / a * (-T*x1 );
        PGMATRIX_DATA(pdf_du, 2, 4) = -sa2 * (0.5) * da_dwy * x2 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwy * (-sx*x1 - sz*x3 + sy*x4) - sa2 / a * ( T*x4 );
        PGMATRIX_DATA(pdf_du, 2, 5) = -sa2 * (0.5) * da_dwz * x2 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwz * (-sx*x1 - sz*x3 + sy*x4) - sa2 / a * (-T*x3 );

        // df_du(3, 3) = -sa2*(1/2)*da_dwx*x3 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*(-T*x4 ); % df_dwx
        // df_du(3, 4) = -sa2*(1/2)*da_dwy*x3 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*(-T*x1 ); % df_dwy
        // df_du(3, 5) = -sa2*(1/2)*da_dwz*x3 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*( T*x2 ); % df_dwz

        PGMATRIX_DATA(pdf_du, 3, 3) = -sa2 * (0.5) * da_dwx * x3 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwx * (-sy*x1 + sz*x2 - sx*x4) - sa2 / a * (-T*x4 );
        PGMATRIX_DATA(pdf_du, 3, 4) = -sa2 * (0.5) * da_dwy * x3 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwy * (-sy*x1 + sz*x2 - sx*x4) - sa2 / a * (-T*x1 );
        PGMATRIX_DATA(pdf_du, 3, 5) = -sa2 * (0.5) * da_dwz * x3 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwz * (-sy*x1 + sz*x2 - sx*x4) - sa2 / a * ( T*x2 );

        // df_du(4, 3) = -sa2*(1/2)*da_dwx*x4 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*( T*x3 ); % df_dwx
        // df_du(4, 4) = -sa2*(1/2)*da_dwy*x4 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*(-T*x2 ); % df_dwy
        // df_du(4, 5) = -sa2*(1/2)*da_dwz*x4 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*(-T*x1 ); % df_dwz

        PGMATRIX_DATA(pdf_du, 4, 3) = -sa2 * (0.5) * da_dwx * x4 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwx * (-sz*x1 - sy*x2 + sx*x3) - sa2 / a * ( T*x3 );
        PGMATRIX_DATA(pdf_du, 4, 4) = -sa2 * (0.5) * da_dwy * x4 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwy * (-sz*x1 - sy*x2 + sx*x3) - sa2 / a * (-T*x2 );
        PGMATRIX_DATA(pdf_du, 4, 5) = -sa2 * (0.5) * da_dwz * x4 + (-ca2 * (0.5) * a + sa2) / (a*a) * da_dwz * (-sz*x1 - sy*x2 + sx*x3) - sa2 / a * (-T*x1 );

    }

    // R = quaternions2dcm(x(1:4));
    QUATERNIONS_Q0(q) = PGMATRIX_DATA(pX_previous, X_q0, 1);
    QUATERNIONS_Q1(q) = PGMATRIX_DATA(pX_previous, X_q1, 1);
    QUATERNIONS_Q2(q) = PGMATRIX_DATA(pX_previous, X_q2, 1);
    QUATERNIONS_Q3(q) = PGMATRIX_DATA(pX_previous, X_q3, 1);
    rotation_quaternions2dcm(&q, &R);

    // df_du(5:7,1:3)  = R*(T^2)/2;
    PGMATRIX_SUBMATRIX_COPY(pdf_du, 5, 1, &R);
    PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_du, 5, 7, 1, 3, (-0.5*T*T));

    // df_du(8:10,1:3)  = R*(T);
    PGMATRIX_SUBMATRIX_COPY(pdf_du, 8, 1, &R);
    PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_du, 8, 10, 1, 3, (T));

    // Retorna
    return 1; 
} 

int localization_filter_process_model_df_dx(PGMATRIX pdf_dx, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias)
{
    double sx, sy, sz, a, ca2, sa2;
    double x1, x2, x3, x4;
    DCM_DECLARE(R);
    QUATERNIONS_DECLARE(q);
    GMATRIX_DECLARE(dR_dxi, 3, 3);
    GMATRIX_DECLARE(da_dx, 3, 4);

    sx = PGMATRIX_DATA(pu_imu, U_wx, 1) * T;
    sy = PGMATRIX_DATA(pu_imu, U_wy, 1) * T;
    sz = PGMATRIX_DATA(pu_imu, U_wz, 1) * T;
    a = sqrt(sx*sx + sy*sy + sz*sz);

    ca2 = cos(a/2.0);
    sa2 = sin(a/2.0);

    // zera df_dx:
    PGMATRIX_ZEROES(pdf_dx);

    // %Atualiza a atitude
    if (a>0.0)
    {
        PGMATRIX_DATA(pdf_dx, 1, 1) = ca2;
        PGMATRIX_DATA(pdf_dx, 1, 2) = -(sa2 / a) * sx;
        PGMATRIX_DATA(pdf_dx, 1, 3) = -(sa2 / a) * sy;
        PGMATRIX_DATA(pdf_dx, 1, 4) = -(sa2 / a) * sz;

        PGMATRIX_DATA(pdf_dx, 2, 1) = (sa2 / a) * sx;
        PGMATRIX_DATA(pdf_dx, 2, 2) = ca2;
        PGMATRIX_DATA(pdf_dx, 2, 3) = (sa2 / a) * sz;
        PGMATRIX_DATA(pdf_dx, 2, 4) = -(sa2 / a) * sy;

        PGMATRIX_DATA(pdf_dx, 3, 1) = -(sa2 / a) * sy;
        PGMATRIX_DATA(pdf_dx, 3, 2) = (sa2 / a) * sz;
        PGMATRIX_DATA(pdf_dx, 3, 3) = ca2;
        PGMATRIX_DATA(pdf_dx, 3, 4) = -(sa2 / a) * sx;

        PGMATRIX_DATA(pdf_dx, 4, 1) = -(sa2 / a) * sz;
        PGMATRIX_DATA(pdf_dx, 4, 2) = -(sa2 / a) * sy;
        PGMATRIX_DATA(pdf_dx, 4, 3) = (sa2 / a) * sx;
        PGMATRIX_DATA(pdf_dx, 4, 4) = ca2;
    }
    else
    {
        PGMATRIX_DATA(pdf_dx, 1, 1) = 1.0;
        PGMATRIX_DATA(pdf_dx, 2, 2) = 1.0;
        PGMATRIX_DATA(pdf_dx, 3, 3) = 1.0;
        PGMATRIX_DATA(pdf_dx, 4, 4) = 1.0;
    }

    // R = quaternions2dcm(x(1:4));
    QUATERNIONS_Q0(q) = PGMATRIX_DATA(pX_previous, X_q0, 1);
    QUATERNIONS_Q1(q) = PGMATRIX_DATA(pX_previous, X_q1, 1);
    QUATERNIONS_Q2(q) = PGMATRIX_DATA(pX_previous, X_q2, 1);
    QUATERNIONS_Q3(q) = PGMATRIX_DATA(pX_previous, X_q3, 1);
    rotation_quaternions2dcm(&q, &R);

    x1 = PGMATRIX_DATA(pX_previous, X_q0, 1);
    x2 = PGMATRIX_DATA(pX_previous, X_q1, 1);
    x3 = PGMATRIX_DATA(pX_previous, X_q2, 1);
    x4 = PGMATRIX_DATA(pX_previous, X_q3, 1);


    // dR_dx1 = 2*[(2*x(1))/2    (-x(4))       (x(3));
    //          ( x(4))       (2*x(1))/2    (-x(2));
    //          (-x(3))       (x(2))        (2*x(1))/2];
    GMATRIX_DATA(dR_dxi, 1, 1) = (2 * x1);
    GMATRIX_DATA(dR_dxi, 1, 2) = 2.0 * (-x4);
    GMATRIX_DATA(dR_dxi, 1, 3) = 2.0 * (x3);
    GMATRIX_DATA(dR_dxi, 2, 1) = 2.0 * ( x4);
    GMATRIX_DATA(dR_dxi, 2, 2) = (2 * x1);
    GMATRIX_DATA(dR_dxi, 2, 3) = 2.0 * (-x2);
    GMATRIX_DATA(dR_dxi, 3, 1) = 2.0 * (-x3);
    GMATRIX_DATA(dR_dxi, 3, 2) = 2.0 * (x2);
    GMATRIX_DATA(dR_dxi, 3, 3) = (2 * x1);
    // da_dx(:, 1) = dR_dx1*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
    GMATRIX_DATA(da_dx, 1, 1) = GMATRIX_DATA(dR_dxi, 1, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 1, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 1, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 2, 1) = GMATRIX_DATA(dR_dxi, 2, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 2, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 2, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 3, 1) = GMATRIX_DATA(dR_dxi, 3, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 3, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 3, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);

    // dR_dx2 = 2*[(2*x(2))/2    (x(3))        (x(4));
    //          (x(3))        (-2*x(2))/2   (-x(1));
    //          (x(4))        (x(1))        (-2*x(2))/2];
    GMATRIX_DATA(dR_dxi, 1, 1) = (2 * x2);
    GMATRIX_DATA(dR_dxi, 1, 2) = 2.0 * (x3);
    GMATRIX_DATA(dR_dxi, 1, 3) = 2.0 * (x4);
    GMATRIX_DATA(dR_dxi, 2, 1) = 2.0 * (x3);
    GMATRIX_DATA(dR_dxi, 2, 2) = (-2 * x2);
    GMATRIX_DATA(dR_dxi, 2, 3) = 2.0 * (-x1);
    GMATRIX_DATA(dR_dxi, 3, 1) = 2.0 * (x4);
    GMATRIX_DATA(dR_dxi, 3, 2) = 2.0 * (x1);
    GMATRIX_DATA(dR_dxi, 3, 3) = (-2 * x2);

    // da_dx(:,2) = dR_dx2*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
    GMATRIX_DATA(da_dx, 1, 2) = GMATRIX_DATA(dR_dxi, 1, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 1, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 1, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 2, 2) = GMATRIX_DATA(dR_dxi, 2, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 2, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 2, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 3, 2) = GMATRIX_DATA(dR_dxi, 3, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 3, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 3, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);

    // dR_dx3 = 2*[(-2*x(3))/2   (x(2))        (x(1));
    //          (x(2))        (2*x(3))/2    (x(4));
    //          (-x(1))       (x(4))        (-2*x(3))/2];
    GMATRIX_DATA(dR_dxi, 1, 1) = (-2 * x3);
    GMATRIX_DATA(dR_dxi, 1, 2) = 2.0 * (x2);
    GMATRIX_DATA(dR_dxi, 1, 3) = 2.0 * (x1);
    GMATRIX_DATA(dR_dxi, 2, 1) = 2.0 * (x2);
    GMATRIX_DATA(dR_dxi, 2, 2) = (2 * x3);
    GMATRIX_DATA(dR_dxi, 2, 3) = 2.0 * (x4);
    GMATRIX_DATA(dR_dxi, 3, 1) = 2.0 * (-x1);
    GMATRIX_DATA(dR_dxi, 3, 2) = 2.0 * (x4);
    GMATRIX_DATA(dR_dxi, 3, 3) = (-2 * x3);

    // da_dx(:,3) = dR_dx3*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
    GMATRIX_DATA(da_dx, 1, 3) = GMATRIX_DATA(dR_dxi, 1, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 1, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 1, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 2, 3) = GMATRIX_DATA(dR_dxi, 2, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 2, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 2, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 3, 3) = GMATRIX_DATA(dR_dxi, 3, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 3, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 3, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);

    // dR_dx4 = 2*[(-2*x(4))/2   -x(1)         (x(2));
    //          (x(1))        (-2*x(4))/2    (x(3));
    //          (x(2))        (x(3))         (2*x(4))/2];
    GMATRIX_DATA(dR_dxi, 1, 1) = (-2 * x4);
    GMATRIX_DATA(dR_dxi, 1, 2) = 2.0 * (-x1);
    GMATRIX_DATA(dR_dxi, 1, 3) = 2.0 * (x2);
    GMATRIX_DATA(dR_dxi, 2, 1) = 2.0 * (x1);
    GMATRIX_DATA(dR_dxi, 2, 2) = (-2 * x4);
    GMATRIX_DATA(dR_dxi, 2, 3) = 2.0 * (x3);
    GMATRIX_DATA(dR_dxi, 3, 1) = 2.0 * (x2);
    GMATRIX_DATA(dR_dxi, 3, 2) = 2.0 * (x3);
    GMATRIX_DATA(dR_dxi, 3, 3) = (2 * x4);

    // da_dx(:,4) = dR_dx4*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
    GMATRIX_DATA(da_dx, 1, 4) = GMATRIX_DATA(dR_dxi, 1, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 1, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 1, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 2, 4) = GMATRIX_DATA(dR_dxi, 2, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 2, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 2, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);
    GMATRIX_DATA(da_dx, 3, 4) = GMATRIX_DATA(dR_dxi, 3, 1) * PGMATRIX_DATA(pu_imu, U_ax, 1) + GMATRIX_DATA(dR_dxi, 3, 2) * PGMATRIX_DATA(pu_imu, U_ay, 1) + GMATRIX_DATA(dR_dxi, 3, 3) * PGMATRIX_DATA(pu_imu, U_az, 1);

    // df_dx(5:7,1:4)  = da_dx*(T^2)/2;
    PGMATRIX_SUBMATRIX_COPY(pdf_dx, 5, 1, &da_dx);
    PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx, 5, 7, 1, 4, (0.5*T*T));

    // df_dx(5:7,5:7)  = eye(3);
    PGMATRIX_DATA(pdf_dx, 5, 5) = 1.0;
    PGMATRIX_DATA(pdf_dx, 6, 6) = 1.0;
    PGMATRIX_DATA(pdf_dx, 7, 7) = 1.0;

    // df_dx(5:7,8:10) = eye(3)*T;
    PGMATRIX_DATA(pdf_dx, 5, 8) = T;
    PGMATRIX_DATA(pdf_dx, 6, 9) = T;
    PGMATRIX_DATA(pdf_dx, 7, 10) = T;

    // if flagestimateaccelerometerbias
    //     df_dx(5:7,11:13) = -R*(T^2)/2;
    // end
    if (FlagEstimateAccelerometerBias)
    {
        PGMATRIX_SUBMATRIX_COPY(pdf_dx, 5, 11, &R);
        PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx, 5, 7, 11, 13, (-0.5*T*T));
    }

    // df_dx(8:10,1:4)  = da_dx*(T);
    PGMATRIX_SUBMATRIX_COPY(pdf_dx, 8, 1, &da_dx);
    PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx, 8, 10, 1, 4, (T));

    // df_dx(8:10,8:10) = eye(3);
    PGMATRIX_DATA(pdf_dx, 8, 8) = 1.0;
    PGMATRIX_DATA(pdf_dx, 9, 9) = 1.0;
    PGMATRIX_DATA(pdf_dx, 10, 10) = 1.0;

    // if flagestimateaccelerometerbias
    //     df_dx(8:10,11:13) = -R*(T);
    // end
    if (FlagEstimateAccelerometerBias)
    {
        PGMATRIX_SUBMATRIX_COPY(pdf_dx, 8, 11, &R);
        PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx, 8, 10, 11, 13, (-T));
    }

    // if flagestimateaccelerometerbias
    //     df_dx(11:13,11:13) = eye(3);
    // end
    if (FlagEstimateAccelerometerBias)
    {
        PGMATRIX_DATA(pdf_dx, 11, 11) = 1.0;
        PGMATRIX_DATA(pdf_dx, 12, 12) = 1.0;
        PGMATRIX_DATA(pdf_dx, 13, 13) = 1.0;
    }

    // Retorna
    return 1; 
} 

int localization_filter_process_model_evaluate(PGMATRIX pX, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias)
{
    double sx, sy, sz, a, ca2, sa2;
    double x1, x2, x3, x4;
    GMATRIX_DECLARE(acc_b, 3, 1);
    GMATRIX_DECLARE(acc_e, 3, 1);
    DCM_DECLARE(R);
    QUATERNIONS_DECLARE(q);

    sx = PGMATRIX_DATA(pu_imu, U_wx, 1) * T;
    sy = PGMATRIX_DATA(pu_imu, U_wy, 1) * T;
    sz = PGMATRIX_DATA(pu_imu, U_wz, 1) * T;
    a = sqrt(sx*sx + sy*sy + sz*sz);

    ca2 = cos(a/2.0);
    sa2 = sin(a/2.0);

    x1 = PGMATRIX_DATA(pX_previous, X_q0, 1);
    x2 = PGMATRIX_DATA(pX_previous, X_q1, 1);
    x3 = PGMATRIX_DATA(pX_previous, X_q2, 1);
    x4 = PGMATRIX_DATA(pX_previous, X_q3, 1);

    // %Atualiza a atitude
    if (a>0.0)
    {
        // q(1) = cos(a/2)*x(1) - sin(a/2)/a*( sx*x(2) + sy*x(3) + sz*x(4));
        PGMATRIX_DATA(pX, X_q0, 1) = ca2*x1 - (sa2/a) * (sx*x2 + sy*x3 + sz*x4);
        // q(2) = cos(a/2)*x(2) - sin(a/2)/a*(-sx*x(1) - sz*x(3) + sy*x(4));
        PGMATRIX_DATA(pX, X_q1, 1) = ca2*x2 - (sa2/a) * (-sx*x1 - sz*x3 + sy*x4);
        // q(3) = cos(a/2)*x(3) - sin(a/2)/a*(-sy*x(1) + sz*x(2) - sx*x(4));
        PGMATRIX_DATA(pX, X_q2, 1) = ca2*x3 - (sa2/a) * (-sy*x1 + sz*x2 - sx*x4);
        // q(4) = cos(a/2)*x(4) - sin(a/2)/a*(-sz*x(1) - sy*x(2) + sx*x(3));
        PGMATRIX_DATA(pX, X_q3, 1) = ca2*x4 - (sa2/a) * (-sz*x1 - sy*x2 + sx*x3);
    }
    else
    {
        PGMATRIX_DATA(pX, X_q0, 1) = x1;
        PGMATRIX_DATA(pX, X_q1, 1) = x2;
        PGMATRIX_DATA(pX, X_q2, 1) = x3;
        PGMATRIX_DATA(pX, X_q3, 1) = x4;
    }
    // %preserva a atitude, posi��o e velocidade pr�via
    // if flagestimateaccelerometerbias
    //     a = quaternions2dcm(x(1:4))*[imumeasure.ax-x(11); imumeasure.ay-x(12); imumeasure.az-x(13)] + G; % gravity compensation
    // else
    //     a = quaternions2dcm(x(1:4))*[imumeasure.ax; imumeasure.ay; imumeasure.az] + G; % gravity compensation
    // end
    QUATERNIONS_Q0(q) = PGMATRIX_DATA(pX_previous, X_q0, 1);
    QUATERNIONS_Q1(q) = PGMATRIX_DATA(pX_previous, X_q1, 1);
    QUATERNIONS_Q2(q) = PGMATRIX_DATA(pX_previous, X_q2, 1);
    QUATERNIONS_Q3(q) = PGMATRIX_DATA(pX_previous, X_q3, 1);
    rotation_quaternions2dcm(&q, &R);

    if (FlagEstimateAccelerometerBias)
    {
        GMATRIX_DATA(acc_b, 1, 1) = PGMATRIX_DATA(pu_imu, U_ax, 1) - PGMATRIX_DATA(pX_previous, X_bax, 1);
        GMATRIX_DATA(acc_b, 2, 1) = PGMATRIX_DATA(pu_imu, U_ay, 1) - PGMATRIX_DATA(pX_previous, X_bay, 1);
        GMATRIX_DATA(acc_b, 3, 1) = PGMATRIX_DATA(pu_imu, U_az, 1) - PGMATRIX_DATA(pX_previous, X_baz, 1);
    }
    else
    {
        GMATRIX_DATA(acc_b, 1, 1) = PGMATRIX_DATA(pu_imu, U_ax, 1);
        GMATRIX_DATA(acc_b, 2, 1) = PGMATRIX_DATA(pu_imu, U_ay, 1);
        GMATRIX_DATA(acc_b, 3, 1) = PGMATRIX_DATA(pu_imu, U_az, 1);
    }
    GMATRIX_DATA(acc_e, 1, 1) = GMATRIX_DATA(R, 1, 1) * GMATRIX_DATA(acc_b, 1, 1) + GMATRIX_DATA(R, 1, 2) * GMATRIX_DATA(acc_b, 2, 1) + GMATRIX_DATA(R, 1, 3) * GMATRIX_DATA(acc_b, 3, 1) + PGMATRIX_DATA(pG, 1, 1);
    GMATRIX_DATA(acc_e, 2, 1) = GMATRIX_DATA(R, 2, 1) * GMATRIX_DATA(acc_b, 1, 1) + GMATRIX_DATA(R, 2, 2) * GMATRIX_DATA(acc_b, 2, 1) + GMATRIX_DATA(R, 2, 3) * GMATRIX_DATA(acc_b, 3, 1) + PGMATRIX_DATA(pG, 2, 1);
    GMATRIX_DATA(acc_e, 3, 1) = GMATRIX_DATA(R, 3, 1) * GMATRIX_DATA(acc_b, 1, 1) + GMATRIX_DATA(R, 3, 2) * GMATRIX_DATA(acc_b, 2, 1) + GMATRIX_DATA(R, 3, 3) * GMATRIX_DATA(acc_b, 3, 1) + PGMATRIX_DATA(pG, 3, 1);

    // p = x(5:7) + x(8:10)*T + a*(T^2)/2;
    PGMATRIX_DATA(pX, X_px, 1) = PGMATRIX_DATA(pX_previous, X_px, 1) + PGMATRIX_DATA(pX_previous, X_vx, 1) * T + GMATRIX_DATA(acc_e, 1, 1) * (T*T) / 2.0;
    PGMATRIX_DATA(pX, X_py, 1) = PGMATRIX_DATA(pX_previous, X_py, 1) + PGMATRIX_DATA(pX_previous, X_vy, 1) * T + GMATRIX_DATA(acc_e, 2, 1) * (T*T) / 2.0;
    PGMATRIX_DATA(pX, X_pz, 1) = PGMATRIX_DATA(pX_previous, X_pz, 1) + PGMATRIX_DATA(pX_previous, X_vz, 1) * T + GMATRIX_DATA(acc_e, 3, 1) * (T*T) / 2.0;

    // v = x(8:10) + T*a;
    PGMATRIX_DATA(pX, X_vx, 1) = PGMATRIX_DATA(pX_previous, X_vx, 1) + GMATRIX_DATA(acc_e, 1, 1) * T;
    PGMATRIX_DATA(pX, X_vy, 1) = PGMATRIX_DATA(pX_previous, X_vy, 1) + GMATRIX_DATA(acc_e, 2, 1) * T;
    PGMATRIX_DATA(pX, X_vz, 1) = PGMATRIX_DATA(pX_previous, X_vz, 1) + GMATRIX_DATA(acc_e, 3, 1) * T;

    // % par�metros:
    // if flagestimateaccelerometerbias
    //     var_out(11:13,1) = x(11:13,1);
    // end
    if (FlagEstimateAccelerometerBias)
    {
        PGMATRIX_DATA(pX, X_bax, 1) = PGMATRIX_DATA(pX_previous, X_bax, 1);
        PGMATRIX_DATA(pX, X_bay, 1) = PGMATRIX_DATA(pX_previous, X_bay, 1);
        PGMATRIX_DATA(pX, X_baz, 1) = PGMATRIX_DATA(pX_previous, X_baz, 1);
    }

    // Retorna
    return 1; 
} 
/*****************************************************************************
*** int localization_triad(PQUATERNIONS pq, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, PQUATERNIONS pq_previous)
*** Entradas: 
*** Saidas:
*****************************************************************************/
int localization_triad(PQUATERNIONS pq, ImuMeasure *imu_measure_ptr, MagnetometerMeasure *magnetometer_measure_ptr, PGMATRIX pM, PGMATRIX pG, PQUATERNIONS pq_previous)
{
    double aux;
    GMATRIX_DECLARE(mag, 3, 1); // valores locais normalizados, para n�o alterar na origem
    GMATRIX_DECLARE(acc, 3, 1); // valores locais normalizados, para n�o alterar na origem
    GMATRIX_DECLARE(M, 3, 1); // valores locais normalizados, para n�o alterar na origem
    GMATRIX_DECLARE(G, 3, 1); // valores locais normalizados, para n�o alterar na origem

    GMATRIX_DECLARE(MatDummy, 3, 1);

    GMATRIX_DECLARE(I_b, 3, 1);
    GMATRIX_DECLARE(J_b, 3, 1);
    GMATRIX_DECLARE(K_b, 3, 1);
    GMATRIX_DECLARE(I_i, 3, 1);
    GMATRIX_DECLARE(J_i, 3, 1);
    GMATRIX_DECLARE(K_i, 3, 1);

    GMATRIX_DECLARE(R_i, 3, 3);
    GMATRIX_DECLARE(R_b, 3, 3);
    GMATRIX_DECLARE(R_b_i, 3, 3);

    // testar se medicoes do magnetometro s�o validas:
    if(!magnetometer_measure_ptr->FlagValidMeasure) return 0;

    // mag = [magnetometermeasure.mx; magnetometermeasure.my; magnetometermeasure.mz];
    GMATRIX_DATA(mag, 1, 1) = magnetometer_measure_ptr->mx;
    GMATRIX_DATA(mag, 2, 1) = magnetometer_measure_ptr->my;
    GMATRIX_DATA(mag, 3, 1) = magnetometer_measure_ptr->mz;
    // mag = (mag/sqrt(mag(1)^2 + mag(2)^2 + mag(3)^2));
    aux = PGMATRIX_NORM(&mag);
    if(aux == 0.0) return 0;

    PGMATRIX_MULTIPLY_CONST(&mag, (1.0/aux));

    // %armazena a acelera��o (for�a espec�fica). A acerela��o � negada, pois leva-se em conta que os
    // %acelerometros medem a for�a espec�fica sob sua massa e n�o a acelera��o
    // acc = zeros(3,1);
    // acc(1) = -imumeasure.ax;
    // acc(2) = -imumeasure.ay;
    // acc(3) = -imumeasure.az;
    GMATRIX_DATA(acc, 1, 1) = -imu_measure_ptr->ax;
    GMATRIX_DATA(acc, 2, 1) = -imu_measure_ptr->ay;
    GMATRIX_DATA(acc, 3, 1) = -imu_measure_ptr->az;
    // acc = acc/sqrt(acc(1)^2 + acc(2)^2 + acc(3)^2);
    aux = PGMATRIX_NORM(&acc);
    if(aux == 0.0) return 0;
    
    PGMATRIX_MULTIPLY_CONST(&acc, (1.0/aux));

    // G = G/sqrt(G(1)^2 + G(2)^2 + G(3)^2);
    PGMATRIX_COPY(&G, pG);
    aux = PGMATRIX_NORM(&G);
    if(aux == 0.0) return 0;

    PGMATRIX_MULTIPLY_CONST(&G, (1.0/aux));

    // M = M/sqrt(M(1)^2 + M(2)^2 + M(3)^2);
    PGMATRIX_COPY(&M, pM);
    aux = PGMATRIX_NORM(&M);
    if(aux == 0.0) return 0;

    PGMATRIX_MULTIPLY_CONST(&M, (1.0/aux));


    // Inicializa as variaveis do corpo. sendo "acc" as medidas do acelerometro e
    // "mag" as medidas do magnetometro. Lembrando que a fun��o cross(,)
    // representa produto vetorial entre dois vetores.
    // aux = acc + mag;
    GMATRIX_ADD_COPY(I_b, acc, mag);
    // I_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
    aux = PGMATRIX_NORM(&I_b);
    if(aux == 0.0) return 0;

    PGMATRIX_MULTIPLY_CONST(&I_b, (1.0/aux));


    // aux = cross(I_b, (acc - mag));
    GMATRIX_SUBTRACT_COPY(MatDummy, acc, mag);
    GMATRIX_CROSS_COPY(J_b, I_b, MatDummy);

    // J_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
    aux = PGMATRIX_NORM(&J_b);
    if(aux == 0.0) return 0;

    PGMATRIX_MULTIPLY_CONST(&J_b, (1.0/aux));

    // K_b = cross(I_b, J_b);
    GMATRIX_CROSS_COPY(K_b, I_b, J_b);

    // %Inicializa as variaveis do sistema inercial
    // aux = G + M;
    GMATRIX_ADD_COPY(I_i, G, M);
    // I_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
    aux = PGMATRIX_NORM(&I_i);
    if(aux == 0.0) return 0;

    PGMATRIX_MULTIPLY_CONST(&I_i, (1.0/aux));

    // aux = cross(I_i,(G - M));
    GMATRIX_SUBTRACT_COPY(MatDummy, G, M);
    GMATRIX_CROSS_COPY(J_i, I_i, MatDummy);

    // J_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
    aux = PGMATRIX_NORM(&J_i);
    if(aux == 0.0) return 0;

    PGMATRIX_MULTIPLY_CONST(&J_i, (1.0/aux));

    // K_i = cross(I_i,J_i);
    GMATRIX_CROSS_COPY(K_i, I_i, J_i);

    // Calculo da matriz de rota��o pelo m�todo TRIAD melhorado
    // R_i_b = ([I_b, J_b, K_b]*([I_i, J_i, K_i]'));
    GMATRIX_SUBMATRIX_COPY(R_b, 1, 1, I_b);
    GMATRIX_SUBMATRIX_COPY(R_b, 1, 2, J_b);
    GMATRIX_SUBMATRIX_COPY(R_b, 1, 3, K_b);
    GMATRIX_SUBMATRIX_COPY(R_i, 1, 1, I_i);
    GMATRIX_SUBMATRIX_COPY(R_i, 1, 2, J_i);
    GMATRIX_SUBMATRIX_COPY(R_i, 1, 3, K_i);
    GMATRIX_MULTIPLY_COPY_EXTENDED(R_b_i, R_i, 0, R_b, 1);

    // %Cria��o do quat�rnio de rota��o. (Deve se observar que a matriz de
    // %rota��o � do sistema de refer�ncia para o sistema do corpo a contr�rio do
    // %que o Padilha escreveu. Dados para a confer�ncia est�o em
    // %http://www.uel.br/proppg/semina/pdf/semina_28_1_22_19.pdf).
    // R_b_i = R_i_b';
    // q = dcm2quaternions(R_b_i);
    rotation_dcm2quaternions(&R_b_i, pq);

    // if FlagConsiderPrediction,
    //  q = quaternions_correctsign(q, q_predicted);
    //end
    if(pq_previous != NULL)
    {
        rotation_quaternionscorrectsign(pq, pq_previous);
    }

    // Retorna
    return 1; 
}                      
