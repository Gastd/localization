/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: localization.c 
Conteudo: Funções em código C relacionadas à localização.
Autor: G. A. Borges.
Atualizações: 
	- 19/09/2008: criação por Geovany A. Borges, modulo sem init/close
*****************************************************************************/
// Cabecalhos de bibliotecas C run-time
#include "mex.h" // apenas quando testar impressão no matlab
#include <math.h>
#include <stdio.h>

// Cabecalhos de modulos do projeto
#include "gmatrix.h"
#include "gmatrix_linalg.h"
#include "imu.h"
#include "gps.h"
#include "sonar.h"
#include "magnetometer.h"
#include "rotation.h"
#include "localization.h"
#include "kalman.h"
#include "unscentedtransform.h"

// Definicoes de uso interno

// Prototipos de funcoes internas ao modulo
int localization_filter_prediction_cekf(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PIMUMEASURE pIMUMeasure, PGMATRIX pG, double T);
int localization_filter_correction_cekf(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGPSMEASURE pGPSMeasure, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PSONARMEASURE pSonarMeasure, PGMATRIX pM, PGMATRIX pG, double T);
int localization_filter_prediction_ekf2(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PIMUMEASURE pIMUMeasure, PGMATRIX pG, double T);
int localization_filter_correction_ekf2(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGPSMEASURE pGPSMeasure, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PSONARMEASURE pSonarMeasure, PGMATRIX pM, PGMATRIX pG, double T);
int localization_filter_process_model_evaluate(PGMATRIX pX, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias);
int localization_filter_process_model_df_dx(PGMATRIX pdf_dx, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias);
int localization_filter_process_model_df_du(PGMATRIX pdf_du, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias);
int localization_converted_measurement_triad(PGMATRIX pYm, PGMATRIX pPy, PGMATRIX pX_predicted, PGMATRIX pP_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias);
int localization_filter_observation_model_triad_evaluate(PGMATRIX pYm, PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGMATRIX pX_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG);
void localization_filter_quaternionscorrectsign(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGMATRIX pX_predicted);
int localization_converted_measurement_triad_dg_du_imu(PGMATRIX pdg_du_imu, PGMATRIX pX_predicted, PGMATRIX pP_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias);

// Variaveis globais do módulo

/*****************************************************************************
******************************************************************************
** FUNCOES COM CHAMADA EXTERNA
******************************************************************************
*****************************************************************************/

/*****************************************************************************
*** int localization_init(int AlgorithmCode, int FlagEstimateAccelerometerBias, PLOCALIZATIONFILTERSTRUCT pFilterStruct)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int localization_init(int AlgorithmCode, int FlagEstimateAccelerometerBias, PLOCALIZATIONFILTERSTRUCT pFilterStruct)
{
	// Iniciar componentes da estrutura do filtro	
	pFilterStruct->AlgorithmCode = AlgorithmCode;

	pFilterStruct->FlagEstimateAccelerometerBias = FlagEstimateAccelerometerBias;

	pFilterStruct->Nstates = 4 + 3 + 3; // quaternions + posição 3D + velocidade 3D
	if (FlagEstimateAccelerometerBias){
		pFilterStruct->Nstates += 3; // + bias dos acelerometros.
	}

	pFilterStruct->pX = PGMATRIX_ALLOC(pFilterStruct->Nstates,1); // %should be set by user.
	pFilterStruct->pP = PGMATRIX_ALLOC(pFilterStruct->Nstates,pFilterStruct->Nstates); // %should be set by user.
	pFilterStruct->pPreset = PGMATRIX_ALLOC(pFilterStruct->Nstates,pFilterStruct->Nstates); // %should be set by user.
	
	if(pFilterStruct->AlgorithmCode==LOCALIZATION_ALGORITHMCODE_UKF2){
		pFilterStruct->pXsigma_ukf = PGMATRIX_ALLOC(pFilterStruct->Nstates,2*pFilterStruct->Nstates+1);
		pFilterStruct->pWsigma_ukf = PGMATRIX_ALLOC(1,2*pFilterStruct->Nstates+1);
	}

	pFilterStruct->pQ = PGMATRIX_ALLOC(pFilterStruct->Nstates,pFilterStruct->Nstates);
	PGMATRIX_ZEROES(pFilterStruct->pQ);
	PGMATRIX_DATA(pFilterStruct->pQ,1,1) = GMATRIXMACRO_SQR(0.01);
	PGMATRIX_DATA(pFilterStruct->pQ,2,2) = GMATRIXMACRO_SQR(0.01);
	PGMATRIX_DATA(pFilterStruct->pQ,3,3) = GMATRIXMACRO_SQR(0.01);
	PGMATRIX_DATA(pFilterStruct->pQ,4,4) = GMATRIXMACRO_SQR(0.01);

	PGMATRIX_DATA(pFilterStruct->pQ,5,5) = GMATRIXMACRO_SQR(0.01);
	PGMATRIX_DATA(pFilterStruct->pQ,6,6) = GMATRIXMACRO_SQR(0.01);
	PGMATRIX_DATA(pFilterStruct->pQ,7,7) = GMATRIXMACRO_SQR(0.01);

	PGMATRIX_DATA(pFilterStruct->pQ,8,8) = GMATRIXMACRO_SQR(0.0001);
	PGMATRIX_DATA(pFilterStruct->pQ,9,9) = GMATRIXMACRO_SQR(0.0001);
	PGMATRIX_DATA(pFilterStruct->pQ,10,10) = GMATRIXMACRO_SQR(0.0001);

	if (pFilterStruct->FlagEstimateAccelerometerBias){
		PGMATRIX_DATA(pFilterStruct->pQ,11,11) = GMATRIXMACRO_SQR(0.001);
		PGMATRIX_DATA(pFilterStruct->pQ,12,12) = GMATRIXMACRO_SQR(0.001);
		PGMATRIX_DATA(pFilterStruct->pQ,13,13) = GMATRIXMACRO_SQR(0.001);
	}

	pFilterStruct->pR_pseudomeasurementnorm = PGMATRIX_ALLOC(1,1);
	PGMATRIX_DATA(pFilterStruct->pR_pseudomeasurementnorm,1,1) = GMATRIXMACRO_SQR(1e-5);

	pFilterStruct->pR_convertedmeasurementtriad = PGMATRIX_ALLOC(4,4);
	PGMATRIX_ZEROES(pFilterStruct->pR_convertedmeasurementtriad);

	if(pFilterStruct->AlgorithmCode==LOCALIZATION_ALGORITHMCODE_CEKF){
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,1,1) = GMATRIXMACRO_SQR(0.001);
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,2,2) = GMATRIXMACRO_SQR(0.001);
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,3,3) = GMATRIXMACRO_SQR(0.001);
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,4,4) = GMATRIXMACRO_SQR(0.001);
	}
	else{
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,1,1) = GMATRIXMACRO_SQR(0.01);
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,2,2) = GMATRIXMACRO_SQR(0.01);
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,3,3) = GMATRIXMACRO_SQR(0.01);
		PGMATRIX_DATA(pFilterStruct->pR_convertedmeasurementtriad,4,4) = GMATRIXMACRO_SQR(0.01);
	}

	if(pFilterStruct->AlgorithmCode==LOCALIZATION_ALGORITHMCODE_CEKF){
		pFilterStruct->pS_previous_cekf = PGMATRIX_ALLOC(pFilterStruct->Nstates,4);
		PGMATRIX_ZEROES(pFilterStruct->pS_previous_cekf);
		pFilterStruct->pR_previous_cekf = PGMATRIX_ALLOC(4,4);
		PGMATRIX_IDENTITY(pFilterStruct->pR_previous_cekf);
		pFilterStruct->pA_imu_P_imu_cekf = PGMATRIX_ALLOC(pFilterStruct->Nstates,6);
		PGMATRIX_ZEROES(pFilterStruct->pA_imu_P_imu_cekf);
		pFilterStruct->pinnovation_previous_cekf = PGMATRIX_ALLOC(4,1);
		PGMATRIX_ZEROES(pFilterStruct->pinnovation_previous_cekf);
	}
	
	// Retorna 
    return 1; 
}

/*****************************************************************************
*** int localization_close(PLOCALIZATIONFILTERSTRUCT pFilterStruct)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int localization_close(PLOCALIZATIONFILTERSTRUCT pFilterStruct)
{
	// Desalocar componentes da estrutura do filtro	
	PGMATRIX_FREE(pFilterStruct->pX);
	PGMATRIX_FREE(pFilterStruct->pP); 
	PGMATRIX_FREE(pFilterStruct->pPreset);
	
	if(pFilterStruct->AlgorithmCode==LOCALIZATION_ALGORITHMCODE_UKF2){
		PGMATRIX_FREE(pFilterStruct->pXsigma_ukf);
		PGMATRIX_FREE(pFilterStruct->pWsigma_ukf);
	}

	PGMATRIX_FREE(pFilterStruct->pQ);

	PGMATRIX_FREE(pFilterStruct->pR_pseudomeasurementnorm);

	PGMATRIX_FREE(pFilterStruct->pR_convertedmeasurementtriad);

	if(pFilterStruct->AlgorithmCode==LOCALIZATION_ALGORITHMCODE_CEKF){
		PGMATRIX_FREE(pFilterStruct->pS_previous_cekf);
		PGMATRIX_FREE(pFilterStruct->pR_previous_cekf);
		PGMATRIX_FREE(pFilterStruct->pA_imu_P_imu_cekf);
		PGMATRIX_FREE(pFilterStruct->pinnovation_previous_cekf);
	}

	// Retorna 
    return 1; 
}

/*****************************************************************************
*** int localization_filter_correction(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGPSMEASURE pGPSMeasure, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PSONARMEASURE pSonarMeasure, PGMATRIX pM, PGMATRIX pG, double T)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int localization_filter_correction(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGPSMEASURE pGPSMeasure, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PSONARMEASURE pSonarMeasure, PGMATRIX pM, PGMATRIX pG, double T)
{
	// Correção conforme a arquitetura
	switch(pFilterStruct->AlgorithmCode){
	case LOCALIZATION_ALGORITHMCODE_EKF2:
		if(!localization_filter_correction_ekf2(pFilterStruct, pGPSMeasure, pIMUMeasure, pMagnetometerMeasure, pSonarMeasure, pM, pG, T)) return 0;
		break;

	case LOCALIZATION_ALGORITHMCODE_CEKF:
		if(!localization_filter_correction_cekf(pFilterStruct, pGPSMeasure, pIMUMeasure, pMagnetometerMeasure, pSonarMeasure, pM, pG, T)) return 0;
		break;

	case LOCALIZATION_ALGORITHMCODE_UKF2:
//		if(!localization_filter_correction_ukf2(pFilterStruct, pGPSMeasure, pIMUMeasure, pMagnetometerMeasure, pSonarMeasure, pM, pG, T)) return 0;
		break;
	}

	// Reset em caso de matriz P negativa:
/*	if(sum(eig(kf_structure.P)<0)>0)
		% filter reset:
		disp('*** Warning: Kalman filter reset');
		kf_structure.P = kf_structure.Preset;
	end
*/
	// Retorna 
    return 1; 
}


void localization_filter_quaternionscorrectsign(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGMATRIX pX_predicted)
{

	QUATERNIONS_DECLARE(qcurrent);
	QUATERNIONS_DECLARE(qpredicted);

	QUATERNIONS_Q0(qcurrent) = PGMATRIX_DATA(pFilterStruct->pX,X_q0,1);
	QUATERNIONS_Q1(qcurrent) = PGMATRIX_DATA(pFilterStruct->pX,X_q1,1);
	QUATERNIONS_Q2(qcurrent) = PGMATRIX_DATA(pFilterStruct->pX,X_q2,1);
	QUATERNIONS_Q3(qcurrent) = PGMATRIX_DATA(pFilterStruct->pX,X_q3,1);
	QUATERNIONS_Q0(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q0,1);
	QUATERNIONS_Q1(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q1,1);
	QUATERNIONS_Q2(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q2,1);
	QUATERNIONS_Q3(qpredicted) = PGMATRIX_DATA(pX_predicted,X_q3,1);
	rotation_quaternionscorrectsign(&qcurrent, &qpredicted);
	PGMATRIX_DATA(pFilterStruct->pX,X_q0,1) = QUATERNIONS_Q0(qcurrent);
	PGMATRIX_DATA(pFilterStruct->pX,X_q1,1) = QUATERNIONS_Q1(qcurrent);
	PGMATRIX_DATA(pFilterStruct->pX,X_q2,1) = QUATERNIONS_Q2(qcurrent);
	PGMATRIX_DATA(pFilterStruct->pX,X_q3,1) = QUATERNIONS_Q3(qcurrent);
}

int localization_filter_correction_cekf(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGPSMEASURE pGPSMeasure, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PSONARMEASURE pSonarMeasure, PGMATRIX pM, PGMATRIX pG, double T)
{
// %%% corrige usando medidas do algoritmo TRIAD: 
// %%% Funçao de mediçao: [eye(4) zeros(4,6)]*X (a medida é dada pelo estimador TRIAD.
// %%% Estimar a mediçao do quaternion assim como da matriz R usando UT;

	GMATRIX_DECLARE(dg_du_imu,4,6);
	GMATRIX_DECLARE(H,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(V,4,1);
	GMATRIX_DECLARE(Ym,4,1);
	GMATRIX_DECLARE(Py,4,4);
	GMATRIX_DECLARE(Rtilde,4,4);
	GMATRIX_DECLARE(Stilde,LOCALIZATION_MAXSTATESIZE,4);
	GMATRIX_DECLARE(R,4,4);
	GMATRIX_DECLARE(X_predicted,LOCALIZATION_MAXSTATESIZE,1);
	GMATRIX_DECLARE(P_predicted,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy1,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy2,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy3,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy4,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	DUMMY_MATRICES MatDummy;
	QUATERNIONS_DECLARE(qcurrent);
	QUATERNIONS_DECLARE(qpredicted);

	double aux;

	MatDummy.pMat1 = &MatDummy1;
	MatDummy.pMat2 = &MatDummy2;
	MatDummy.pMat3 = &MatDummy3;
	MatDummy.pMat4 = &MatDummy4;

	// Definição dos tamanhos das matrizes deve ser feito aqui, o que evita alocação dimâmica
	GMATRIX_SETSIZE(X_predicted,pFilterStruct->Nstates,1);
	GMATRIX_SETSIZE(P_predicted,pFilterStruct->Nstates,pFilterStruct->Nstates);
	GMATRIX_SETSIZE(Stilde,pFilterStruct->Nstates,1);

	// Medição pelo TRIAD
	if(pMagnetometerMeasure->FlagValidMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);

		// if kf_structure.flagestimateaccelerometerbias
		//     H = [eye(4) zeros(4,9)];
		// else
		//     H = [eye(4) zeros(4,6)];
		// end
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,4,13);
			GMATRIX_ZEROES(H);
		}
		else{	
			GMATRIX_SETSIZE(H,4,10);
			GMATRIX_ZEROES(H);
		}
		GMATRIX_DATA(H,1,1) = 1.0;
		GMATRIX_DATA(H,2,2) = 1.0;
		GMATRIX_DATA(H,3,3) = 1.0;
		GMATRIX_DATA(H,4,4) = 1.0;

	    // [Ym,Py] = ConvertedMeasurementTRIAD('cekf',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
		GMATRIX_SETSIZE(Ym,4,1);
		GMATRIX_SETSIZE(Py,4,4);
		localization_converted_measurement_triad(&Ym, &Py, &X_predicted, &P_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG, pFilterStruct->FlagEstimateAccelerometerBias);
//		GMATRIX_PRINT_MATLABFORM(Ym);
//		GMATRIX_PRINT_MATLABFORM(Py);
		
		// v = (Ym - H*X);
		PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
		PGMATRIX_SUBTRACT_COPY(&V, &Ym, MatDummy.pMat1);

		// Rtilde = Py + kf_structure.R_convertedmeasurementtriad_cekf;
		PGMATRIX_ADD_COPY(&Rtilde, &Py, pFilterStruct->pR_convertedmeasurementtriad);

		// Stilde = kf_structure.A_imu_P_imu_cekf * localization_filter_model_gtriad('dg_du_imu',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias,1)';
		localization_converted_measurement_triad_dg_du_imu(&dg_du_imu, &X_predicted, &P_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG, pFilterStruct->FlagEstimateAccelerometerBias);
		PGMATRIX_MULTIPLY_COPY_EXTENDED(&Stilde, pFilterStruct->pA_imu_P_imu_cekf, 0, &dg_du_imu, 1);
	
		// K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad_ekf);
		// kf_structure.X = X + K*v;
		// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
		PGMATRIX_ADD_COPY(&R, &Py, pFilterStruct->pR_convertedmeasurementtriad);
		kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

		// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
		localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);

		// kf_structure.R_previous_cekf = Rtilde;
		PGMATRIX_COPY(pFilterStruct->pR_previous_cekf,&Rtilde);

		// kf_structure.S_previous_cekf = Stilde;
		PGMATRIX_COPY(pFilterStruct->pS_previous_cekf,&Stilde);

		// kf_structure.innovation_previous_cekf = (Ym - H*kf_structure.X);
		PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, pFilterStruct->pX);
		PGMATRIX_SUBTRACT_COPY(&V, &Ym, MatDummy.pMat1);
		PGMATRIX_COPY(pFilterStruct->pinnovation_previous_cekf,&V);
	}

	// %%% corrige usando medidas do GPS:
	// %%% Funçao de mediçao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
	if(pGPSMeasure->FlagValidPositionMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
		// if kf_structure.flagestimateaccelerometerbias
		//     H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
		// else
		//     H = [zeros(3,4) eye(3,3) zeros(3,3)];
		// end
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,3,13);
			GMATRIX_ZEROES(H);
		}
		else{
			GMATRIX_SETSIZE(H,3,10);
			GMATRIX_ZEROES(H);
		}
		GMATRIX_DATA(H,1,5) = 1.0;
		GMATRIX_DATA(H,2,6) = 1.0;
		GMATRIX_DATA(H,3,7) = 1.0;
		// v = (gpsmeasure.p - H*X);
		PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
		PGMATRIX_SUBTRACT_COPY(&V, pGPSMeasure->pPosition, MatDummy.pMat1);

		// K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
		// kf_structure.X = X + K*v;
		// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
		PGMATRIX_COPY(&R,pGPSMeasure->pPPosition);
		kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

		// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
		localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);
	}

// %%% corrige usando medidas do GPS:
// %%% Funçao de mediçao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
	if(pGPSMeasure->FlagValidVelocityMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
		// if kf_structure.flagestimateaccelerometerbias
		//     H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
		// else
		//     H = [zeros(3,4) zeros(3,3) eye(3,3)];
		// end
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,3,13);
			GMATRIX_ZEROES(H);
		}
		else{
			GMATRIX_SETSIZE(H,3,10);
			GMATRIX_ZEROES(H);
		}
		GMATRIX_DATA(H,1,8) = 1.0;
		GMATRIX_DATA(H,2,9) = 1.0;
		GMATRIX_DATA(H,3,10) = 1.0;
		// v = (gpsmeasure.v - H*X);
		PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
		PGMATRIX_SUBTRACT_COPY(&V, pGPSMeasure->pVelocity, MatDummy.pMat1);

		// K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
		// kf_structure.X = X + K*v;
		// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
		PGMATRIX_COPY(&R,pGPSMeasure->pPVelocity);
		kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

		// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
		localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);
	}

// %%% corrige usando medidas do sonar:
// %%% Funçao de mediçao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
// %%% Funçao de mediçao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
	if(pSonarMeasure->FlagValidMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
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
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,1,13);
			GMATRIX_ZEROES(H);
		}
		else{
			GMATRIX_SETSIZE(H,1,10);
			GMATRIX_ZEROES(H);
		}

		aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,1,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,2,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,3,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,4,1));

		if(aux != 0.0){ // situação em que a medição seria diferente de infinito
			GMATRIX_DATA(H,1,1) = -GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,1,1);
			GMATRIX_DATA(H,1,2) =  GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,2,1);
			GMATRIX_DATA(H,1,3) =  GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,3,1);
			GMATRIX_DATA(H,1,4) = -GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,4,1);
			GMATRIX_DATA(H,1,7) =  1.0/(aux);

			// v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
			GMATRIX_SETSIZE(V,1,1);
			GMATRIX_DATA(V,1,1) = pSonarMeasure->range - GMATRIX_DATA(X_predicted,7,1)/aux;

			// K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
			// kf_structure.X = X + K*v;
			// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
			GMATRIX_SETSIZE(R,1,1);
			GMATRIX_DATA(R,1,1) = pSonarMeasure->rangevariance;
			kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

			// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
			localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);
		}
	}
// %%% corrige a norma do quaternion usando pseudo-mediçoes:
// %%% Funçao de mediçao: 1 = q0^2 + q1^2 + q2^2 + q3^2
// %%% Funçao de mediçao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
	// X = kf_structure.X;
	PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
	
	// P = kf_structure.P;
	PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
	// dHdX1 = 2*X(1);
	// dHdX2 = 2*X(2);
	// dHdX3 = 2*X(3);
	// dHdX4 = 2*X(4);
	// if kf_structure.flagestimateaccelerometerbias
	//     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
	// else
	//     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
	// end
	if (pFilterStruct->FlagEstimateAccelerometerBias){
		GMATRIX_SETSIZE(H,1,13);
		GMATRIX_ZEROES(H);
	}
	else{
		GMATRIX_SETSIZE(H,1,10);
		GMATRIX_ZEROES(H);
	}
	GMATRIX_DATA(H,1,1) = 2.0*GMATRIX_DATA(X_predicted,1,1);
	GMATRIX_DATA(H,1,2) = 2.0*GMATRIX_DATA(X_predicted,2,1);
	GMATRIX_DATA(H,1,3) = 2.0*GMATRIX_DATA(X_predicted,3,1);
	GMATRIX_DATA(H,1,4) = 2.0*GMATRIX_DATA(X_predicted,4,1);
	// v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
	GMATRIX_SETSIZE(V,1,1);
	aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,1,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,2,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,3,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,4,1));
	GMATRIX_DATA(V,1,1) = 1.0 - aux;


	// K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
	// kf_structure.X = X + K*v;
	// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
	GMATRIX_SETSIZE(R,1,1);
	PGMATRIX_COPY(&R,pFilterStruct->pR_pseudomeasurementnorm);
	kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);
	localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);

	// Retorna 
    return 1; 
}

int localization_filter_correction_ekf2(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGPSMEASURE pGPSMeasure, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PSONARMEASURE pSonarMeasure, PGMATRIX pM, PGMATRIX pG, double T)
{
// %%% corrige usando medidas do algoritmo TRIAD: 
// %%% Funçao de mediçao: [eye(4) zeros(4,6)]*X (a medida é dada pelo estimador TRIAD.
// %%% Estimar a mediçao do quaternion assim como da matriz R usando UT;

	GMATRIX_DECLARE(H,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(V,4,1);
	GMATRIX_DECLARE(Ym,4,1);
	GMATRIX_DECLARE(Py,4,4);
	GMATRIX_DECLARE(R,4,4);
	GMATRIX_DECLARE(X_predicted,LOCALIZATION_MAXSTATESIZE,1);
	GMATRIX_DECLARE(P_predicted,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy1,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy2,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy3,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy4,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	DUMMY_MATRICES MatDummy;
	QUATERNIONS_DECLARE(qcurrent);
	QUATERNIONS_DECLARE(qpredicted);
	double aux;

	MatDummy.pMat1 = &MatDummy1;
	MatDummy.pMat2 = &MatDummy2;
	MatDummy.pMat3 = &MatDummy3;
	MatDummy.pMat4 = &MatDummy4;

	// Definição dos tamanhos das matrizes deve ser feito aqui, o que evita alocação dimâmica
	GMATRIX_SETSIZE(X_predicted,pFilterStruct->Nstates,1);
	GMATRIX_SETSIZE(P_predicted,pFilterStruct->Nstates,pFilterStruct->Nstates);

	// Medição pelo TRIAD
	if(pMagnetometerMeasure->FlagValidMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);

		// if kf_structure.flagestimateaccelerometerbias
		//     H = [eye(4) zeros(4,9)];
		// else
		//     H = [eye(4) zeros(4,6)];
		// end
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,4,13);
			GMATRIX_ZEROES(H);
		}
		else{
			GMATRIX_SETSIZE(H,4,10);
			GMATRIX_ZEROES(H);
		}
		GMATRIX_DATA(H,1,1) = 1.0;
		GMATRIX_DATA(H,2,2) = 1.0;
		GMATRIX_DATA(H,3,3) = 1.0;
		GMATRIX_DATA(H,4,4) = 1.0;

		// [Ym,Py] = ConvertedMeasurementTRIAD('ekf2',X, P, imumeasure, magnetometermeasure, M, G, kf_structure.flagestimateaccelerometerbias);
		GMATRIX_SETSIZE(Ym,4,1);
		GMATRIX_SETSIZE(Py,4,4);
		localization_converted_measurement_triad(&Ym, &Py, &X_predicted, &P_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG, pFilterStruct->FlagEstimateAccelerometerBias);
//		GMATRIX_PRINT_MATLABFORM(Ym);
//		GMATRIX_PRINT_MATLABFORM(Py);
		
		// v = (Ym - H*X);
		PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
		PGMATRIX_SUBTRACT_COPY(&V, &Ym, MatDummy.pMat1);
		
		// K = P*(H')*inv(H*P*(H') + Py + kf_structure.R_convertedmeasurementtriad_ekf);
		// kf_structure.X = X + K*v;
		// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
		PGMATRIX_ADD_COPY(&R, &Py, pFilterStruct->pR_convertedmeasurementtriad);
		kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

		// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
		localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);
	}

	// %%% corrige usando medidas do GPS:
	// %%% Funçao de mediçao: [zeros(3,4) eye(3,3) zeros(3,3)]*X
	if(pGPSMeasure->FlagValidPositionMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
		// if kf_structure.flagestimateaccelerometerbias
		//     H = [zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3)];
		// else
		//     H = [zeros(3,4) eye(3,3) zeros(3,3)];
		// end
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,3,13);
			GMATRIX_ZEROES(H);
		}
		else{
			GMATRIX_SETSIZE(H,3,10);
			GMATRIX_ZEROES(H);
		}
		GMATRIX_DATA(H,1,5) = 1.0;
		GMATRIX_DATA(H,2,6) = 1.0;
		GMATRIX_DATA(H,3,7) = 1.0;
		// v = (gpsmeasure.p - H*X);
		PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
		PGMATRIX_SUBTRACT_COPY(&V, pGPSMeasure->pPosition, MatDummy.pMat1);

		// K = P*(H')*inv(H*P*(H') + gpsmeasure.P_p);
		// kf_structure.X = X + K*v;
		// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
		PGMATRIX_COPY(&R,pGPSMeasure->pPPosition);
		kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

		// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
		localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);
	}

// %%% corrige usando medidas do GPS:
// %%% Funçao de mediçao: [zeros(3,4) zeros(3,3) eye(3,3)]*X
	if(pGPSMeasure->FlagValidVelocityMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
		// if kf_structure.flagestimateaccelerometerbias
		//     H = [zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3)];
		// else
		//     H = [zeros(3,4) zeros(3,3) eye(3,3)];
		// end
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,3,13);
			GMATRIX_ZEROES(H);
		}
		else{
			GMATRIX_SETSIZE(H,3,10);
			GMATRIX_ZEROES(H);
		}
		GMATRIX_DATA(H,1,8) = 1.0;
		GMATRIX_DATA(H,2,9) = 1.0;
		GMATRIX_DATA(H,3,10) = 1.0;
		// v = (gpsmeasure.v - H*X);
		PGMATRIX_MULTIPLY_COPY(MatDummy.pMat1, &H, &X_predicted);
		PGMATRIX_SUBTRACT_COPY(&V, pGPSMeasure->pVelocity, MatDummy.pMat1);

		// K = P*(H')*inv(H*P*(H') + gpsmeasure.P_v);
		// kf_structure.X = X + K*v;
		// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
		PGMATRIX_COPY(&R,pGPSMeasure->pPVelocity);
		kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

		// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
		localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);
	}

// %%% corrige usando medidas do sonar:
// %%% Funçao de mediçao: 2z/(q0^2 - q1^2 - q2^2 + q3^2)
// %%% Funçao de mediçao codificada em estado: X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2);
	if(pSonarMeasure->FlagValidMeasure){
		// X = kf_structure.X;
		PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
		
		// P = kf_structure.P;
		PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
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
		if (pFilterStruct->FlagEstimateAccelerometerBias){
			GMATRIX_SETSIZE(H,1,13);
			GMATRIX_ZEROES(H);
		}
		else{
			GMATRIX_SETSIZE(H,1,10);
			GMATRIX_ZEROES(H);
		}

		aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,1,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,2,1)) - GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,3,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,4,1));

		if(aux != 0.0){ // situação em que a medição seria diferente de infinito
			GMATRIX_DATA(H,1,1) = -GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,1,1);
			GMATRIX_DATA(H,1,2) =  GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,2,1);
			GMATRIX_DATA(H,1,3) =  GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,3,1);
			GMATRIX_DATA(H,1,4) = -GMATRIX_DATA(X_predicted,7,1)*(1.0/(GMATRIXMACRO_SQR(aux)))*2.0*GMATRIX_DATA(X_predicted,4,1);
			GMATRIX_DATA(H,1,7) =  1.0/(aux);

			// v = (sonarmeasure.range - (X(7)/(X(1)^2 - X(2)^2 - X(3)^2 + X(4)^2)));
			GMATRIX_SETSIZE(V,1,1);
			GMATRIX_DATA(V,1,1) = pSonarMeasure->range - GMATRIX_DATA(X_predicted,7,1)/aux;

			// K = P*(H')*inv(H*P*(H') + sonarmeasure.rangevariance);
			// kf_structure.X = X + K*v;
			// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
			GMATRIX_SETSIZE(R,1,1);
			GMATRIX_DATA(R,1,1) = pSonarMeasure->rangevariance;
			kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);

			// kf_structure.X(1:4) = quaternions_correctsign(kf_structure.X(1:4), X(1:4));
			localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);
		}
	}
// %%% corrige a norma do quaternion usando pseudo-mediçoes:
// %%% Funçao de mediçao: 1 = q0^2 + q1^2 + q2^2 + q3^2
// %%% Funçao de mediçao codificada em estado: X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2;
	// X = kf_structure.X;
	PGMATRIX_COPY(&X_predicted, pFilterStruct->pX);
	
	// P = kf_structure.P;
	PGMATRIX_COPY(&P_predicted, pFilterStruct->pP);
	// dHdX1 = 2*X(1);
	// dHdX2 = 2*X(2);
	// dHdX3 = 2*X(3);
	// dHdX4 = 2*X(4);
	// if kf_structure.flagestimateaccelerometerbias
	//     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0 0 0 0];
	// else
	//     H = [dHdX1 dHdX2 dHdX3 dHdX4 0 0 0 0 0 0];
	// end
	if (pFilterStruct->FlagEstimateAccelerometerBias){
		GMATRIX_SETSIZE(H,1,13);
		GMATRIX_ZEROES(H);
	}
	else{
		GMATRIX_SETSIZE(H,1,10);
		GMATRIX_ZEROES(H);
	}
	GMATRIX_DATA(H,1,1) = 2.0*GMATRIX_DATA(X_predicted,1,1);
	GMATRIX_DATA(H,1,2) = 2.0*GMATRIX_DATA(X_predicted,2,1);
	GMATRIX_DATA(H,1,3) = 2.0*GMATRIX_DATA(X_predicted,3,1);
	GMATRIX_DATA(H,1,4) = 2.0*GMATRIX_DATA(X_predicted,4,1);
	// v = (1 - (X(1)^2 + X(2)^2 + X(3)^2 + X(4)^2));
	GMATRIX_SETSIZE(V,1,1);
	aux = GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,1,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,2,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,3,1)) + GMATRIXMACRO_SQR(GMATRIX_DATA(X_predicted,4,1));
	GMATRIX_DATA(V,1,1) = 1.0 - aux;


	// K = P*(H')*inv(H*P*(H') + kf_structure.R_pseudomeasurementnorm);
	// kf_structure.X = X + K*v;
	// kf_structure.P = (eye(kf_structure.Nstates) - K*H)*P;
	GMATRIX_SETSIZE(R,1,1);
	PGMATRIX_COPY(&R,pFilterStruct->pR_pseudomeasurementnorm);
	kalman_EKF_update_innovationform(pFilterStruct->pX, &X_predicted, &V, pFilterStruct->pP, &P_predicted, &R, &H, &MatDummy, 1);
	localization_filter_quaternionscorrectsign(pFilterStruct, &X_predicted);

	// Retorna 
    return 1; 
}



void localization_converted_measurement_triad_UT_nonbiased(PGMATRIX pY, PGMATRIX pX, PGMATRIX pParameters)
{
	QUATERNIONS_DECLARE(q_predicted);
	IMUMEASURE IMUMeasure;
	MAGNETOMETERMEASURE MagnetometerMeasure;
	GMATRIX_DECLARE(M,3,1);
	GMATRIX_DECLARE(G,3,1);

	QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX,X_q0,1);
	QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX,X_q1,1);
	QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX,X_q2,1);
	QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX,X_q3,1);

	IMUMeasure.ax = PGMATRIX_DATA(pX,X_vz+1,1);
	IMUMeasure.ay = PGMATRIX_DATA(pX,X_vz+2,1);
	IMUMeasure.az = PGMATRIX_DATA(pX,X_vz+3,1);

	MagnetometerMeasure.mx = PGMATRIX_DATA(pX,X_vz+4,1);
	MagnetometerMeasure.my = PGMATRIX_DATA(pX,X_vz+5,1);
	MagnetometerMeasure.mz = PGMATRIX_DATA(pX,X_vz+6,1);

	GMATRIX_DATA(M,1,1) = PGMATRIX_DATA(pParameters,1,1);
	GMATRIX_DATA(M,2,1) = PGMATRIX_DATA(pParameters,2,1);
	GMATRIX_DATA(M,3,1) = PGMATRIX_DATA(pParameters,3,1);
	GMATRIX_DATA(G,1,1) = PGMATRIX_DATA(pParameters,4,1);
	GMATRIX_DATA(G,2,1) = PGMATRIX_DATA(pParameters,5,1);
	GMATRIX_DATA(G,3,1) = PGMATRIX_DATA(pParameters,6,1);

	localization_triad(pY, &IMUMeasure, &MagnetometerMeasure, &M, &G, &q_predicted);
}

void localization_converted_measurement_triad_UT_biased(PGMATRIX pY, PGMATRIX pX, PGMATRIX pParameters)
{
	QUATERNIONS_DECLARE(q_predicted);
	IMUMEASURE IMUMeasure;
	MAGNETOMETERMEASURE MagnetometerMeasure;
	GMATRIX_DECLARE(M,3,1);
	GMATRIX_DECLARE(G,3,1);

	QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX,X_q0,1);
	QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX,X_q1,1);
	QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX,X_q2,1);
	QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX,X_q3,1);

	IMUMeasure.ax = PGMATRIX_DATA(pX,X_baz+1,1) - PGMATRIX_DATA(pX,X_bax,1);
	IMUMeasure.ay = PGMATRIX_DATA(pX,X_baz+2,1) - PGMATRIX_DATA(pX,X_bay,1);
	IMUMeasure.az = PGMATRIX_DATA(pX,X_baz+3,1) - PGMATRIX_DATA(pX,X_baz,1);

	MagnetometerMeasure.mx = PGMATRIX_DATA(pX,X_baz+4,1);
	MagnetometerMeasure.my = PGMATRIX_DATA(pX,X_baz+5,1);
	MagnetometerMeasure.mz = PGMATRIX_DATA(pX,X_baz+6,1);

	GMATRIX_DATA(M,1,1) = PGMATRIX_DATA(pParameters,1,1);
	GMATRIX_DATA(M,2,1) = PGMATRIX_DATA(pParameters,2,1);
	GMATRIX_DATA(M,3,1) = PGMATRIX_DATA(pParameters,3,1);
	GMATRIX_DATA(G,1,1) = PGMATRIX_DATA(pParameters,4,1);
	GMATRIX_DATA(G,2,1) = PGMATRIX_DATA(pParameters,5,1);
	GMATRIX_DATA(G,3,1) = PGMATRIX_DATA(pParameters,6,1);

	localization_triad(pY, &IMUMeasure, &MagnetometerMeasure, &M, &G, &q_predicted);
}

int localization_converted_measurement_triad(PGMATRIX pYm, PGMATRIX pPy, PGMATRIX pX_predicted, PGMATRIX pP_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias)
{
	/* Nesse caso, será aplicada a unscented transform para o seguinte problema:
			Y = f(Z)
	   com Z = [X; u_imu_acc; u_mag] e Pz = diag[Px, Pu_imu_acc, Pu_mag];

	*/

	int i;
	UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM_INIT(4,LOCALIZATION_MAXSTATESIZE+6);
	QUATERNIONS_DECLARE(q_predicted);
	GMATRIX_DECLARE(Parameters,6,1);

	QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q0,1);
	QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q1,1);
	QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q2,1);
	QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q3,1);

	// Definição dos tamanhos das matrizes deve ser feito aqui, o que evita alocação dimâmica
	GMATRIX_SETSIZE(XMean ,pX_predicted->Nr+6,1);
	GMATRIX_SETSIZE(XCov ,(pX_predicted->Nr+6) ,(pX_predicted->Nr+6));
	GMATRIX_SETSIZE(_utX ,(pX_predicted->Nr+6) ,1);
	GMATRIX_SETSIZE(_utXSamples	,(pX_predicted->Nr+6),(2*(pX_predicted->Nr+6)+1));
	GMATRIX_SETSIZE(_utYSamples	,(4) ,(2*(pX_predicted->Nr+6)+1));
	GMATRIX_SETSIZE(_utWSamples	,(1) ,(2*(pX_predicted->Nr+6)+1));
	GMATRIX_SETSIZE(_utQSamples	,(pX_predicted->Nr+6), (pX_predicted->Nr+6));

	for(i=1;i<=pX_predicted->Nr;++i){
		GMATRIX_DATA(XMean,i,1) = PGMATRIX_DATA(pX_predicted,i,1);
	}
	GMATRIX_DATA(XMean,pX_predicted->Nr+1,1) = pIMUMeasure->ax;
	GMATRIX_DATA(XMean,pX_predicted->Nr+2,1) = pIMUMeasure->ay;
	GMATRIX_DATA(XMean,pX_predicted->Nr+3,1) = pIMUMeasure->az;
	GMATRIX_DATA(XMean,pX_predicted->Nr+4,1) = pMagnetometerMeasure->mx;
	GMATRIX_DATA(XMean,pX_predicted->Nr+5,1) = pMagnetometerMeasure->my;
	GMATRIX_DATA(XMean,pX_predicted->Nr+6,1) = pMagnetometerMeasure->mz;
	GMATRIX_ZEROES(XCov);
	PGMATRIX_SUBMATRIX_COPY(&XCov,1,1,pP_predicted);
	GMATRIX_DATA(XCov,pX_predicted->Nr+1,pX_predicted->Nr+1) = pIMUMeasure->axvariance;
	GMATRIX_DATA(XCov,pX_predicted->Nr+2,pX_predicted->Nr+2) = pIMUMeasure->ayvariance;
	GMATRIX_DATA(XCov,pX_predicted->Nr+3,pX_predicted->Nr+3) = pIMUMeasure->azvariance;
	GMATRIX_DATA(XCov,pX_predicted->Nr+4,pX_predicted->Nr+4) = pMagnetometerMeasure->mxvariance;
	GMATRIX_DATA(XCov,pX_predicted->Nr+5,pX_predicted->Nr+5) = pMagnetometerMeasure->myvariance;
	GMATRIX_DATA(XCov,pX_predicted->Nr+6,pX_predicted->Nr+6) = pMagnetometerMeasure->mzvariance;

	GMATRIX_DATA(Parameters,1,1) = PGMATRIX_DATA(pM,1,1);
	GMATRIX_DATA(Parameters,2,1) = PGMATRIX_DATA(pM,2,1);
	GMATRIX_DATA(Parameters,3,1) = PGMATRIX_DATA(pM,3,1);
	GMATRIX_DATA(Parameters,4,1) = PGMATRIX_DATA(pG,1,1);
	GMATRIX_DATA(Parameters,5,1) = PGMATRIX_DATA(pG,2,1);
	GMATRIX_DATA(Parameters,6,1) = PGMATRIX_DATA(pG,3,1);

	if (FlagEstimateAccelerometerBias){
		UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM(localization_converted_measurement_triad_UT_biased);
	}
	else{
		UNCERTAINTYPROPAGATION_USCENTEDTRANSFORM(localization_converted_measurement_triad_UT_nonbiased);
	}

	PGMATRIX_COPY(pYm,&YMean);
	PGMATRIX_COPY(pPy,&YCov);

	rotation_quaternionsnormalize(pYm);
	rotation_quaternionscorrectsign(pYm, &q_predicted);

	
	// Retorna 
    return 1; 
}

int localization_converted_measurement_triad_dg_du_imu(PGMATRIX pdg_du_imu, PGMATRIX pX_predicted, PGMATRIX pP_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, int FlagEstimateAccelerometerBias)
{
	IMUMEASURE IMUMeasure;
	QUATERNIONS_DECLARE(q_predicted);
	QUATERNIONS_DECLARE(y);
	QUATERNIONS_DECLARE(y_ax);
	QUATERNIONS_DECLARE(y_ay);
	QUATERNIONS_DECLARE(y_az);
	double delta;


	// Define q_predicted
	QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q0,1);
	QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q1,1);
	QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q2,1);
	QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q3,1);

	// Considera polarização dos acelerômetros
	IMUMeasure.wx = pIMUMeasure->wx;
	IMUMeasure.wy = pIMUMeasure->wy;
	IMUMeasure.wz = pIMUMeasure->wz;
	if (FlagEstimateAccelerometerBias){
		IMUMeasure.ax = pIMUMeasure->ax - PGMATRIX_DATA(pX_predicted,X_bax,1);
		IMUMeasure.ay = pIMUMeasure->ay - PGMATRIX_DATA(pX_predicted,X_bay,1);
		IMUMeasure.az = pIMUMeasure->az - PGMATRIX_DATA(pX_predicted,X_baz,1);
	}
	else{
		IMUMeasure.ax = pIMUMeasure->ax;
		IMUMeasure.ay = pIMUMeasure->ay;
		IMUMeasure.az = pIMUMeasure->az;
	}

	// Método aproximado de calculo de dg_du_imu para o caso da triad
	delta = 1e-5;

	// y = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
	localization_triad(&y, &IMUMeasure, pMagnetometerMeasure, pM, pG, &q_predicted);

    // imumeasure.ax = imumeasure.ax + delta;
    // y_ax = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    // imumeasure.ax = imumeasure.ax - delta;
	IMUMeasure.ax += delta;
	localization_triad(&y_ax, &IMUMeasure, pMagnetometerMeasure, pM, pG, &q_predicted);
	IMUMeasure.ax -= delta;

    // imumeasure.ay = imumeasure.ay + delta;
    // y_ay = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    // imumeasure.ay = imumeasure.ay - delta;
	IMUMeasure.ay += delta;
	localization_triad(&y_ay, &IMUMeasure, pMagnetometerMeasure, pM, pG, &q_predicted);
	IMUMeasure.ay -= delta;

    // imumeasure.az = imumeasure.az + delta;
    // y_az = localization_triad(imumeasure, magnetometermeasure, M, G, q_predicted);
    // imumeasure.az = imumeasure.az - delta;
	IMUMeasure.az += delta;
	localization_triad(&y_az, &IMUMeasure, pMagnetometerMeasure, pM, pG, &q_predicted);
	IMUMeasure.az -= delta;
    
    // dg_du_imu = [(y_ax-y)/delta (y_ay-y)/delta (y_az-y)/delta zeros(4,3)];
	PGMATRIX_SETSIZE(pdg_du_imu,4,6);
	PGMATRIX_ZEROES(pdg_du_imu);
	PGMATRIX_DATA(pdg_du_imu,1,1) = (GMATRIX_DATA(y_ax,1,1) - GMATRIX_DATA(y,1,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,2,1) = (GMATRIX_DATA(y_ax,2,1) - GMATRIX_DATA(y,2,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,3,1) = (GMATRIX_DATA(y_ax,3,1) - GMATRIX_DATA(y,3,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,4,1) = (GMATRIX_DATA(y_ax,4,1) - GMATRIX_DATA(y,4,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,1,2) = (GMATRIX_DATA(y_ay,1,1) - GMATRIX_DATA(y,1,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,2,2) = (GMATRIX_DATA(y_ay,2,1) - GMATRIX_DATA(y,2,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,3,2) = (GMATRIX_DATA(y_ay,3,1) - GMATRIX_DATA(y,3,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,4,2) = (GMATRIX_DATA(y_ay,4,1) - GMATRIX_DATA(y,4,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,1,3) = (GMATRIX_DATA(y_az,1,1) - GMATRIX_DATA(y,1,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,2,3) = (GMATRIX_DATA(y_az,2,1) - GMATRIX_DATA(y,2,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,3,3) = (GMATRIX_DATA(y_az,3,1) - GMATRIX_DATA(y,3,1))/delta;
	PGMATRIX_DATA(pdg_du_imu,4,3) = (GMATRIX_DATA(y_az,4,1) - GMATRIX_DATA(y,4,1))/delta;

}

/*														
int localization_converted_measurement_triad(PGMATRIX pYm, PGMATRIX pPy, PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGMATRIX pX_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG)
{
	int currentcode;

	switch(pFilterStruct->AlgorithmCode){
	case LOCALIZATION_ALGORITHMCODE_EKF2:
		currentcode = pFilterStruct->AlgorithmCode;
		pFilterStruct->AlgorithmCode = LOCALIZATION_ALGORITHMCODE_CEKF;
		if(!localization_converted_measurement_triad(pYm, pPy, pFilterStruct, pX_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG)){
			pFilterStruct->AlgorithmCode = currentcode;
			return 0;
		}
		pFilterStruct->AlgorithmCode = currentcode;
		break;

	case LOCALIZATION_ALGORITHMCODE_CEKF:
		// Ym =  localization_filter_model_gtriad('evaluate',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
		localization_filter_observation_model_triad_evaluate(pYm, pFilterStruct, pX_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG);
		// dg_du_imu = localization_filter_model_gtriad('dg_du_imu',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
		localization_filter_observation_model_triad_dg_du_imu(&dg_du_imu, pFilterStruct, pX_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG);
		// dg_du_mag = localization_filter_model_gtriad('dg_du_mag',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
		localization_filter_observation_model_triad_dg_du_mag(&dg_du_mag, pFilterStruct, pX_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG);
		// dg_dx = localization_filter_model_gtriad('dg_dx',X, Px, imumeasure, magnetometermeasure, M, G, flagestimateaccelerometerbias,1);
		localization_filter_observation_model_triad_dg_du_x(&dg_du_x, pFilterStruct, pX_predicted, pIMUMeasure, pMagnetometerMeasure, pM, pG);
		// P_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
		// P_mag = diag([magnetometermeasure.mxvariance; magnetometermeasure.myvariance; magnetometermeasure.mzvariance]);
		// Py = dg_du_imu * P_imu * dg_du_imu' + dg_du_mag * P_mag * dg_du_mag' + dg_dx * Px * dg_dx';
		break;

	case LOCALIZATION_ALGORITHMCODE_UKF2:
		break;
	default:
	    return 0; 
	}

   // Retorna 
    return 1; 
}

int localization_filter_observation_model_triad_evaluate(PGMATRIX pYm, PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGMATRIX pX_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG)
{
	double ax, ay, az;

	QUATERNIONS_DECLARE(q_predicted);

	QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q0,1);
	QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q1,1);
	QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q2,1);
	QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q3,1);

	if (pFilterStruct->FlagEstimateAccelerometerBias){
		ax = pIMUMeasure->ax;
		ay = pIMUMeasure->ay;
		az = pIMUMeasure->az;
		pIMUMeasure->ax = pIMUMeasure->ax - PGMATRIX_DATA(pX_predicted,X_bax,1);
		pIMUMeasure->ay = pIMUMeasure->ay - PGMATRIX_DATA(pX_predicted,X_bay,1);
		pIMUMeasure->az = pIMUMeasure->az - PGMATRIX_DATA(pX_predicted,X_baz,1);
	}

	localization_triad(pYm, pIMUMeasure, pMagnetometerMeasure, pM, pG, &q_predicted);

	if (pFilterStruct->FlagEstimateAccelerometerBias){
		pIMUMeasure->ax = ax;
		pIMUMeasure->ay = ay;
		pIMUMeasure->az = az;
	}
}

int localization_filter_observation_model_triad_dg_du_imu(PGMATRIX pdg_du_imu, PLOCALIZATIONFILTERSTRUCT pFilterStruct, PGMATRIX pX_predicted, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG)
{

	QUATERNIONS_DECLARE(q_predicted);

	QUATERNIONS_Q0(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q0,1);
	QUATERNIONS_Q1(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q1,1);
	QUATERNIONS_Q2(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q2,1);
	QUATERNIONS_Q3(q_predicted) = PGMATRIX_DATA(pX_predicted,X_q3,1);

	if (pFilterStruct->FlagEstimateAccelerometerBias){
		ax = pIMUMeasure->ax;
		ay = pIMUMeasure->ay;
		az = pIMUMeasure->az;
		pIMUMeasure->ax = pIMUMeasure->ax - PGMATRIX_DATA(pX_predicted,X_bax,1);
		pIMUMeasure->ay = pIMUMeasure->ay - PGMATRIX_DATA(pX_predicted,X_bay,1);
		pIMUMeasure->az = pIMUMeasure->az - PGMATRIX_DATA(pX_predicted,X_baz,1);
	}

	localization_triad(pYm, pIMUMeasure, pMagnetometerMeasure, pM, pG, &q_predicted);

	if (pFilterStruct->FlagEstimateAccelerometerBias){
		pIMUMeasure->ax = ax;
		pIMUMeasure->ay = ay;
		pIMUMeasure->az = az;
	}

}

/*****************************************************************************
*** int localization_filter_prediction(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PIMUMEASURE pIMUMeasure, PGMATRIX pG, double T)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int localization_filter_prediction(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PIMUMEASURE pIMUMeasure, PGMATRIX pG, double T)
{

	// Predição conforme a arquitetura
	switch(pFilterStruct->AlgorithmCode){
	case LOCALIZATION_ALGORITHMCODE_EKF2:
		if(!localization_filter_prediction_ekf2(pFilterStruct, pIMUMeasure, pG, T)) return 0;
		break;

	case LOCALIZATION_ALGORITHMCODE_CEKF:
		if(!localization_filter_prediction_cekf(pFilterStruct, pIMUMeasure, pG, T)) return 0;
		break;

	case LOCALIZATION_ALGORITHMCODE_UKF2:
//		if(!localization_filter_prediction_ukf2(pFilterStruct, pIMUMeasure, pG, T)) return 0;
		break;
	}

	// Reset em caso de matriz P negativa:
/*	if(sum(eig(kf_structure.P)<0)>0)
		% filter reset:
		disp('*** Warning: Kalman filter reset');
		kf_structure.P = kf_structure.Preset;
	end
*/
	// Retorna 
    return 1; 
}                      

int localization_filter_prediction_cekf(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PIMUMEASURE pIMUMeasure, PGMATRIX pG, double T)
{
	GMATRIX_DECLARE(X_previous,LOCALIZATION_MAXSTATESIZE,1);
	GMATRIX_DECLARE(P_previous,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(u_imu,6,1);
	GMATRIX_DECLARE(P_u_imu,6,6);
	GMATRIX_DECLARE(df_dx,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(df_du,LOCALIZATION_MAXSTATESIZE,6);
	GMATRIX_DECLARE(K_previous,LOCALIZATION_MAXSTATESIZE,4);
	GMATRIX_DECLARE(H,4,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(Qtilde,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy1,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(MatDummy2,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);

	/*pFilterStruct->pS_previous_cekf = PGMATRIX_ALLOC(pFilterStruct->Nstates,4);
		PGMATRIX_ZEROES(pFilterStruct->pS_previous_cekf);
		pFilterStruct->pR_previous_cekf = PGMATRIX_ALLOC(4,4);
		PGMATRIX_IDENTITY(pFilterStruct->pR_previous_cekf);
		pFilterStruct->pA_imu_P_imu_cekf = PGMATRIX_ALLOC(pFilterStruct->Nstates,6);
		PGMATRIX_ZEROES(pFilterStruct->pA_imu_P_imu_cekf);
		pFilterStruct->pinnovation_previous_cekf = PGMATRIX_ALLOC(4,1);
		PGMATRIX_ZEROES(pFilterStruct->pinnovation_previous_cekf);*/

	// Definição dos tamanhos das matrizes deve ser feito aqui, o que evita alocação dimâmica
	GMATRIX_SETSIZE(X_previous,pFilterStruct->Nstates,1);
	GMATRIX_SETSIZE(P_previous,pFilterStruct->Nstates,pFilterStruct->Nstates);
	GMATRIX_SETSIZE(df_dx,pFilterStruct->Nstates,pFilterStruct->Nstates);
	GMATRIX_SETSIZE(df_du,pFilterStruct->Nstates,6);
	GMATRIX_SETSIZE(K_previous,pFilterStruct->Nstates,4);
	GMATRIX_SETSIZE(Qtilde,pFilterStruct->Nstates,pFilterStruct->Nstates);

	// X_previous = kf_structure.X;
	PGMATRIX_COPY(&X_previous, pFilterStruct->pX);
	
	// P_previous = kf_structure.P;
	PGMATRIX_COPY(&P_previous, pFilterStruct->pP);
	
	// u_imu = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
	GMATRIX_DATA(u_imu,1,1) = pIMUMeasure->ax;
	GMATRIX_DATA(u_imu,2,1) = pIMUMeasure->ay;
	GMATRIX_DATA(u_imu,3,1) = pIMUMeasure->az;
	GMATRIX_DATA(u_imu,4,1) = pIMUMeasure->wx;
	GMATRIX_DATA(u_imu,5,1) = pIMUMeasure->wy;
	GMATRIX_DATA(u_imu,6,1) = pIMUMeasure->wz;

	// Pu_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
	GMATRIX_ZEROES(P_u_imu);
	GMATRIX_DATA(P_u_imu,1,1) = pIMUMeasure->axvariance;
	GMATRIX_DATA(P_u_imu,2,2) = pIMUMeasure->ayvariance;
	GMATRIX_DATA(P_u_imu,3,3) = pIMUMeasure->azvariance;
	GMATRIX_DATA(P_u_imu,4,4) = pIMUMeasure->wxvariance;
	GMATRIX_DATA(P_u_imu,5,5) = pIMUMeasure->wyvariance;
	GMATRIX_DATA(P_u_imu,6,6) = pIMUMeasure->wzvariance;

	// K_previous = kf_structure.S_previous_cekf * inv(kf_structure.R_previous_cekf);
	PGMATRIX_INVERSE_COPY(&MatDummy1, pFilterStruct->pR_previous_cekf);
	PGMATRIX_MULTIPLY_COPY(&K_previous, pFilterStruct->pS_previous_cekf, &MatDummy1);


	// kf_structure.X = localization_filter_model_f('evaluate',x_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias) + K_previous*kf_structure.innovation_previous_cekf;
	localization_filter_process_model_evaluate(pFilterStruct->pX, &X_previous, &u_imu, pG, T, pFilterStruct->FlagEstimateAccelerometerBias);
	PGMATRIX_MULTIPLY_ADD(pFilterStruct->pX,&K_previous,pFilterStruct->pinnovation_previous_cekf);

	// df_dx = localization_filter_model_f('df_dx',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
	localization_filter_process_model_df_dx(&df_dx, &X_previous, &u_imu, pG, T, pFilterStruct->FlagEstimateAccelerometerBias);
	
	// df_du = localization_filter_model_f('df_du',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
	localization_filter_process_model_df_du(&df_du, &X_previous, &u_imu, pG, T, pFilterStruct->FlagEstimateAccelerometerBias);

	// Qtilde = kf_structure.Q_cekf + df_du_imu*P_u_imu*df_du_imu'; 
	PGMATRIX_TRIPLEMULTIPLY_COPY_EXTENDED(&Qtilde, &df_du, 0, &P_u_imu, 0, &df_du, 1, &MatDummy1);
	PGMATRIX_ADD(&Qtilde, pFilterStruct->pQ);

	// if kf_structure.flagestimateaccelerometerbias
	//     H = [eye(4) zeros(4,9)];
	// else
	//     H = [eye(4) zeros(4,6)];
	// end
	if (pFilterStruct->FlagEstimateAccelerometerBias){
		GMATRIX_SETSIZE(H,4,13);
		GMATRIX_ZEROES(H);
	}
	else{
		GMATRIX_SETSIZE(H,4,10);
		GMATRIX_ZEROES(H);
	}
	GMATRIX_DATA(H,1,1) = 1.0;
	GMATRIX_DATA(H,2,2) = 1.0;
	GMATRIX_DATA(H,3,3) = 1.0;
	GMATRIX_DATA(H,4,4) = 1.0;

	// kf_structure.P = (df_dx-K_previous*H)*P_previous*(df_dx-K_previous*H)' + Qtilde - K_previous*kf_structure.R_previous_cekf*K_previous';
	PGMATRIX_TRIPLEMULTIPLY_COPY_EXTENDED(pFilterStruct->pP, &K_previous, 0, pFilterStruct->pR_previous_cekf, 0, &K_previous, 1, &MatDummy1);
	PGMATRIX_MULTIPLY_CONST(pFilterStruct->pP, -1.0); // kf_structure.P = - K_previous*kf_structure.R_previous_cekf*K_previous';
	PGMATRIX_ADD(pFilterStruct->pP, &Qtilde); // kf_structure.P = Qtilde - K_previous*kf_structure.R_previous_cekf*K_previous';
	PGMATRIX_MULTIPLY_COPY(&MatDummy1, &K_previous, &H); 
	PGMATRIX_SUBTRACT_COPY(&MatDummy2, &df_dx, &MatDummy1); // MatDummy2 = df_dx-K_previous*H
	PGMATRIX_TRIPLEMULTIPLY_ADD_EXTENDED(pFilterStruct->pP, &MatDummy2, 0, &P_previous, 0, &MatDummy2, 1, &MatDummy1);

	// kf_structure.A_imu_P_imu_cekf = df_du_imu*P_u_imu;
	PGMATRIX_MULTIPLY_COPY(pFilterStruct->pA_imu_P_imu_cekf, &df_du, &P_u_imu); 
	
	// Retorna 
    return 1; 
}  

int localization_filter_prediction_ekf2(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PIMUMEASURE pIMUMeasure, PGMATRIX pG, double T)
{
	GMATRIX_DECLARE(X_previous,LOCALIZATION_MAXSTATESIZE,1);
	GMATRIX_DECLARE(P_previous,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(u_imu,6,1);
	GMATRIX_DECLARE(P_u_imu,6,6);
	GMATRIX_DECLARE(df_dx,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);
	GMATRIX_DECLARE(df_du,LOCALIZATION_MAXSTATESIZE,6);
	GMATRIX_DECLARE(MatDummy,LOCALIZATION_MAXSTATESIZE,LOCALIZATION_MAXSTATESIZE);

	// Definição dos tamanhos das matrizes deve ser feito aqui, o que evita alocação dimâmica
	GMATRIX_SETSIZE(X_previous,pFilterStruct->Nstates,1);
	GMATRIX_SETSIZE(P_previous,pFilterStruct->Nstates,pFilterStruct->Nstates);
	GMATRIX_SETSIZE(df_dx,pFilterStruct->Nstates,pFilterStruct->Nstates);
	GMATRIX_SETSIZE(df_du,pFilterStruct->Nstates,6);

	// X_previous = kf_structure.X;
	PGMATRIX_COPY(&X_previous, pFilterStruct->pX);
	
	// P_previous = kf_structure.P;
	PGMATRIX_COPY(&P_previous, pFilterStruct->pP);
	
	// u_imu = [imumeasure.ax; imumeasure.ay; imumeasure.az; imumeasure.wx; imumeasure.wy; imumeasure.wz];
	GMATRIX_DATA(u_imu,1,1) = pIMUMeasure->ax;
	GMATRIX_DATA(u_imu,2,1) = pIMUMeasure->ay;
	GMATRIX_DATA(u_imu,3,1) = pIMUMeasure->az;
	GMATRIX_DATA(u_imu,4,1) = pIMUMeasure->wx;
	GMATRIX_DATA(u_imu,5,1) = pIMUMeasure->wy;
	GMATRIX_DATA(u_imu,6,1) = pIMUMeasure->wz;

	// Pu_imu = diag([imumeasure.axvariance; imumeasure.ayvariance; imumeasure.azvariance; imumeasure.wxvariance; imumeasure.wyvariance; imumeasure.wzvariance]);
	GMATRIX_ZEROES(P_u_imu);
	GMATRIX_DATA(P_u_imu,1,1) = pIMUMeasure->axvariance;
	GMATRIX_DATA(P_u_imu,2,2) = pIMUMeasure->ayvariance;
	GMATRIX_DATA(P_u_imu,3,3) = pIMUMeasure->azvariance;
	GMATRIX_DATA(P_u_imu,4,4) = pIMUMeasure->wxvariance;
	GMATRIX_DATA(P_u_imu,5,5) = pIMUMeasure->wyvariance;
	GMATRIX_DATA(P_u_imu,6,6) = pIMUMeasure->wzvariance;

	// kf_structure.X = localization_filter_model_f('evaluate',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
	localization_filter_process_model_evaluate(pFilterStruct->pX, &X_previous, &u_imu, pG, T, pFilterStruct->FlagEstimateAccelerometerBias);
	
	// df_dx = localization_filter_model_f('df_dx',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
	localization_filter_process_model_df_dx(&df_dx, &X_previous, &u_imu, pG, T, pFilterStruct->FlagEstimateAccelerometerBias);
	
	// df_du = localization_filter_model_f('df_du',X_previous,u_imu,G,T,kf_structure.flagestimateaccelerometerbias);
	localization_filter_process_model_df_du(&df_du, &X_previous, &u_imu, pG, T, pFilterStruct->FlagEstimateAccelerometerBias);
	
	// kf_structure.P = df_dx*P_previous*df_dx' + df_du*P_u_imu*df_du' + kf_structure.Q_ekf;
	PGMATRIX_TRIPLEMULTIPLY_COPY_EXTENDED(pFilterStruct->pP, &df_dx, 0, &P_previous, 0, &df_dx, 1, &MatDummy);
	PGMATRIX_TRIPLEMULTIPLY_ADD_EXTENDED(pFilterStruct->pP, &df_du, 0, &P_u_imu, 0, &df_du, 1, &MatDummy);
	PGMATRIX_ADD(pFilterStruct->pP, pFilterStruct->pQ);

	// Retorna 
    return 1; 
}                      

int localization_filter_process_model_df_du(PGMATRIX pdf_du, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias)
{
	double sx, sy, sz, a, ca2, sa2;
	double x1, x2, x3, x4;
	double da_dwx, da_dwy, da_dwz;
	GMATRIX_DECLARE(acc_b,3,1);
	GMATRIX_DECLARE(acc_e,3,1);
	DCM_DECLARE(R);
	QUATERNIONS_DECLARE(q);
	GMATRIX_DECLARE(dR_dxi,3,3);
	GMATRIX_DECLARE(da_dx,3,4);

	sx = PGMATRIX_DATA(pu_imu,U_wx,1)*T;
	sy = PGMATRIX_DATA(pu_imu,U_wy,1)*T;
	sz = PGMATRIX_DATA(pu_imu,U_wz,1)*T;
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
		da_dwx = (0.5)*(1.0/a)*(2.0*sx*T);
		da_dwy = (0.5)*(1.0/a)*(2.0*sy*T);
		da_dwz = (0.5)*(1.0/a)*(2.0*sz*T);

		x1 = PGMATRIX_DATA(pX_previous,X_q0,1);
		x2 = PGMATRIX_DATA(pX_previous,X_q1,1);
		x3 = PGMATRIX_DATA(pX_previous,X_q2,1);
		x4 = PGMATRIX_DATA(pX_previous,X_q3,1);

		// df_du(1,3) = -sa2*(1/2)*da_dwx*x1 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x2 ); % df_dwx
		// df_du(1,4) = -sa2*(1/2)*da_dwy*x1 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x3 ); % df_dwy
		// df_du(1,5) = -sa2*(1/2)*da_dwz*x1 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x4 ); % df_dwz

		PGMATRIX_DATA(pdf_du,1,3) = -sa2*(0.5)*da_dwx*x1 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwx * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x2 ); 
		PGMATRIX_DATA(pdf_du,1,4) = -sa2*(0.5)*da_dwy*x1 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwy * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x3 ); 
		PGMATRIX_DATA(pdf_du,1,5) = -sa2*(0.5)*da_dwz*x1 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwz * ( sx*x2 + sy*x3 + sz*x4) - sa2/a*( T*x4 ); 

		// df_du(2,3) = -sa2*(1/2)*da_dwx*x2 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*(-T*x1 ); % df_dwx
		// df_du(2,4) = -sa2*(1/2)*da_dwy*x2 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*( T*x4 ); % df_dwy
		// df_du(2,5) = -sa2*(1/2)*da_dwz*x2 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*(-T*x3 ); % df_dwz

		PGMATRIX_DATA(pdf_du,2,3) = -sa2*(0.5)*da_dwx*x2 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwx * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*(-T*x1 );
		PGMATRIX_DATA(pdf_du,2,4) = -sa2*(0.5)*da_dwy*x2 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwy * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*( T*x4 );
		PGMATRIX_DATA(pdf_du,2,5) = -sa2*(0.5)*da_dwz*x2 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwz * (-sx*x1 - sz*x3 + sy*x4) - sa2/a*(-T*x3 );

		// df_du(3,3) = -sa2*(1/2)*da_dwx*x3 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*(-T*x4 ); % df_dwx
		// df_du(3,4) = -sa2*(1/2)*da_dwy*x3 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*(-T*x1 ); % df_dwy
		// df_du(3,5) = -sa2*(1/2)*da_dwz*x3 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*( T*x2 ); % df_dwz

		PGMATRIX_DATA(pdf_du,3,3) = -sa2*(0.5)*da_dwx*x3 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwx * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*(-T*x4 );
		PGMATRIX_DATA(pdf_du,3,4) = -sa2*(0.5)*da_dwy*x3 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwy * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*(-T*x1 );
		PGMATRIX_DATA(pdf_du,3,5) = -sa2*(0.5)*da_dwz*x3 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwz * (-sy*x1 + sz*x2 - sx*x4) - sa2/a*( T*x2 );

		// df_du(4,3) = -sa2*(1/2)*da_dwx*x4 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwx * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*( T*x3 ); % df_dwx
		// df_du(4,4) = -sa2*(1/2)*da_dwy*x4 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwy * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*(-T*x2 ); % df_dwy
		// df_du(4,5) = -sa2*(1/2)*da_dwz*x4 + (-ca2*(1/2)*a + sa2)/(a*a)* da_dwz * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*(-T*x1 ); % df_dwz

		PGMATRIX_DATA(pdf_du,4,3) = -sa2*(0.5)*da_dwx*x4 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwx * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*( T*x3 );
		PGMATRIX_DATA(pdf_du,4,4) = -sa2*(0.5)*da_dwy*x4 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwy * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*(-T*x2 );
		PGMATRIX_DATA(pdf_du,4,5) = -sa2*(0.5)*da_dwz*x4 + (-ca2*(0.5)*a + sa2)/(a*a)* da_dwz * (-sz*x1 - sy*x2 + sx*x3) - sa2/a*(-T*x1 );

	}

    // R = quaternions2dcm(x(1:4));
	QUATERNIONS_Q0(q) = PGMATRIX_DATA(pX_previous,X_q0,1);
	QUATERNIONS_Q1(q) = PGMATRIX_DATA(pX_previous,X_q1,1);
	QUATERNIONS_Q2(q) = PGMATRIX_DATA(pX_previous,X_q2,1);
	QUATERNIONS_Q3(q) = PGMATRIX_DATA(pX_previous,X_q3,1);
	rotation_quaternions2dcm(&q, &R);

	// df_du(5:7 ,1:3)  = R*(T^2)/2;
	PGMATRIX_SUBMATRIX_COPY(pdf_du, 5, 1, &R);
	PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_du,5,7,1,3,(-0.5*T*T));

	// df_du(8:10,1:3)  = R*(T);
	PGMATRIX_SUBMATRIX_COPY(pdf_du, 8, 1, &R);
	PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_du,8,10,1,3,(T));

	// Retorna 
    return 1; 
} 

int localization_filter_process_model_df_dx(PGMATRIX pdf_dx, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias)
{
	double sx, sy, sz, a, ca2, sa2;
	double x1, x2, x3, x4;
	DCM_DECLARE(R);
	QUATERNIONS_DECLARE(q);
	GMATRIX_DECLARE(dR_dxi,3,3);
	GMATRIX_DECLARE(da_dx,3,4);

	sx = PGMATRIX_DATA(pu_imu,U_wx,1)*T;
	sy = PGMATRIX_DATA(pu_imu,U_wy,1)*T;
	sz = PGMATRIX_DATA(pu_imu,U_wz,1)*T;
	a = sqrt(sx*sx + sy*sy + sz*sz);

	ca2 = cos(a/2.0);
	sa2 = sin(a/2.0);

	// zera df_dx:
	PGMATRIX_ZEROES(pdf_dx);

	// %Atualiza a atitude
	if (a>0.0){
        PGMATRIX_DATA(pdf_dx,1,1) = ca2;
        PGMATRIX_DATA(pdf_dx,1,2) = -(sa2/a)*sx;
        PGMATRIX_DATA(pdf_dx,1,3) = -(sa2/a)*sy;
        PGMATRIX_DATA(pdf_dx,1,4) = -(sa2/a)*sz;

        PGMATRIX_DATA(pdf_dx,2,1) = (sa2/a)*sx;
        PGMATRIX_DATA(pdf_dx,2,2) = ca2;
        PGMATRIX_DATA(pdf_dx,2,3) = (sa2/a)*sz;
        PGMATRIX_DATA(pdf_dx,2,4) = -(sa2/a)*sy;

        PGMATRIX_DATA(pdf_dx,3,1) = -(sa2/a)*sy;
        PGMATRIX_DATA(pdf_dx,3,2) = (sa2/a)*sz;
        PGMATRIX_DATA(pdf_dx,3,3) = ca2;
        PGMATRIX_DATA(pdf_dx,3,4) = -(sa2/a)*sx;

        PGMATRIX_DATA(pdf_dx,4,1) = -(sa2/a)*sz;
        PGMATRIX_DATA(pdf_dx,4,2) = -(sa2/a)*sy;
        PGMATRIX_DATA(pdf_dx,4,3) = (sa2/a)*sx;
        PGMATRIX_DATA(pdf_dx,4,4) = ca2;
	}
	else{
        PGMATRIX_DATA(pdf_dx,1,1) = 1.0;
        PGMATRIX_DATA(pdf_dx,2,2) = 1.0;
        PGMATRIX_DATA(pdf_dx,3,3) = 1.0;
        PGMATRIX_DATA(pdf_dx,4,4) = 1.0;
	}

    // R = quaternions2dcm(x(1:4));
	QUATERNIONS_Q0(q) = PGMATRIX_DATA(pX_previous,X_q0,1);
	QUATERNIONS_Q1(q) = PGMATRIX_DATA(pX_previous,X_q1,1);
	QUATERNIONS_Q2(q) = PGMATRIX_DATA(pX_previous,X_q2,1);
	QUATERNIONS_Q3(q) = PGMATRIX_DATA(pX_previous,X_q3,1);
	rotation_quaternions2dcm(&q, &R);

	x1 = PGMATRIX_DATA(pX_previous,X_q0,1);
	x2 = PGMATRIX_DATA(pX_previous,X_q1,1);
	x3 = PGMATRIX_DATA(pX_previous,X_q2,1);
	x4 = PGMATRIX_DATA(pX_previous,X_q3,1);


	// dR_dx1 = 2*[(2*x(1))/2    (-x(4))       (x(3));
	// 			( x(4))       (2*x(1))/2    (-x(2));
	// 			(-x(3))       (x(2))        (2*x(1))/2];
	GMATRIX_DATA(dR_dxi,1,1) = (2*x1);    
	GMATRIX_DATA(dR_dxi,1,2) = 2.0*(-x4);      
	GMATRIX_DATA(dR_dxi,1,3) = 2.0*(x3);
	GMATRIX_DATA(dR_dxi,2,1) = 2.0*( x4);       
	GMATRIX_DATA(dR_dxi,2,2) = (2*x1);    
	GMATRIX_DATA(dR_dxi,2,3) = 2.0*(-x2);
	GMATRIX_DATA(dR_dxi,3,1) = 2.0*(-x3);     
	GMATRIX_DATA(dR_dxi,3,2) = 2.0*(x2);      
	GMATRIX_DATA(dR_dxi,3,3) = (2*x1);
	// da_dx(:,1) = dR_dx1*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
	GMATRIX_DATA(da_dx,1,1) = GMATRIX_DATA(dR_dxi,1,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,1,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,1,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,2,1) = GMATRIX_DATA(dR_dxi,2,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,2,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,2,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,3,1) = GMATRIX_DATA(dR_dxi,3,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,3,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,3,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	
	// dR_dx2 = 2*[(2*x(2))/2    (x(3))        (x(4));
	// 			(x(3))        (-2*x(2))/2   (-x(1));
	// 			(x(4))        (x(1))        (-2*x(2))/2];
	GMATRIX_DATA(dR_dxi,1,1) = (2*x2);
	GMATRIX_DATA(dR_dxi,1,2) = 2.0*(x3);   
	GMATRIX_DATA(dR_dxi,1,3) = 2.0*(x4);
	GMATRIX_DATA(dR_dxi,2,1) = 2.0*(x3);   
	GMATRIX_DATA(dR_dxi,2,2) = (-2*x2);   
	GMATRIX_DATA(dR_dxi,2,3) = 2.0*(-x1);
	GMATRIX_DATA(dR_dxi,3,1) = 2.0*(x4);        
	GMATRIX_DATA(dR_dxi,3,2) = 2.0*(x1);        
	GMATRIX_DATA(dR_dxi,3,3) = (-2*x2);

	// da_dx(:,2) = dR_dx2*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
	GMATRIX_DATA(da_dx,1,2) = GMATRIX_DATA(dR_dxi,1,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,1,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,1,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,2,2) = GMATRIX_DATA(dR_dxi,2,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,2,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,2,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,3,2) = GMATRIX_DATA(dR_dxi,3,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,3,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,3,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	
	// dR_dx3 = 2*[(-2*x(3))/2   (x(2))        (x(1));
	// 			(x(2))        (2*x(3))/2    (x(4));
	// 			(-x(1))       (x(4))        (-2*x(3))/2];
	GMATRIX_DATA(dR_dxi,1,1) = (-2*x3);
	GMATRIX_DATA(dR_dxi,1,2) = 2.0*(x2);
	GMATRIX_DATA(dR_dxi,1,3) = 2.0*(x1);
	GMATRIX_DATA(dR_dxi,2,1) = 2.0*(x2);  
	GMATRIX_DATA(dR_dxi,2,2) = (2*x3);   
	GMATRIX_DATA(dR_dxi,2,3) = 2.0*(x4);
	GMATRIX_DATA(dR_dxi,3,1) = 2.0*(-x1); 
	GMATRIX_DATA(dR_dxi,3,2) = 2.0*(x4);
	GMATRIX_DATA(dR_dxi,3,3) = (-2*x3);

	// da_dx(:,3) = dR_dx3*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
	GMATRIX_DATA(da_dx,1,3) = GMATRIX_DATA(dR_dxi,1,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,1,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,1,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,2,3) = GMATRIX_DATA(dR_dxi,2,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,2,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,2,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,3,3) = GMATRIX_DATA(dR_dxi,3,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,3,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,3,3)*PGMATRIX_DATA(pu_imu,U_az,1);

	// dR_dx4 = 2*[(-2*x(4))/2   -x(1)         (x(2));
	// 			(x(1))        (-2*x(4))/2    (x(3));
	// 			(x(2))        (x(3))         (2*x(4))/2];
	GMATRIX_DATA(dR_dxi,1,1) = (-2*x4);
	GMATRIX_DATA(dR_dxi,1,2) = 2.0*(-x1);
	GMATRIX_DATA(dR_dxi,1,3) = 2.0*(x2);
	GMATRIX_DATA(dR_dxi,2,1) = 2.0*(x1);      
	GMATRIX_DATA(dR_dxi,2,2) = (-2*x4);  
	GMATRIX_DATA(dR_dxi,2,3) = 2.0*(x3);
	GMATRIX_DATA(dR_dxi,3,1) = 2.0*(x2);     
	GMATRIX_DATA(dR_dxi,3,2) = 2.0*(x3);       
	GMATRIX_DATA(dR_dxi,3,3) = (2*x4);

	// da_dx(:,4) = dR_dx4*[imumeasure.ax; imumeasure.ay; imumeasure.az]; % somente quaternions
	GMATRIX_DATA(da_dx,1,4) = GMATRIX_DATA(dR_dxi,1,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,1,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,1,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,2,4) = GMATRIX_DATA(dR_dxi,2,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,2,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,2,3)*PGMATRIX_DATA(pu_imu,U_az,1);
	GMATRIX_DATA(da_dx,3,4) = GMATRIX_DATA(dR_dxi,3,1)*PGMATRIX_DATA(pu_imu,U_ax,1) + GMATRIX_DATA(dR_dxi,3,2)*PGMATRIX_DATA(pu_imu,U_ay,1) + GMATRIX_DATA(dR_dxi,3,3)*PGMATRIX_DATA(pu_imu,U_az,1);

	// df_dx(5:7,1:4)  = da_dx*(T^2)/2;
	PGMATRIX_SUBMATRIX_COPY(pdf_dx, 5, 1, &da_dx);
	PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx,5,7,1,4,(0.5*T*T));

	// df_dx(5:7,5:7)  = eye(3);
    PGMATRIX_DATA(pdf_dx,5,5) = 1.0;
    PGMATRIX_DATA(pdf_dx,6,6) = 1.0;
    PGMATRIX_DATA(pdf_dx,7,7) = 1.0;
	
	// df_dx(5:7,8:10) = eye(3)*T;
    PGMATRIX_DATA(pdf_dx,5,8) = T;
    PGMATRIX_DATA(pdf_dx,6,9) = T;
    PGMATRIX_DATA(pdf_dx,7,10) = T;
	
	// if flagestimateaccelerometerbias
	//     df_dx(5:7,11:13) = -R*(T^2)/2;
	// end
	if (FlagEstimateAccelerometerBias){
		PGMATRIX_SUBMATRIX_COPY(pdf_dx, 5, 11, &R);
		PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx,5,7,11,13,(-0.5*T*T));
	}

	// df_dx(8:10,1:4)  = da_dx*(T);
	PGMATRIX_SUBMATRIX_COPY(pdf_dx, 8, 1, &da_dx);
	PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx,8,10,1,4,(T));

	// df_dx(8:10,8:10) = eye(3);
    PGMATRIX_DATA(pdf_dx,8,8) = 1.0;
    PGMATRIX_DATA(pdf_dx,9,9) = 1.0;
    PGMATRIX_DATA(pdf_dx,10,10) = 1.0;

	// if flagestimateaccelerometerbias
	//     df_dx(8:10,11:13) = -R*(T);
	// end
	if (FlagEstimateAccelerometerBias){
		PGMATRIX_SUBMATRIX_COPY(pdf_dx, 8, 11, &R);
		PGMATRIX_SUBMATRIX_MULTIPLY_CONST(pdf_dx,8,10,11,13,(-T));
	}

	// if flagestimateaccelerometerbias
	//     df_dx(11:13,11:13) = eye(3);
	// end
	if (FlagEstimateAccelerometerBias){
		PGMATRIX_DATA(pdf_dx,11,11) = 1.0;
		PGMATRIX_DATA(pdf_dx,12,12) = 1.0;
		PGMATRIX_DATA(pdf_dx,13,13) = 1.0;
	}

	// Retorna 
    return 1; 
} 

int localization_filter_process_model_evaluate(PGMATRIX pX, PGMATRIX pX_previous, PGMATRIX pu_imu, PGMATRIX pG, double T, int FlagEstimateAccelerometerBias)
{
	double sx, sy, sz, a, ca2, sa2;
	double x1, x2, x3, x4;
	GMATRIX_DECLARE(acc_b,3,1);
	GMATRIX_DECLARE(acc_e,3,1);
	DCM_DECLARE(R);
	QUATERNIONS_DECLARE(q);

	sx = PGMATRIX_DATA(pu_imu,U_wx,1)*T;
	sy = PGMATRIX_DATA(pu_imu,U_wy,1)*T;
	sz = PGMATRIX_DATA(pu_imu,U_wz,1)*T;
	a = sqrt(sx*sx + sy*sy + sz*sz);

	ca2 = cos(a/2.0);
	sa2 = sin(a/2.0);

	x1 = PGMATRIX_DATA(pX_previous,X_q0,1);
	x2 = PGMATRIX_DATA(pX_previous,X_q1,1);
	x3 = PGMATRIX_DATA(pX_previous,X_q2,1);
	x4 = PGMATRIX_DATA(pX_previous,X_q3,1);

	// %Atualiza a atitude
	if (a>0.0){
		// q(1) = cos(a/2)*x(1) - sin(a/2)/a*( sx*x(2) + sy*x(3) + sz*x(4));
		PGMATRIX_DATA(pX,X_q0,1) = ca2*x1 - (sa2/a)*(sx*x2 + sy*x3 + sz*x4);
		// q(2) = cos(a/2)*x(2) - sin(a/2)/a*(-sx*x(1) - sz*x(3) + sy*x(4));
		PGMATRIX_DATA(pX,X_q1,1) = ca2*x2 - (sa2/a)*(-sx*x1 - sz*x3 + sy*x4);
		// q(3) = cos(a/2)*x(3) - sin(a/2)/a*(-sy*x(1) + sz*x(2) - sx*x(4));
		PGMATRIX_DATA(pX,X_q2,1) = ca2*x3 - (sa2/a)*(-sy*x1 + sz*x2 - sx*x4);
		// q(4) = cos(a/2)*x(4) - sin(a/2)/a*(-sz*x(1) - sy*x(2) + sx*x(3));
		PGMATRIX_DATA(pX,X_q3,1) = ca2*x4 - (sa2/a)*(-sz*x1 - sy*x2 + sx*x3);
	}
	else{
		PGMATRIX_DATA(pX,X_q0,1) = x1;
		PGMATRIX_DATA(pX,X_q1,1) = x2;
		PGMATRIX_DATA(pX,X_q2,1) = x3;
		PGMATRIX_DATA(pX,X_q3,1) = x4;
	}
	// %preserva a atitude, posição e velocidade prévia
	// if flagestimateaccelerometerbias
	//     a = quaternions2dcm(x(1:4))*[imumeasure.ax-x(11); imumeasure.ay-x(12); imumeasure.az-x(13)] + G; % gravity compensation
	// else
	//     a = quaternions2dcm(x(1:4))*[imumeasure.ax; imumeasure.ay; imumeasure.az] + G; % gravity compensation
	// end
	QUATERNIONS_Q0(q) = PGMATRIX_DATA(pX_previous,X_q0,1);
	QUATERNIONS_Q1(q) = PGMATRIX_DATA(pX_previous,X_q1,1);
	QUATERNIONS_Q2(q) = PGMATRIX_DATA(pX_previous,X_q2,1);
	QUATERNIONS_Q3(q) = PGMATRIX_DATA(pX_previous,X_q3,1);
	rotation_quaternions2dcm(&q, &R);

	if (FlagEstimateAccelerometerBias){
		GMATRIX_DATA(acc_b,1,1) = PGMATRIX_DATA(pu_imu,U_ax,1) - PGMATRIX_DATA(pX_previous,X_bax,1);
		GMATRIX_DATA(acc_b,2,1) = PGMATRIX_DATA(pu_imu,U_ay,1) - PGMATRIX_DATA(pX_previous,X_bay,1);
		GMATRIX_DATA(acc_b,3,1) = PGMATRIX_DATA(pu_imu,U_az,1) - PGMATRIX_DATA(pX_previous,X_baz,1);
	}
	else{
		GMATRIX_DATA(acc_b,1,1) = PGMATRIX_DATA(pu_imu,U_ax,1);
		GMATRIX_DATA(acc_b,2,1) = PGMATRIX_DATA(pu_imu,U_ay,1);
		GMATRIX_DATA(acc_b,3,1) = PGMATRIX_DATA(pu_imu,U_az,1);
	}
	GMATRIX_DATA(acc_e,1,1) = GMATRIX_DATA(R,1,1)*GMATRIX_DATA(acc_b,1,1) + GMATRIX_DATA(R,1,2)*GMATRIX_DATA(acc_b,2,1) + GMATRIX_DATA(R,1,3)*GMATRIX_DATA(acc_b,3,1) + PGMATRIX_DATA(pG,1,1); 
	GMATRIX_DATA(acc_e,2,1) = GMATRIX_DATA(R,2,1)*GMATRIX_DATA(acc_b,1,1) + GMATRIX_DATA(R,2,2)*GMATRIX_DATA(acc_b,2,1) + GMATRIX_DATA(R,2,3)*GMATRIX_DATA(acc_b,3,1) + PGMATRIX_DATA(pG,2,1); 
	GMATRIX_DATA(acc_e,3,1) = GMATRIX_DATA(R,3,1)*GMATRIX_DATA(acc_b,1,1) + GMATRIX_DATA(R,3,2)*GMATRIX_DATA(acc_b,2,1) + GMATRIX_DATA(R,3,3)*GMATRIX_DATA(acc_b,3,1) + PGMATRIX_DATA(pG,3,1); 

	// p = x(5:7) + x(8:10)*T + a*(T^2)/2;
	PGMATRIX_DATA(pX,X_px,1) = PGMATRIX_DATA(pX_previous,X_px,1) + PGMATRIX_DATA(pX_previous,X_vx,1)*T + GMATRIX_DATA(acc_e,1,1)*(T*T)/2.0;
	PGMATRIX_DATA(pX,X_py,1) = PGMATRIX_DATA(pX_previous,X_py,1) + PGMATRIX_DATA(pX_previous,X_vy,1)*T + GMATRIX_DATA(acc_e,2,1)*(T*T)/2.0;
	PGMATRIX_DATA(pX,X_pz,1) = PGMATRIX_DATA(pX_previous,X_pz,1) + PGMATRIX_DATA(pX_previous,X_vz,1)*T + GMATRIX_DATA(acc_e,3,1)*(T*T)/2.0;

	// v = x(8:10) + T*a;
	PGMATRIX_DATA(pX,X_vx,1) = PGMATRIX_DATA(pX_previous,X_vx,1) + GMATRIX_DATA(acc_e,1,1)*T;
	PGMATRIX_DATA(pX,X_vy,1) = PGMATRIX_DATA(pX_previous,X_vy,1) + GMATRIX_DATA(acc_e,2,1)*T;
	PGMATRIX_DATA(pX,X_vz,1) = PGMATRIX_DATA(pX_previous,X_vz,1) + GMATRIX_DATA(acc_e,3,1)*T;

	// % parâmetros:
	// if flagestimateaccelerometerbias
	//     var_out(11:13,1) = x(11:13,1);
	// end
	if (FlagEstimateAccelerometerBias){
		PGMATRIX_DATA(pX,X_bax,1) = PGMATRIX_DATA(pX_previous,X_bax,1);
		PGMATRIX_DATA(pX,X_bay,1) = PGMATRIX_DATA(pX_previous,X_bay,1);
		PGMATRIX_DATA(pX,X_baz,1) = PGMATRIX_DATA(pX_previous,X_baz,1);
	}

	// Retorna 
    return 1; 
} 
/*****************************************************************************
*** int localization_triad(PQUATERNIONS pq, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, PQUATERNIONS pq_previous)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int localization_triad(PQUATERNIONS pq, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, PQUATERNIONS pq_previous)
{
	double aux;
	GMATRIX_DECLARE(mag,3,1); // valores locais normalizados, para não alterar na origem
	GMATRIX_DECLARE(acc,3,1); // valores locais normalizados, para não alterar na origem
	GMATRIX_DECLARE(M,3,1); // valores locais normalizados, para não alterar na origem
	GMATRIX_DECLARE(G,3,1); // valores locais normalizados, para não alterar na origem

	GMATRIX_DECLARE(MatDummy,3,1); 

	GMATRIX_DECLARE(I_b,3,1); 
	GMATRIX_DECLARE(J_b,3,1); 
	GMATRIX_DECLARE(K_b,3,1); 
	GMATRIX_DECLARE(I_i,3,1); 
	GMATRIX_DECLARE(J_i,3,1); 
	GMATRIX_DECLARE(K_i,3,1); 

	GMATRIX_DECLARE(R_i,3,3); 
	GMATRIX_DECLARE(R_b,3,3); 
	GMATRIX_DECLARE(R_b_i,3,3); 

	// testar se medicoes do magnetometro são validas:
	if(!pMagnetometerMeasure->FlagValidMeasure) return 0;

	// mag = [magnetometermeasure.mx; magnetometermeasure.my; magnetometermeasure.mz];
	GMATRIX_DATA(mag,1,1) = pMagnetometerMeasure->mx;
	GMATRIX_DATA(mag,2,1) = pMagnetometerMeasure->my;
	GMATRIX_DATA(mag,3,1) = pMagnetometerMeasure->mz;
	// mag = (mag/sqrt(mag(1)^2 + mag(2)^2 + mag(3)^2));
	aux = PGMATRIX_NORM(&mag);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&mag,(1.0/aux));

	// %armazena a aceleração (força específica). A acerelação é negada, pois leva-se em conta que os
	// %acelerometros medem a força específica sob sua massa e não a aceleração
	// acc = zeros(3,1);
	// acc(1) = -imumeasure.ax;
	// acc(2) = -imumeasure.ay;
	// acc(3) = -imumeasure.az;
	GMATRIX_DATA(acc,1,1) = -pIMUMeasure->ax;
	GMATRIX_DATA(acc,2,1) = -pIMUMeasure->ay;
	GMATRIX_DATA(acc,3,1) = -pIMUMeasure->az;
	// acc = acc/sqrt(acc(1)^2 + acc(2)^2 + acc(3)^2);
	aux = PGMATRIX_NORM(&acc);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&acc,(1.0/aux));

	// G = G/sqrt(G(1)^2 + G(2)^2 + G(3)^2);
	PGMATRIX_COPY(&G,pG);
	aux = PGMATRIX_NORM(&G);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&G,(1.0/aux));

	// M = M/sqrt(M(1)^2 + M(2)^2 + M(3)^2);
	PGMATRIX_COPY(&M,pM);
	aux = PGMATRIX_NORM(&M);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&M,(1.0/aux));


	// Inicializa as variaveis do corpo. sendo "acc" as medidas do acelerometro e
	// "mag" as medidas do magnetometro. Lembrando que a função cross(,)
	// representa produto vetorial entre dois vetores.
	// aux = acc + mag;
	GMATRIX_ADD_COPY(I_b,acc,mag);
	// I_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
	aux = PGMATRIX_NORM(&I_b);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&I_b,(1.0/aux));


	// aux = cross(I_b,(acc - mag));
	GMATRIX_SUBTRACT_COPY(MatDummy,acc,mag);
	GMATRIX_CROSS_COPY(J_b,I_b,MatDummy);
	
	// J_b = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
	aux = PGMATRIX_NORM(&J_b);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&J_b,(1.0/aux));

	// K_b = cross(I_b,J_b);
	GMATRIX_CROSS_COPY(K_b,I_b,J_b);

	// %Inicializa as variaveis do sistema inercial
	// aux = G + M;
	GMATRIX_ADD_COPY(I_i,G,M);
	// I_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
	aux = PGMATRIX_NORM(&I_i);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&I_i,(1.0/aux));

	// aux = cross(I_i,(G - M));
	GMATRIX_SUBTRACT_COPY(MatDummy,G,M);
	GMATRIX_CROSS_COPY(J_i,I_i,MatDummy);

	// J_i = aux/sqrt(aux(1)^2 + aux(2)^2 + aux(3)^2);
	aux = PGMATRIX_NORM(&J_i);
	if(aux==0.0) return 0;
	PGMATRIX_MULTIPLY_CONST(&J_i,(1.0/aux));

	// K_i = cross(I_i,J_i);
	GMATRIX_CROSS_COPY(K_i,I_i,J_i);

	// Calculo da matriz de rotação pelo método TRIAD melhorado
	// R_i_b = ([I_b, J_b, K_b]*([I_i, J_i, K_i]'));
	GMATRIX_SUBMATRIX_COPY(R_b,1,1,I_b);
	GMATRIX_SUBMATRIX_COPY(R_b,1,2,J_b);
	GMATRIX_SUBMATRIX_COPY(R_b,1,3,K_b);
	GMATRIX_SUBMATRIX_COPY(R_i,1,1,I_i);
	GMATRIX_SUBMATRIX_COPY(R_i,1,2,J_i);
	GMATRIX_SUBMATRIX_COPY(R_i,1,3,K_i);
	GMATRIX_MULTIPLY_COPY_EXTENDED(R_b_i, R_i, 0, R_b, 1);

	// %Criação do quatérnio de rotação. (Deve se observar que a matriz de
	// %rotação é do sistema de referência para o sistema do corpo a contrário do
	// %que o Padilha escreveu. Dados para a conferência estão em
	// %http://www.uel.br/proppg/semina/pdf/semina_28_1_22_19.pdf).
	// R_b_i = R_i_b';
	// q = dcm2quaternions(R_b_i);
	rotation_dcm2quaternions(&R_b_i, pq);

	// if FlagConsiderPrediction,
	//	q = quaternions_correctsign(q, q_predicted);
	//end
	if(pq_previous != NULL){
		rotation_quaternionscorrectsign(pq, pq_previous);
	}

	// Retorna 
    return 1; 
}                      

