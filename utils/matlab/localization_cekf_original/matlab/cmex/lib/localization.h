/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: localization.h 
Conteudo: Cabeçalho de funções em código C relacionadas à localização.
Autor: G. A. Borges.
Atualizações: 
	- 28/07/2008: criação do exemplo, por Geovany A. Borges
*****************************************************************************/

#ifndef LOCALIZATION_H
#define LOCALIZATION_H

// Definicoes de uso externo
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
	PGMATRIX pS_previous_cekf; // para uso com o CEKF
	PGMATRIX pR_previous_cekf; // para uso com o CEKF
	PGMATRIX pA_imu_P_imu_cekf; // para uso com o CEKF
	PGMATRIX pinnovation_previous_cekf; // para uso com o CEKF

}LOCALIZATIONFILTERSTRUCT, *PLOCALIZATIONFILTERSTRUCT;

#define X_q0		1
#define X_q1		2
#define X_q2		3
#define X_q3		4
#define X_px		5
#define X_py		6
#define X_pz		7
#define X_vx		8
#define X_vy		9
#define X_vz		10
#define X_bax		11
#define X_bay		12
#define X_baz		13

#define U_ax		1
#define U_ay		2
#define U_az		3
#define U_wx		4
#define U_wy		5
#define U_wz		6

#define LOCALIZATION_MAXSTATESIZE			13

#define LOCALIZATION_ALGORITHMCODE_EKF2		1
#define LOCALIZATION_ALGORITHMCODE_CEKF		2
#define LOCALIZATION_ALGORITHMCODE_UKF2		3

// Prototipos externos:
int localization_init(int AlgorithmCode, int FlagEstimateAccelerometerBias, PLOCALIZATIONFILTERSTRUCT pFilterStruct);
int localization_close(PLOCALIZATIONFILTERSTRUCT pFilterStruct);
int localization_triad(PQUATERNIONS pq, PIMUMEASURE pIMUMeasure, PMAGNETOMETERMEASURE pMagnetometerMeasure, PGMATRIX pM, PGMATRIX pG, PQUATERNIONS pq_previous);
int localization_filter_prediction(PLOCALIZATIONFILTERSTRUCT pFilterStruct, PIMUMEASURE pIMUMeasure, PGMATRIX pG, double T);

#endif //LOCALIZATION_H
