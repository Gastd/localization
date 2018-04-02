/**************************************************************
 Programa: mexlocalization.c
 Objetivo: Versao Cmex de localization
**************************************************************/
/* A ne pas oublier: x[i][j] (Format C) <=> x[i+j*M] (Format CMEX) */

#include "mex.h"
#include <stdio.h>
#include <math.h>

#ifdef SYSTEM_WINDOWS
  #include "WINPERFTime.h"
#endif 
#ifdef SYSTEM_LINUX
  #include <sys/time.h>
  #include <unistd.h> /* for libc5 */
  #include <sys/io.h> /* for glibc */
#endif 

#define GMATRIX_PRINTCOMMAND mexPrintf

#include "gmatrix.h"
#include "gmatrix_matlab.h"
#include "imu.h"
#include "magnetometer.h"
#include "gps.h"
#include "sonar.h"
#include "rotation.h"
#include "localization.h"

mxArray *SetFilterStruct(PLOCALIZATIONFILTERSTRUCT pFilterStruct)
{
    const char *field_names_ekf2[] = {"flagestimateaccelerometerbias","Nstates","X","P","Preset" ,"Q_ekf",
                                     "R_pseudomeasurementnorm","R_convertedmeasurementtriad_ekf"};
    const char *field_names_ukf2[] = {"flagestimateaccelerometerbias","Nstates","X","P","Preset" ,"Xsigma" ,"Wsigma","Q_ukf",
                                     "R_pseudomeasurementnorm","R_convertedmeasurementtriad_ukf"};
    const char *field_names_cekf[] = {"flagestimateaccelerometerbias","Nstates","X","P","Preset" ,"Q_cekf",
                                     "R_pseudomeasurementnorm","R_convertedmeasurementtriad_cekf",
                                     "S_previous_cekf","R_previous_cekf","A_imu_P_imu_cekf","innovation_previous_cekf"};
    mxArray *pmxarray;
    int dims[2] = {1, 1};

    switch(pFilterStruct->AlgorithmCode){
    case LOCALIZATION_ALGORITHMCODE_EKF2:
        pmxarray = mxCreateStructArray(2, dims, (sizeof(field_names_ekf2)/sizeof(*field_names_ekf2)), field_names_ekf2);
        break;
    case LOCALIZATION_ALGORITHMCODE_UKF2:
        pmxarray = mxCreateStructArray(2, dims, (sizeof(field_names_ukf2)/sizeof(*field_names_ukf2)), field_names_ukf2);
        break;
    case LOCALIZATION_ALGORITHMCODE_CEKF:
        pmxarray = mxCreateStructArray(2, dims, (sizeof(field_names_cekf)/sizeof(*field_names_cekf)), field_names_cekf);
        break;
    }

    mxSetField(pmxarray, 0, "flagestimateaccelerometerbias", mxCreateDoubleScalar((double)(pFilterStruct->FlagEstimateAccelerometerBias)));  
    mxSetField(pmxarray, 0, "Nstates", mxCreateDoubleScalar((double)(pFilterStruct->Nstates)));  
    mxSetField(pmxarray, 0, "X", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pX));  
    mxSetField(pmxarray, 0, "P", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pP));  
    mxSetField(pmxarray, 0, "Preset", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pPreset));  
    
    switch(pFilterStruct->AlgorithmCode){
    case LOCALIZATION_ALGORITHMCODE_EKF2:
        mxSetField(pmxarray, 0, "Q_ekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pQ));  
        mxSetField(pmxarray, 0, "R_pseudomeasurementnorm", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pR_pseudomeasurementnorm));  
        mxSetField(pmxarray, 0, "R_convertedmeasurementtriad_ekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pR_convertedmeasurementtriad));  
        break;
    case LOCALIZATION_ALGORITHMCODE_UKF2:
        mxSetField(pmxarray, 0, "Q_ukf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pQ));  
        mxSetField(pmxarray, 0, "Xsigma", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pXsigma_ukf));  
        mxSetField(pmxarray, 0, "Wsigma", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pWsigma_ukf));  
        mxSetField(pmxarray, 0, "R_pseudomeasurementnorm", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pR_pseudomeasurementnorm));  
        mxSetField(pmxarray, 0, "R_convertedmeasurementtriad_ukf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pR_convertedmeasurementtriad));  
        break;
    case LOCALIZATION_ALGORITHMCODE_CEKF:
        mxSetField(pmxarray, 0, "Q_cekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pQ));  
        mxSetField(pmxarray, 0, "R_pseudomeasurementnorm", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pR_pseudomeasurementnorm));  
        mxSetField(pmxarray, 0, "R_convertedmeasurementtriad_cekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pR_convertedmeasurementtriad));  
        mxSetField(pmxarray, 0, "S_previous_cekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pS_previous_cekf));  
        mxSetField(pmxarray, 0, "R_previous_cekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pR_previous_cekf));  
        mxSetField(pmxarray, 0, "A_imu_P_imu_cekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pA_imu_P_imu_cekf));  
        mxSetField(pmxarray, 0, "innovation_previous_cekf", PGMATRIX_COPY_TO_MXARRAY(pFilterStruct->pinnovation_previous_cekf));  
        break;
    }
    
    return pmxarray;
}

void GetFilterStruct(int AlgorithmCode, PLOCALIZATIONFILTERSTRUCT pFilterStruct, const mxArray *pmxarray)
{
    int FlagPrint = 0;
//    mxArray *pmxdataarray;
    
    if(!mxIsStruct(pmxarray)){
        mexErrMsgTxt("pmxarray should be a structure.");
    }
   	pFilterStruct->AlgorithmCode = AlgorithmCode;
    pFilterStruct->FlagEstimateAccelerometerBias = (int) mxGetScalar(mxGetField(pmxarray, 0, "flagestimateaccelerometerbias")); 
    pFilterStruct->Nstates = (int) mxGetScalar(mxGetField(pmxarray, 0, "Nstates")); 
    if (FlagPrint){
        printf("\npFilterStruct->AlgorithmCode = %i",pFilterStruct->AlgorithmCode);
        printf("\npFilterStruct->FlagEstimateAccelerometerBias = %i",pFilterStruct->FlagEstimateAccelerometerBias);
        printf("\npFilterStruct->Nstates = %i",pFilterStruct->Nstates);
    }

    // aloca a estrutura
    pFilterStruct->pX = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "X")); 
    pFilterStruct->pP = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "P")); 
    pFilterStruct->pPreset = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "Preset")); 
    if (FlagPrint){
        PGMATRIX_PRINT_MATLABFORM(pFilterStruct->pX);
        PGMATRIX_PRINT_MATLABFORM(pFilterStruct->pP);
        PGMATRIX_PRINT_MATLABFORM(pFilterStruct->pPreset);
    }
    if(AlgorithmCode==LOCALIZATION_ALGORITHMCODE_EKF2){
        pFilterStruct->pQ = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "Q_ekf")); 
        pFilterStruct->pR_convertedmeasurementtriad = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "R_convertedmeasurementtriad_ekf")); 
    }
    if(AlgorithmCode==LOCALIZATION_ALGORITHMCODE_CEKF){
        pFilterStruct->pQ = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "Q_cekf")); 
        pFilterStruct->pR_convertedmeasurementtriad = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "R_convertedmeasurementtriad_cekf")); 
        pFilterStruct->pS_previous_cekf = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "S_previous_cekf")); 
        pFilterStruct->pR_previous_cekf = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "R_previous_cekf")); 
        pFilterStruct->pA_imu_P_imu_cekf = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "A_imu_P_imu_cekf")); 
        pFilterStruct->pinnovation_previous_cekf = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "innovation_previous_cekf")); 
    }
    if(AlgorithmCode==LOCALIZATION_ALGORITHMCODE_UKF2){
        pFilterStruct->pQ = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "Q_ukf")); 
        pFilterStruct->pXsigma_ukf = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "Xsigma")); 
        pFilterStruct->pWsigma_ukf = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "Wsigma")); 
        pFilterStruct->pR_convertedmeasurementtriad = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "R_convertedmeasurementtriad_ukf")); 
    }
    pFilterStruct->pR_pseudomeasurementnorm = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "R_pseudomeasurementnorm")); 
    if (FlagPrint){
        PGMATRIX_PRINT_MATLABFORM(pFilterStruct->pQ);
        PGMATRIX_PRINT_MATLABFORM(pFilterStruct->pR_convertedmeasurementtriad);
        PGMATRIX_PRINT_MATLABFORM(pFilterStruct->pR_pseudomeasurementnorm);
    }
}

void GetGPSMeasure(PGPSMEASURE pGPSMeasure, const mxArray *pmxarray)
{
    if(!mxIsStruct(pmxarray)){
        mexErrMsgTxt("pmxarray should be a structure.");
    }

    pGPSMeasure->FlagValidPositionMeasure = (int) mxGetScalar(mxGetField(pmxarray, 0, "flagvalidpmeasure")); 
    pGPSMeasure->FlagValidVelocityMeasure = (int) mxGetScalar(mxGetField(pmxarray, 0, "flagvalidvmeasure")); 
    pGPSMeasure->pPosition = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "p")); 
    pGPSMeasure->pVelocity = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "v")); 
    pGPSMeasure->pPPosition = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "P_p")); 
    pGPSMeasure->pPVelocity = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "P_v")); 
}

void GetSonarMeasure(PSONARMEASURE pSonarMeasure, const mxArray *pmxarray)
{
    PGMATRIX pMatrix;
    
    if(!mxIsStruct(pmxarray)){
        mexErrMsgTxt("pmxarray should be a structure.");
    }

    pSonarMeasure->FlagValidMeasure = (int) mxGetScalar(mxGetField(pmxarray, 0, "flagvalidmeasure")); 
    pSonarMeasure->range = mxGetScalar(mxGetField(pmxarray, 0, "range")); 
    pSonarMeasure->rangevariance = mxGetScalar(mxGetField(pmxarray, 0, "rangevariance"));
    pSonarMeasure->pR_s2b = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "R_s2b")); 
    pSonarMeasure->pt_s2b = PGMATRIX_ALLOC_FROM_MXARRAY(mxGetField(pmxarray, 0, "t_s2b")); 
}

void GetIMUMeasure(PIMUMEASURE pIMUMeasure, const mxArray *pmxarray)
{
    if(!mxIsStruct(pmxarray)){
        mexErrMsgTxt("pmxarray should be a structure.");
    }
    pIMUMeasure->ax = mxGetScalar(mxGetField(pmxarray, 0, "ax")); 
    pIMUMeasure->ay = mxGetScalar(mxGetField(pmxarray, 0, "ay")); 
    pIMUMeasure->az = mxGetScalar(mxGetField(pmxarray, 0, "az")); 
    pIMUMeasure->wx = mxGetScalar(mxGetField(pmxarray, 0, "wx")); 
    pIMUMeasure->wy = mxGetScalar(mxGetField(pmxarray, 0, "wy")); 
    pIMUMeasure->wz = mxGetScalar(mxGetField(pmxarray, 0, "wz")); 
    pIMUMeasure->axvariance = mxGetScalar(mxGetField(pmxarray, 0, "axvariance")); 
    pIMUMeasure->ayvariance = mxGetScalar(mxGetField(pmxarray, 0, "ayvariance")); 
    pIMUMeasure->azvariance = mxGetScalar(mxGetField(pmxarray, 0, "azvariance")); 
    pIMUMeasure->wxvariance = mxGetScalar(mxGetField(pmxarray, 0, "wxvariance")); 
    pIMUMeasure->wyvariance = mxGetScalar(mxGetField(pmxarray, 0, "wyvariance")); 
    pIMUMeasure->wzvariance = mxGetScalar(mxGetField(pmxarray, 0, "wzvariance")); 
}

void GetMagnetometerMeasure(PMAGNETOMETERMEASURE pMagnetometerMeasure, const mxArray *pmxarray)
{
    if(!mxIsStruct(pmxarray)){
        mexErrMsgTxt("pmxarray should be a structure.");
    }
    pMagnetometerMeasure->mx = mxGetScalar(mxGetField(pmxarray, 0, "mx")); 
    pMagnetometerMeasure->my = mxGetScalar(mxGetField(pmxarray, 0, "my")); 
    pMagnetometerMeasure->mz = mxGetScalar(mxGetField(pmxarray, 0, "mz")); 
    pMagnetometerMeasure->mxvariance = mxGetScalar(mxGetField(pmxarray, 0, "mxvariance")); 
    pMagnetometerMeasure->myvariance = mxGetScalar(mxGetField(pmxarray, 0, "myvariance")); 
    pMagnetometerMeasure->mzvariance = mxGetScalar(mxGetField(pmxarray, 0, "mzvariance")); 
    pMagnetometerMeasure->FlagValidMeasure = (int) mxGetScalar(mxGetField(pmxarray, 0, "flagvalidmeasure")); 
}

PGMATRIX GetM(const mxArray *pmxarray)
{
   	if ( (mxGetM(pmxarray)!=3) ){
		mexErrMsgTxt("M must have 3 rows.");
	}
   	if ( (mxGetN(pmxarray)!=1) ){
		mexErrMsgTxt("M must have 1 column.");
	}
    return(PGMATRIX_ALLOC_FROM_MXARRAY(pmxarray));
}

PGMATRIX GetG(const mxArray *pmxarray)
{
   	if ( (mxGetM(pmxarray)!=3) ){
		mexErrMsgTxt("G must have 3 rows.");
	}
   	if ( (mxGetN(pmxarray)!=1) ){
		mexErrMsgTxt("G must have 1 column.");
	}
    return(PGMATRIX_ALLOC_FROM_MXARRAY(pmxarray));
}

PQUATERNIONS Getq(const mxArray *pmxarray)
{
   	if ( (mxGetM(pmxarray)!=4) ){
		mexErrMsgTxt("q must have 4 rows.");
	}
   	if ( (mxGetN(pmxarray)!=1) ){
		mexErrMsgTxt("q must have 1 column.");
	}
    return(PGMATRIX_ALLOC_FROM_MXARRAY(pmxarray));
}

double GetT(const mxArray *pmxarray)
{
   	if ( (mxGetM(pmxarray)!=1) ){
		mexErrMsgTxt("T must have 1 rows.");
	}
   	if ( (mxGetN(pmxarray)!=1) ){
		mexErrMsgTxt("T must have 1 column.");
	}
    return(mxGetScalar(pmxarray));
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int		strsize;
	char *pLocalizationFunctionName;
    LOCALIZATIONFILTERSTRUCT FilterStruct;
    IMUMEASURE IMUMeasure;
    GPSMEASURE GPSMeasure;
    MAGNETOMETERMEASURE MagnetometerMeasure;
    SONARMEASURE SonarMeasure;
    PGMATRIX pM, pG, pq_previous;
    double T,Telapsed;
    QUATERNIONS_DECLARE(q);
#ifdef SYSTEM_WINDOWS
    WINPERF_COUNTER WinperfCounter;
#endif 
#ifdef SYSTEM_LINUX
    struct timeval     timereset;
    struct timeval     time;
#endif

	/* Processar as entradas */
   	strsize = (mxGetM(prhs[0])*mxGetN(prhs[0])) + 1;
	pLocalizationFunctionName = (char *) mxCalloc(strsize, sizeof(char));
	if(mxGetString(prhs[0],pLocalizationFunctionName,strsize)!=0){
		mexErrMsgTxt("Not enough space to allocate buffer for pLocalizationFunctionName.");
	}

    /* Algoritmo TRIAD Improved */
    if (strcmp(pLocalizationFunctionName,"TRIAD")==0){
        // imumeasure, magnetometermeasure, M, G, q_previous
        GetIMUMeasure(&IMUMeasure, prhs[1]);
        GetMagnetometerMeasure(&MagnetometerMeasure, prhs[2]);
        pM = GetM(prhs[3]);
        pG = GetG(prhs[4]);
        pq_previous = Getq(prhs[5]);
        
        // funï¿½ï¿½o principal:
#ifdef SYSTEM_WINDOWS
        WINPERF_ResetCounter(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&timereset, NULL);
#endif
        localization_triad(&q, &IMUMeasure, &MagnetometerMeasure, pM, pG, pq_previous);
#ifdef SYSTEM_WINDOWS
        Telapsed = WINPERF_GetElapsedTime(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&time, NULL);
	Telapsed = (time.tv_sec - timereset.tv_sec) + 1e-6 * (time.tv_usec - timereset.tv_usec);
#endif
        
        // saï¿½da:
        plhs[0] = PGMATRIX_COPY_TO_MXARRAY(&q);
        if(nlhs==2){
            plhs[1] = mxCreateDoubleScalar(Telapsed);
        }
        
        // liberar memoria:
        PGMATRIX_FREE(pG);
        PGMATRIX_FREE(pM);
        PGMATRIX_FREE(pq_previous);
    }

   /* % Algoritmo FKE aplicado ï¿½ arquitetura correlata: etapa de prediï¿½ï¿½o */
    if (strcmp(pLocalizationFunctionName,"FILTER_EKF2_PREDICTION")==0){
        // ekf_structure, imumeasure, G, T
        GetFilterStruct(LOCALIZATION_ALGORITHMCODE_EKF2, &FilterStruct, prhs[1]);
        GetIMUMeasure(&IMUMeasure, prhs[2]);
        pG = GetG(prhs[3]);
        T = GetT(prhs[4]);

        // funï¿½ï¿½o principal:
#ifdef SYSTEM_WINDOWS
        WINPERF_ResetCounter(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&timereset, NULL);
#endif
        localization_filter_prediction(&FilterStruct, &IMUMeasure, pG, T);
#ifdef SYSTEM_WINDOWS
        Telapsed = WINPERF_GetElapsedTime(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&time, NULL);
	Telapsed = (time.tv_sec - timereset.tv_sec) + 1e-6 * (time.tv_usec - timereset.tv_usec);
#endif
        
        // saï¿½da:
        plhs[0] = SetFilterStruct(&FilterStruct);
        if(nlhs==2){
            plhs[1] = mxCreateDoubleScalar(Telapsed);
        }
        
        // liberar memoria:
        PGMATRIX_FREE(pG);
        localization_close(&FilterStruct);
    }

   /* % Algoritmo FKEC aplicado ï¿½ arquitetura correlata: etapa de prediï¿½ï¿½o */
    if (strcmp(pLocalizationFunctionName,"FILTER_CEKF_PREDICTION")==0){
        // ekf_structure, imumeasure, G, T
        GetFilterStruct(LOCALIZATION_ALGORITHMCODE_CEKF, &FilterStruct, prhs[1]);
        GetIMUMeasure(&IMUMeasure, prhs[2]);
        pG = GetG(prhs[3]);
        T = GetT(prhs[4]);

        // funï¿½ï¿½o principal:
#ifdef SYSTEM_WINDOWS
        WINPERF_ResetCounter(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&timereset, NULL);
#endif
        localization_filter_prediction(&FilterStruct, &IMUMeasure, pG, T);
#ifdef SYSTEM_WINDOWS
        Telapsed = WINPERF_GetElapsedTime(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&time, NULL);
	Telapsed = (time.tv_sec - timereset.tv_sec) + 1e-6 * (time.tv_usec - timereset.tv_usec);
#endif
        
        // saï¿½da:
        plhs[0] = SetFilterStruct(&FilterStruct);
        if(nlhs==2){
            plhs[1] = mxCreateDoubleScalar(Telapsed);
        }
        
        // liberar memoria:
        PGMATRIX_FREE(pG);
        localization_close(&FilterStruct);
    }
    
    
   /* % Algoritmo FKE aplicado ï¿½ arquitetura correlata: etapa de correï¿½ï¿½o */
    if (strcmp(pLocalizationFunctionName,"FILTER_EKF2_CORRECTION")==0){
        // ekf_structure, imumeasure, G, T
        GetFilterStruct(LOCALIZATION_ALGORITHMCODE_EKF2, &FilterStruct, prhs[1]);
        GetGPSMeasure(&GPSMeasure, prhs[2]);
        GetIMUMeasure(&IMUMeasure, prhs[3]);
        GetMagnetometerMeasure(&MagnetometerMeasure, prhs[4]);
        GetSonarMeasure(&SonarMeasure, prhs[5]);
        pM = GetM(prhs[6]);
        pG = GetG(prhs[7]);
        T = GetT(prhs[8]);

        // funï¿½ï¿½o principal:
#ifdef SYSTEM_WINDOWS
        WINPERF_ResetCounter(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&timereset, NULL);
#endif
        localization_filter_correction(&FilterStruct, &GPSMeasure, &IMUMeasure, &MagnetometerMeasure, &SonarMeasure, pM, pG, T);
#ifdef SYSTEM_WINDOWS
        Telapsed = WINPERF_GetElapsedTime(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&time, NULL);
	Telapsed = (time.tv_sec - timereset.tv_sec) + 1e-6 * (time.tv_usec - timereset.tv_usec);
#endif
        
        // saï¿½da:
        plhs[0] = SetFilterStruct(&FilterStruct);
        if(nlhs==2){
            plhs[1] = mxCreateDoubleScalar(Telapsed);
        }
        
        // liberar memoria:
        PGMATRIX_FREE(pG);
        localization_close(&FilterStruct);
    }

   /* % Algoritmo FKEC aplicado ï¿½ arquitetura correlata: etapa de correï¿½ï¿½o */
    if (strcmp(pLocalizationFunctionName,"FILTER_CEKF_CORRECTION")==0){
        // ekf_structure, imumeasure, G, T
        GetFilterStruct(LOCALIZATION_ALGORITHMCODE_CEKF, &FilterStruct, prhs[1]);
        GetGPSMeasure(&GPSMeasure, prhs[2]);
        GetIMUMeasure(&IMUMeasure, prhs[3]);
        GetMagnetometerMeasure(&MagnetometerMeasure, prhs[4]);
        GetSonarMeasure(&SonarMeasure, prhs[5]);
        pM = GetM(prhs[6]);
        pG = GetG(prhs[7]);
        T = GetT(prhs[8]);

        // funï¿½ï¿½o principal:
#ifdef SYSTEM_WINDOWS
        WINPERF_ResetCounter(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&timereset, NULL);
#endif
        localization_filter_correction(&FilterStruct, &GPSMeasure, &IMUMeasure, &MagnetometerMeasure, &SonarMeasure, pM, pG, T);
#ifdef SYSTEM_WINDOWS
        Telapsed = WINPERF_GetElapsedTime(&WinperfCounter);
#endif 
#ifdef SYSTEM_LINUX
        gettimeofday(&time, NULL);
	Telapsed = (time.tv_sec - timereset.tv_sec) + 1e-6 * (time.tv_usec - timereset.tv_usec);
#endif
        
        // saï¿½da:
        plhs[0] = SetFilterStruct(&FilterStruct);
        if(nlhs==2){
            plhs[1] = mxCreateDoubleScalar(Telapsed);
        }
        
        // liberar memoria:
        PGMATRIX_FREE(pG);
        localization_close(&FilterStruct);
    }
}
	
