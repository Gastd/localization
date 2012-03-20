#include <math.h>
#include <stdio.h>

#include "gmatrix.h"
#include "gmatrix_linalg.h"

void kalman_KF_update_innovationform(PGMATRIX pX, PGMATRIX pXpredicted, PGMATRIX pV, PGMATRIX pP, PGMATRIX pPpredicted, PGMATRIX pR, PDUMMY_MATRICES pDummy, int FlagUpdateCovariances)
{
	//  K = Ppredicted*inv(Ppredicted+R)
	//	X = Xpredicted + K*(v);
	//  P = Ppredicted - K*Ppredicted;

	PGMATRIX_ADD_COPY(pDummy->pMat2,pPpredicted,pR);
	PGMATRIX_INVERSE(pDummy->pMat2);
	PGMATRIX_MULTIPLY_COPY(pDummy->pMat1,pPpredicted,pDummy->pMat2); // pDummy->pMat1 = Kalman Gain.
	PGMATRIX_MULTIPLY_COPY(pDummy->pMat2,pDummy->pMat1,pV);
	PGMATRIX_ADD_COPY(pX,pXpredicted,pDummy->pMat2);

	if (!FlagUpdateCovariances){
		PGMATRIX_COPY(pP,pPpredicted);
	}
	else{
		PGMATRIX_MULTIPLY_COPY(pDummy->pMat2,pDummy->pMat1,pPpredicted);
		PGMATRIX_SUBTRACT_COPY(pP,pPpredicted,pDummy->pMat2);
	}

	return;
}

void kalman_EKF_update_innovationform(PGMATRIX pX, PGMATRIX pXpredicted, PGMATRIX pV, PGMATRIX pP, PGMATRIX pPpredicted, PGMATRIX pR, PGMATRIX pH, PDUMMY_MATRICES pDummy, int FlagUpdateCovariances)
{
	//  K = Ppredicted*H'*inv(H*Ppredicted*H'+R)
	//	X = Xpredicted + K*(V);
	//  P = Ppredicted - K*H*Ppredicted;

	PGMATRIX_MULTIPLY_COPY(pDummy->pMat1,pH,pPpredicted);
	PGMATRIX_MULTIPLY_COPY_EXTENDED(pDummy->pMat2, pDummy->pMat1, 0, pH, 1);
	PGMATRIX_ADD(pDummy->pMat2,pR);
	PGMATRIX_INVERSE(pDummy->pMat2);
	PGMATRIX_MULTIPLY_COPY_EXTENDED(pDummy->pMat3,pH,1,pDummy->pMat2,0);
	PGMATRIX_MULTIPLY_COPY(pDummy->pMat1,pPpredicted,pDummy->pMat3); // pDummy->pMat1 = Kalman Gain.

//	PGMATRIX_PRINT_MATLABFORM(pDummy->pMat1); // pDummy->pMat1 = Kalman Gain.

	PGMATRIX_MULTIPLY_COPY(pDummy->pMat2,pDummy->pMat1,pV);
	PGMATRIX_ADD_COPY(pX,pXpredicted,pDummy->pMat2);

	if (!FlagUpdateCovariances){
		PGMATRIX_COPY(pP,pPpredicted);
	}
	else{
		PGMATRIX_MULTIPLY_COPY(pDummy->pMat2,pDummy->pMat1,pH);
		PGMATRIX_MULTIPLY_COPY(pDummy->pMat3,pDummy->pMat2,pPpredicted);
		PGMATRIX_SUBTRACT_COPY(pP,pPpredicted,pDummy->pMat3);
	}
}
