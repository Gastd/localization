/**********************************************************************************************
***********************************************************************************************
***  File:			GMatrix_linalg.h
***	 Author:		Geovany Araujo Borges
***	 Contents:		GMatrix_linalg header file.
***********************************************************************************************
**********************************************************************************************/

#ifndef GMATRIX_LINALG_H
#define GMATRIX_LINALG_H

/**********************************************************************************************
***** GMatrix: Function prototypes
**********************************************************************************************/

#define GMATRIX_LUDCMP(MatLU, Mat) PGMATRIX_LUDCMP(&MatLU, &Mat)
double PGMATRIX_LUDCMP(PGMATRIX pMatLU, PGMATRIX pMat);

#define GMATRIX_GAUSSJORDAN(AInv, X, A, B) PGMATRIX_GAUSSJORDAN(&AInv, &X, &A, &B)
void PGMATRIX_GAUSSJORDAN(PGMATRIX pAInverse, PGMATRIX pX, PGMATRIX pA, PGMATRIX pB);

#define GMATRIX_CHOLESKY(L, A) PGMATRIX_CHOLESKY(&L, &A)
void PGMATRIX_CHOLESKY(PGMATRIX pL, PGMATRIX pA);

#define GMATRIX_DETERMINANT2(Mat) PGMATRIX_DETERMINANT2(&Mat)
double PGMATRIX_DETERMINANT2(PGMATRIX pMat);

#define GMATRIX_DETERMINANT3(Mat) PGMATRIX_DETERMINANT3(&Mat)
double PGMATRIX_DETERMINANT3(PGMATRIX pMat);

#define GMATRIX_DETERMINANT(Mat, MatDummy) PGMATRIX_DETERMINANT(&Mat, &MatDummy)
double PGMATRIX_DETERMINANT(PGMATRIX pMat, PGMATRIX pMatDummy);

#define GMATRIX_NORM(Mat) PGMATRIX_NORM(&Mat)
double PGMATRIX_NORM(PGMATRIX pMat);

#define GMATRIX_INVERSE2_COPY(MatInverse, Mat) PGMATRIX_INVERSE2_COPY(&MatInverse, &Mat)
void PGMATRIX_INVERSE2_COPY(PGMATRIX pMatInverse, PGMATRIX pMat); 

#define GMATRIX_INVERSE3_COPY(MatInverse, Mat) PGMATRIX_INVERSE3_COPY(&MatInverse, &Mat)
void PGMATRIX_INVERSE3_COPY(PGMATRIX pMatInverse, PGMATRIX pMat); 

#define	GMATRIX_INVERSE_COPY(MatInverse, Mat) PGMATRIX_INVERSE_COPY(&MatInverse, &Mat)
void PGMATRIX_INVERSE_COPY(PGMATRIX pMatInverse, PGMATRIX pMat); 

#define	GMATRIX_INVERSE(Mat) PGMATRIX_INVERSE(&Mat)
void PGMATRIX_INVERSE(PGMATRIX pMat); 

#define	GMATRIX_SVD(U,S,V,Mat,FlagSorted) PGMATRIX_SVD(&U,&S,&V,&Mat,FlagSorted)
void PGMATRIX_SVD(PGMATRIX pU,PGMATRIX pS,PGMATRIX pV,PGMATRIX pMat, unsigned char FlagSorted);

#define	GMATRIX_PSEUDOINVERSE(Apinv,A,MatDummy) PGMATRIX_PSEUDOINVERSE(&Apinv,&A,&MatDummy)
void PGMATRIX_PSEUDOINVERSE(PGMATRIX pApinv, PGMATRIX pA, PGMATRIX pMatDummy);

#define GMATRIX_LEFT_PSEUDOINVERSE_COPY(Apinv,A,MatDummy) PGMATRIX_LEFT_PSEUDOINVERSE_COPY(&Apinv,&A,&MatDummy)
void PGMATRIX_LEFT_PSEUDOINVERSE_COPY(PGMATRIX pApinv, PGMATRIX pA, PGMATRIX pMatDummy);

#define GMATRIX_DERIVATIVE_LEFT_PSEUDOINVERSE(dApinvdx,dAdx,A,Apinv,MatDummy1,MatDummy2) PGMATRIX_DERIVATIVE_LEFT_PSEUDOINVERSE(&dApinvdx,&dAdx,&A,&Apinv,&MatDummy1,&MatDummy2)
void PGMATRIX_DERIVATIVE_LEFT_PSEUDOINVERSE(PGMATRIX pdApinvdx, PGMATRIX pdAdx, PGMATRIX pA, PGMATRIX pApinv, PGMATRIX pMatDummy1, PGMATRIX pMatDummy2);

#define GMATRIX_RANK_FROMSVD(U,S,V) PGMATRIX_RANK_FROMSVD(&U,&S,&V)
int PGMATRIX_RANK_FROMSVD(PGMATRIX pU, PGMATRIX pS, PGMATRIX pV);

#define GMATRIX_RANK(Mat) PGMATRIX_RANK(&Mat)
int PGMATRIX_RANK(PGMATRIX pMat);

#endif //GMATRIX_LINALG_H
