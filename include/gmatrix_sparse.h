/**********************************************************************************************
***********************************************************************************************
***  File:			Gmatrix_sparse.h
***	 Author:		Geovany Araujo Borges
***	 Contents:		Gmatrix_sparse header file.
***********************************************************************************************
**********************************************************************************************/

#ifndef GMATRIX_SPARSE_H
#define GMATRIX_SPARSE_H

/**********************************************************************************************
***** GMatrix: Type Definitions for specific sparse representations
**********************************************************************************************/

/*! 
 *  GMATRIX_SPARSE is the data structure of a matrix or a vector. 
 *  PGMATRIX_SPARSE is a pointer for GMATRIX_SPARSE.
 */

// Main GMATRIX_SPARSE type structure and pointer
typedef struct{	
		int N;		 /* Matrix dimension N x N */
		int MaxSize; /* Number of entries of *Data. */
		double *Data;  /* Pointer to entries of GMATRIX_SPARSE. The entries are row organized */	
		int *Index; /* Pointer to entries of GMATRIX_SPARSE. The entries are row organized */	
} GMATRIX_SPARSE, *PGMATRIX_SPARSE;

/**********************************************************************************************
***** GMatrix: Function prototypes
**********************************************************************************************/

#define GMATRIX_SPARSITY_COEFFICIENT(Mat,tolerance) PGMATRIX_SPARSITY_COEFFICIENT(&Mat,tolerance)
double PGMATRIX_SPARSITY_COEFFICIENT(PGMATRIX pMat, double tolerance);
PGMATRIX_PRINT_RAWFORM_NAMED(char* NameString, PGMATRIX_SPARSE pMat);
PGMATRIX_SPARSE PGMATRIX_SPARSE_ALLOC(int N);
PGMATRIX_SPARSE PGMATRIX_SPARSE_ALLOC_FROM_GMATRIX(PGMATRIX pMat, double threshold);
void PGMATRIX_SPARSE_FREE(PGMATRIX_SPARSE pMatrixSparse);

#endif //GMATRIX_SPARSE_H
