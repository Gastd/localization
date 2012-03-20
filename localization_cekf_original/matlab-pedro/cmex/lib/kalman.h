/*****************************************************************************
Arquivo: kalman.h 
Conteudo: Cabe�alho de fun��es em c�digo C relacionadas ao filtro de kalman.
Autor: G. A. Borges.
Atualiza��es: 
	- 28/07/2008: cria��o do exemplo, por Geovany A. Borges
*****************************************************************************/

#ifndef KALMAN_H
#define KALMAN_H

// Definicoes de uso externo

// Prototipos externos:
void kalman_KF_update_innovationform(PGMATRIX pX, PGMATRIX pXpredicted, PGMATRIX pV, PGMATRIX pP, PGMATRIX pPpredicted, PGMATRIX pR, PDUMMY_MATRICES pDummy, int FlagUpdateCovariances);
void kalman_EKF_update_innovationform(PGMATRIX pX, PGMATRIX pXpredicted, PGMATRIX pV, PGMATRIX pP, PGMATRIX pPpredicted, PGMATRIX pR, PGMATRIX pH, PDUMMY_MATRICES pDummy, int FlagUpdateCovariances);

#endif //KALMAN_H
