/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: gps.c 
Conteudo: Funções em código C relacionadas ao GPS.
Autor: G. A. Borges.
Atualizações: 
	- 19/09/2008: criação por Geovany A. Borges, modulo sem init/close
*****************************************************************************/
// Cabecalhos de bibliotecas C run-time
#include <math.h>
#include <stdio.h>

// Cabecalhos de modulos do projeto
#include "gmatrix.h"

// Definicoes de uso interno

// Prototipos de funcoes internas ao modulo

// Variaveis globais do módulo

/*****************************************************************************
******************************************************************************
** FUNCOES COM CHAMADA EXTERNA
******************************************************************************
*****************************************************************************/

/*****************************************************************************
*** int gps_init(PGPSMEASURE pGPSMeasure);
*** Entradas: 
***	Saidas:
*****************************************************************************/
int gps_init(PGPSMEASURE pGPSMeasure)
{
	pGPSMeasure->FlagValidPositionMeasure = 0;
	pGPSMeasure->FlagValidVelocityMeasure = 0;
	pGPSMeasure->pPosition = PGMATRIX_ALLOC(3,1);
	pGPSMeasure->pVelocity = PGMATRIX_ALLOC(3,1);
	pGPSMeasure->pPPosition = PGMATRIX_ALLOC(3,3);
	pGPSMeasure->pPVelocity = PGMATRIX_ALLOC(3,3);

	PGMATRIX_ZEROES(pGPSMeasure->pPPosition);
	PGMATRIX_DATA(pGPSMeasure->pPPosition,1,1) = GMATRIXMACRO_SQR(2);
	PGMATRIX_DATA(pGPSMeasure->pPPosition,2,2) = GMATRIXMACRO_SQR(2);
	PGMATRIX_DATA(pGPSMeasure->pPPosition,3,3) = GMATRIXMACRO_SQR(3);

	PGMATRIX_ZEROES(pGPSMeasure->pPVelocity);
	PGMATRIX_DATA(pGPSMeasure->pPVelocity,1,1) = GMATRIXMACRO_SQR(0.2);
	PGMATRIX_DATA(pGPSMeasure->pPVelocity,2,2) = GMATRIXMACRO_SQR(0.2);
	PGMATRIX_DATA(pGPSMeasure->pPVelocity,3,3) = GMATRIXMACRO_SQR(0.3);

	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int gps_close(PGPSMEASURE pGPSMeasure);
*** Entradas: 
***	Saidas:
*****************************************************************************/
int gps_close(PGPSMEASURE pGPSMeasure)
{
	PGMATRIX_FREE(pGPSMeasure->p);
	PGMATRIX_FREE(pGPSMeasure->v);
	PGMATRIX_FREE(pGPSMeasure->Pp);
	PGMATRIX_FREE(pGPSMeasure->Pv);

	// Retorna 
    return 1; 
}                      
