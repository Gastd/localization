/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: sonar.c 
Conteudo: Fun��es em c�digo C relacionadas ao SONAR.
Autor: G. A. Borges.
Atualiza��es: 
	- 19/09/2008: cria��o por Geovany A. Borges, modulo sem init/close
*****************************************************************************/
// Cabecalhos de bibliotecas C run-time
#include <math.h>
#include <stdio.h>

// Cabecalhos de modulos do projeto
#include "gmatrix.h"

// Definicoes de uso interno

// Prototipos de funcoes internas ao modulo

// Variaveis globais do m�dulo

/*****************************************************************************
******************************************************************************
** FUNCOES COM CHAMADA EXTERNA
******************************************************************************
*****************************************************************************/

/*****************************************************************************
*** int sonar_init(void)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int sonar_init(PSONARMEASURE pSonarMeasure)
{
	pSonarMeasure->FlagValidMeasure = 0;
	pSonarMeasure->range = 0;
	pSonarMeasure->rangevariance = GMATRIXMACRO_SQR(0.01);

	pSonarMeasure->pR_s2b = PGMATRIX_ALLOC(3,3);
	pSonarMeasure->pt_s2b = PGMATRIX_ALLOC(3,1);

	PGMATRIX_IDENTITY(pSonarMeasure->pR_s2b);
	PGMATRIX_ZEROES(pSonarMeasure->pt_s2b);

	// Retorna 
    return 1; 
}                      

/*****************************************************************************
*** int sonar_close(PSONARMEASURE pSonarMeasure)
*** Entradas: 
***	Saidas:
*****************************************************************************/
int sonar_close(PSONARMEASURE pSonarMeasure)
{
	PGMATRIX_FREE(pSonarMeasure->pR_s2b);
	PGMATRIX_FREE(pSonarMeasure->pt_s2b);

	// Retorna 
    return 1; 
}                      

