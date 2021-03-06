/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: GPS.h 
Conteudo: Cabe�alho de fun��es em c�digo C relacionadas ao GPS.
Autor: G. A. Borges.
Atualiza��es: 
	- 28/07/2008: cria��o do exemplo, por Geovany A. Borges
*****************************************************************************/

#ifndef GPS_H
#define GPS_H

// Definicoes de uso externo
typedef struct{
	int FlagValidPositionMeasure;
	int FlagValidVelocityMeasure;
	PGMATRIX pPosition;
	PGMATRIX pVelocity;
	PGMATRIX pPPosition;
	PGMATRIX pPVelocity;
} GPSMEASURE, *PGPSMEASURE;

// Prototipos externos:
int gps_init(PGPSMEASURE pGPSMeasure);
int gps_close(PGPSMEASURE pGPSMeasure);


#endif //GPS_H
