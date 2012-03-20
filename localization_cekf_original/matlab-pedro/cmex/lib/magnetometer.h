/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: magnetometer.h 
Conteudo: Cabe�alho de fun��es em c�digo C relacionadas ao magnet�metro.
Autor: G. A. Borges.
Atualiza��es: 
	- 28/07/2008: cria��o do exemplo, por Geovany A. Borges
*****************************************************************************/

#ifndef MAGNETOMETER_H
#define MAGNETOMETER_H

// Definicoes de uso externo
typedef struct{
	double mx;
	double my;
	double mz;
	double mxvariance;
	double myvariance;
	double mzvariance;
	int    FlagValidMeasure;
} MAGNETOMETERMEASURE, *PMAGNETOMETERMEASURE;

// Prototipos externos:
int magnetometer_init(void);

#endif //MAGNETOMETER_H
