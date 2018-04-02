/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: magnetometer.h 
Conteudo: Cabeçalho de funções em código C relacionadas ao magnetômetro.
Autor: G. A. Borges.
Atualizações: 
	- 28/07/2008: criação do exemplo, por Geovany A. Borges
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
