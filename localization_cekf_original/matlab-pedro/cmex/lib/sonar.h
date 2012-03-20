/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: sonar.h 
Conteudo: Cabeçalho de funções em código C relacionadas ao sonar.
Autor: G. A. Borges.
Atualizações: 
	- 28/07/2008: criação do exemplo, por Geovany A. Borges
*****************************************************************************/

#ifndef SONAR_H
#define SONAR_H

// Definicoes de uso externo
typedef struct{
	int FlagValidMeasure;
	double range;
	double rangevariance;
	PGMATRIX pR_s2b; // Matriz de rotação relacionando a medida rs = [range 0 0]' ao sistema de coordenadas do corpo: rb = R_s2b * rs + t_s2b; 
	PGMATRIX pt_s2b; // Vetor de tranlação relacionando a medida rs = [range 0 0]' ao sistema de coordenadas do corpo.
} SONARMEASURE, *PSONARMEASURE;

// Prototipos externos:
int sonar_init(void);

#endif //SONAR_H
