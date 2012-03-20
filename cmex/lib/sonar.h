#ifndef SONAR_H
#define SONAR_H

typedef struct{
	int FlagValidMeasure;
	double range;
	double rangevariance;
	PGMATRIX pR_s2b; // Matriz de rota��o relacionando a medida rs = [range 0 0]' ao sistema de coordenadas do corpo: rb = R_s2b * rs + t_s2b; 
	PGMATRIX pt_s2b; // Vetor de tranla��o relacionando a medida rs = [range 0 0]' ao sistema de coordenadas do corpo.
} SONARMEASURE, *PSONARMEASURE;

int sonar_init(PSONARMEASURE pSonarMeasure);
int sonar_close(PSONARMEASURE pSonarMeasure);

#endif //SONAR_H
