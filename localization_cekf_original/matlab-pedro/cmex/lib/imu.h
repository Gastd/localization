/*****************************************************************************
Projeto CARCARAH (UnB-Expansion)
Arquivo: imu.h 
Conteudo: Cabeçalho de funções em código C relacionadas ao IMU.
Autor: G. A. Borges.
Atualizações: 
	- 28/07/2008: criação do exemplo, por Geovany A. Borges
*****************************************************************************/

#ifndef IMU_H
#define IMU_H

// Definicoes de uso externo
typedef struct{
    double ax;
    double ay;
    double az;
    double wx;
    double wy;
    double wz;
    double axvariance;
    double ayvariance;
    double azvariance;
    double wxvariance;
    double wyvariance;
    double wzvariance;
} IMUMEASURE, *PIMUMEASURE;


// Prototipos externos:
int imu_init(void);

#endif //IMU_H
