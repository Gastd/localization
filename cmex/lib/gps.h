#ifndef GPS_H
#define GPS_H

typedef struct{
	int FlagValidPositionMeasure;
	int FlagValidVelocityMeasure;
	PGMATRIX pPosition;
	PGMATRIX pVelocity;
	PGMATRIX pPPosition;
	PGMATRIX pPVelocity;
} GPSMEASURE, *PGPSMEASURE;

int gps_init(PGPSMEASURE pGPSMeasure);
int gps_close(PGPSMEASURE pGPSMeasure);

#endif //GPS_H
