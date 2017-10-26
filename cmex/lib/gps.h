#ifndef GPS_H
#define GPS_H

typedef struct{
	int FlagValidPositionMeasure;
	int FlagValidVelocityMeasure;
	PGMATRIX pPosition;
	PGMATRIX pVelocity;
	PGMATRIX pPPosition;
	PGMATRIX pPVelocity;
} GpsMeasure;

int gps_init(GpsMeasure *gps_measure_ptr);
int gps_close(GpsMeasure *gps_measure_ptr);

#endif //GPS_H
