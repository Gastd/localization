#ifndef KALMAN_H
#define KALMAN_H

void kalman_KF_update_innovationform(PGMATRIX pX, PGMATRIX pXpredicted, PGMATRIX pV, PGMATRIX pP, PGMATRIX pPpredicted, PGMATRIX pR, PDUMMY_MATRICES pDummy, int FlagUpdateCovariances);
void kalman_EKF_update_innovationform(PGMATRIX pX, PGMATRIX pXpredicted, PGMATRIX pV, PGMATRIX pP, PGMATRIX pPpredicted, PGMATRIX pR, PGMATRIX pH, PDUMMY_MATRICES pDummy, int FlagUpdateCovariances);

#endif //KALMAN_H
