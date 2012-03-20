#include <math.h>
#include <stdio.h>
#include "gmatrix.h"
#include "gmatrix_linalg.h"
#include "rotation.h"

int rotation_quaternionscorrectsign(PQUATERNIONS pqcurrent, PQUATERNIONS pqprevious)
{
	QUATERNIONS_DECLARE(plus);
	QUATERNIONS_DECLARE(minus);
	
	// now we check if we want q or -q
	// plus  = q_predicted - q;
	PGMATRIX_SUBTRACT_COPY(&plus, pqprevious, pqcurrent);
	// minus = q_predicted + q;
	PGMATRIX_ADD_COPY(&minus, pqprevious, pqcurrent);
	// if (norm(plus) > norm(minus))
	//		q = - q;
	// end
	if( GMATRIX_NORM(plus) > GMATRIX_NORM(minus) ){
		PGMATRIX_MULTIPLY_CONST(pqcurrent, -1.0);
	}
	
	// Retorna 
    return 1; 
}                      

int rotation_quaternions2euler(PQUATERNIONS pq, double *proll, double *ppitch, double *pyaw)
{
	double norm;
	DCM_DECLARE(R);

	// make quaternion unit length: done in quaternions2dcm
	//	R = quaternions2dcm(q);
	if(!rotation_quaternions2dcm(pq, &R)) return 0;

	//	[roll,pitch,yaw] = dcm2euler(R);
	if(!rotation_dcm2euler(&R, proll, ppitch, pyaw)) return 0;

	
	// Retorna 
    return 1; 
}                      

int rotation_euler2quaternions(double roll, double pitch, double yaw, PQUATERNIONS pq)
{
	DCM_DECLARE(R);

	// R = euler2dcm(roll,pitch,yaw);
	if(!rotation_euler2dcm(roll, pitch, yaw, &R)) return 0;

	// q = dcm2quaternions(R);
	if(!rotation_dcm2quaternions(&R, pq)) return 0;
	
	// Retorna 
    return 1; 
}                      

int rotation_quaternions2dcm(PQUATERNIONS pq, PDCM pR)
{
	double a;

	// make quaternion unit length:
	rotation_quaternionsnormalize(pq);

	// in acordance with Titterton & Weston, pg. 48, and Coutsias & Romero, last pg
	// R = 2*[(q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2)/2    (q(2)*q(3)-q(1)*q(4))                   (q(2)*q(4)+q(1)*q(3));
	PDCM_R11(pR) = 1.0*(PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq) + PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq) - PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq) - PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R12(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q2(pq) - PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R13(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q3(pq) + PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q2(pq));

	//	   (q(2)*q(3)+q(1)*q(4))                    (q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2)/2   (q(3)*q(4)-q(1)*q(2));
	PDCM_R21(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q2(pq) + PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R22(pR) = 1.0*(PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq) - PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq) + PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq) - PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq));
	PDCM_R23(pR) = 2.0*(PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q3(pq) - PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q1(pq));

	//	   (q(2)*q(4)-q(1)*q(3))                    (q(3)*q(4)+q(1)*q(2))                   (q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2)/2];
	PDCM_R31(pR) = 2.0*(PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q3(pq) - PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q2(pq));
	PDCM_R32(pR) = 2.0*(PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q3(pq) + PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q1(pq));
	PDCM_R33(pR) = 1.0*(PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq) - PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq) - PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq) + PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq));

	// Retorna 
    return 1; 
}                      

int rotation_dcm2quaternions(PDCM pR, PQUATERNIONS pq)
{
	double a;
//	% improved code: Antonio Padilha
//	q0 = .5*sqrt(max(0,1+R(1,1)+R(2,2)+R(3,3)));
	a = 1.0 + PGMATRIX_DATA(pR,1,1) + PGMATRIX_DATA(pR,2,2) + PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q0(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q0(pq) = 0.0;

//	q1 = .5*sqrt(max(0,1+R(1,1)-R(2,2)-R(3,3)));
	a = 1.0 + PGMATRIX_DATA(pR,1,1) - PGMATRIX_DATA(pR,2,2) - PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q1(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q1(pq) = 0.0;

//	q2 = .5*sqrt(max(0,1-R(1,1)+R(2,2)-R(3,3)));
	a = 1.0 - PGMATRIX_DATA(pR,1,1) + PGMATRIX_DATA(pR,2,2) - PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q2(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q2(pq) = 0.0;

//	q3 = .5*sqrt(max(0,1-R(1,1)-R(2,2)+R(3,3)));
	a = 1.0 - PGMATRIX_DATA(pR,1,1) - PGMATRIX_DATA(pR,2,2) + PGMATRIX_DATA(pR,3,3);
	if (a > 0) PQUATERNIONS_Q3(pq) = 0.5*sqrt(a);
	else PQUATERNIONS_Q3(pq) = 0.0;

//	q1 = abs(q1)*sign(R(3,2)-R(2,3));
	if ((PGMATRIX_DATA(pR,3,2)-PGMATRIX_DATA(pR,2,3)) >= 0.0)
		PQUATERNIONS_Q1(pq) =  fabs(PQUATERNIONS_Q1(pq));
	else
		PQUATERNIONS_Q1(pq) = -fabs(PQUATERNIONS_Q1(pq));

//	q2 = abs(q2)*sign(R(1,3)-R(3,1));
	if ((PGMATRIX_DATA(pR,1,3)-PGMATRIX_DATA(pR,3,1)) >= 0.0)
		PQUATERNIONS_Q2(pq) =  fabs(PQUATERNIONS_Q2(pq));
	else
		PQUATERNIONS_Q2(pq) = -fabs(PQUATERNIONS_Q2(pq));

//	q3 = abs(q3)*sign(R(2,1)-R(1,2));
	if ((PGMATRIX_DATA(pR,2,1)-PGMATRIX_DATA(pR,1,2)) >= 0.0)
		PQUATERNIONS_Q3(pq) =  fabs(PQUATERNIONS_Q3(pq));
	else
		PQUATERNIONS_Q3(pq) = -fabs(PQUATERNIONS_Q3(pq));

	// Retorna 
    return 1; 
}             

int rotation_quaternionsnormalize(PQUATERNIONS pq)
{
	double norm;
	DCM_DECLARE(R);

	// % make quaternion unit length:
	// q = q / norm(q);
	norm  = PQUATERNIONS_Q0(pq)*PQUATERNIONS_Q0(pq);
	norm += PQUATERNIONS_Q1(pq)*PQUATERNIONS_Q1(pq);
	norm += PQUATERNIONS_Q2(pq)*PQUATERNIONS_Q2(pq);
	norm += PQUATERNIONS_Q3(pq)*PQUATERNIONS_Q3(pq);
	norm =  sqrt(norm);

	if(norm == 0.0) return 0;

	PQUATERNIONS_Q0(pq) = PQUATERNIONS_Q0(pq) / norm;
	PQUATERNIONS_Q1(pq) = PQUATERNIONS_Q1(pq) / norm;
	PQUATERNIONS_Q2(pq) = PQUATERNIONS_Q2(pq) / norm;
	PQUATERNIONS_Q3(pq) = PQUATERNIONS_Q3(pq) / norm;
	
	// Retorna 
    return 1; 
}                      

int rotation_dcm2euler(PDCM pR, double *proll, double *ppitch, double *pyaw)
{
	double cosPitch;

	// pitch = asin(-R(3,1));
	*ppitch = asin(-PGMATRIX_DATA(pR,3,1));
	
	cosPitch = cos(*ppitch);
	if (cosPitch == 0)
		return 0;
	
	//roll = atan2(R(3,2)/cos(pitch),R(3,3)/cos(pitch));
	*proll  = atan2(PGMATRIX_DATA(pR,3,2)/cosPitch,PGMATRIX_DATA(pR,3,3)/cosPitch);
	
	// yaw = atan2(R(2,1)/cos(pitch),R(1,1)/cos(pitch));
	*pyaw  = atan2(PGMATRIX_DATA(pR,2,1)/cosPitch,PGMATRIX_DATA(pR,1,1)/cosPitch);
	
	// Retorna 
    return 1; 
}                      

int rotation_euler2dcm(double roll, double pitch, double yaw, PDCM pR)
{

	if (pR->Nr !=3){
		return 0;
	}
	if (pR->Nc !=3){
		return 0;
	}

	// R(1,:) = [cos(pitch)*cos(yaw) -cos(roll)*sin(yaw)+sin(roll)*sin(pitch)*cos(yaw)  sin(roll)*sin(yaw)+cos(roll)*sin(pitch)*cos(yaw)];
	PGMATRIX_DATA(pR,1,1) = cos(pitch)*cos(yaw);
	PGMATRIX_DATA(pR,1,2) = -cos(roll)*sin(yaw)+sin(roll)*sin(pitch)*cos(yaw);
	PGMATRIX_DATA(pR,1,3) = sin(roll)*sin(yaw)+cos(roll)*sin(pitch)*cos(yaw);

	// R(2,:) = [cos(pitch)*sin(yaw)  cos(roll)*cos(yaw)+sin(roll)*sin(pitch)*sin(yaw) -sin(roll)*cos(yaw)+cos(roll)*sin(pitch)*sin(yaw)];
	PGMATRIX_DATA(pR,2,1) = cos(pitch)*sin(yaw);
	PGMATRIX_DATA(pR,2,2) = cos(roll)*cos(yaw)+sin(roll)*sin(pitch)*sin(yaw); 
	PGMATRIX_DATA(pR,2,3) = -sin(roll)*cos(yaw)+cos(roll)*sin(pitch)*sin(yaw);
	
	// R(3,:) = [-sin(pitch)          sin(roll)*cos(pitch)                              cos(roll)*cos(pitch)];
	PGMATRIX_DATA(pR,3,1) = -sin(pitch);
	PGMATRIX_DATA(pR,3,2) = sin(roll)*cos(pitch);
	PGMATRIX_DATA(pR,3,3) = cos(roll)*cos(pitch);

	// Retorna 
    return 1; 
}                      
