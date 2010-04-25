/*
 *  CQRMTest.c
 *  
 *
 *  Created by Herbert J. Bernstein on 2/20/09.
 *  Copyright 2009 Herbert J. Bernstein. All rights reserved.
 *
 */

/*  Work supported in part by NIH NIGMS under grant 1R15GM078077-01 and DOE 
 under grant ER63601-1021466-0009501.  Any opinions, findings, and 
 conclusions or recommendations expressed in this material are those of the   
 author(s) and do not necessarily reflect the views of the funding agencies.
 */


/**********************************************************************
 *                                                                    *
 * YOU MAY REDISTRIBUTE THE CQRlib API UNDER THE TERMS OF THE LGPL    *
 *                                                                    *
 **********************************************************************/

/************************* LGPL NOTICES *******************************
 *                                                                    *
 * This library is free software; you can redistribute it and/or      *
 * modify it under the terms of the GNU Lesser General Public         *
 * License as published by the Free Software Foundation; either       *
 * version 2.1 of the License, or (at your option) any later version. *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
 * Lesser General Public License for more details.                    *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License along with this library; if not, write to the Free         *
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
 * MA  02110-1301  USA                                                *
 *                                                                    *
 **********************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#ifndef USE_LOCAL_HEADERS
#include <cqrlib.h>
#else
#include "cqrlib.h" 
#endif

#ifdef USE_MINGW_RAND
#define rand(x) random(x)
#define srand(x) srandom(x)
#endif


int main(int argc, char ** argv) {
    
    CQRQuaternionHandle q1, q2, q3;
    CQRQuaternion q4, qx, qy, qz;
    double normsq;
    double Matx[3][3], Maty[3][3], Matz[3][3];
    double PI;
    double vx[3] = {1.,0.,0.};
    double vy[3] = {0.,1.,0.};
    double vz[3] = {0.,0.,1.};
    double EXX,EXY,EXZ;
    double EYX,EYY,EYZ;
    double EZX,EZY,EZZ;
    
    PI = 4.*atan2(1.,1.);
    
    if (CQRCreateEmptyQuaternion(&q1)) {
        fprintf(stderr," CQRCreatEmptyQuaternion for q1 failed\n");
    }
    
    if (q1->w !=0.  ||  q1->x !=0. || q1->y !=0. || q1->z !=0.) {
        fprintf(stderr," CQRCreateEmptyQuaternion for q1 non-zero [ %g, %g, %g, %g ]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRCreateQuaternion(&q2, 1.,2.,3.,4.)) {
        fprintf(stderr," CQRCreateQuaternion for q2 failed\n");
    }
    
    if (q2->w !=1.  ||  q2->x !=2. || q2->y !=3. || q2->z !=4.) {
        fprintf(stderr," CQRCreateQuaternion for q2 wrong value [ %g, %g, %g, %g ] != [1.,2,3.,4.]\n",q2->w,q2->x,q2->y,q2->z);
    }
    
    if (CQRCreateEmptyQuaternion(&q3)) {
        fprintf(stderr," CQRCreateEmptyQuaternion for q3 failed\n");
    }
    
    if (CQRAdd(q3,q1,q2)) {
        fprintf(stderr," CQRAdd(q3,q1,q2) failed\n");
    }
    
    if (q3->w !=1.  ||  q3->x !=2. || q3->y !=3. || q3->z !=4.) {
        fprintf(stderr," CQRAdd(q3,q1,q2) q3 wrong value [ %g, %g, %g, %g ] != [1.,2.,3.,4.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQRSetQuaternion(q1,-9999.,-9998.,-9997.,-9996.)) {
        fprintf(stderr," CQRSetQuaternion(q1,-9999.,-9998.,-9997.,-9996.)\n");
    }
    
    if (q1->w !=-9999.  ||  q1->x !=-9998. || q1->y !=-9997. || q1->z !=-9996.) {
        fprintf(stderr," CQRSetQuaternion q1 wrong value [ %g, %g, %g, %g ] != [-9999.,-9998.,-9997.,-9996.]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRSubtract(q1,q3,q2)) {
        fprintf(stderr," CQRSubtract(q1,q3,q2) failed\n");
    }
    
    if (q1->w !=0.  ||  q1->x !=0. || q1->y !=0. || q1->z !=0.) {
        fprintf(stderr," CQR Subtract(q1,q3,q2) for q1 non-zero [ %g, %g, %g, %g ]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRSetQuaternion(&q4,1.,-2.,-3.,-4.)) {
        fprintf(stderr," CQRSetQuaternion(&q4,1.,-2.,-3.,-4.)\n");
    }
    
    if (q4.w != 1.  ||  q4.x !=-2. || q4.y !=-3. || q4.z !=-4.) {
        fprintf(stderr," CQRSetQuaternion &q4 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q4.w,q4.x,q4.y,q4.z);
    }
    
    if (CQRMultiply(q1,&q4,q2)) {
        fprintf(stderr," CQRMultiply(q1,&q4,q2) failed\n");
    }
    
    if (fabs(q1->w-30.)>DBL_MIN || fabs(q1->x)>DBL_MIN || fabs(q1->y)>DBL_MIN ||fabs(q1->z)>DBL_MIN )  {
        fprintf(stderr,"  CQRMultiply(q1,&q4,q2)  q1 wrong value [ %g, %g, %g, %g ] != [30.,0.,0.,0.]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRDivide(q3,q1,q2)) {
        fprintf(stderr," CQRDivide(q3,q1,q2) failed\n");
    }
    
    if (fabs(q3->w-1.)>DBL_MIN || fabs(q3->x+2.)>DBL_MIN || fabs(q3->y+3.)>DBL_MIN ||fabs(q3->z+4.)>DBL_MIN )  {
        fprintf(stderr,"  CQRDivide(q3,q1,q2)  q3 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQRScalarMultiply(q3,q2,3.)) {
        fprintf(stderr," CQRScalarMultiply(q3,q1,3.) failed\n");
    }
    
    if (fabs(q3->w-3.)>DBL_MIN || fabs(q3->x-6.)>DBL_MIN || fabs(q3->y-9.)>DBL_MIN ||fabs(q3->z-12.)>DBL_MIN )  {
        fprintf(stderr,"  CQRScalarMultiply(q3,q2,3.)  q3 wrong value [ %g, %g, %g, %g ] != [3.,6.,9.,12.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQRConjugate(q3,q2)) {
        fprintf(stderr," CQRConjugate(q3,q2) failed\n");
    }
    
    if (q3->w != 1.  ||  q3->x !=-2. || q3->y !=-3. || q3->z !=-4.) {
        fprintf(stderr," CQRConjugate(q3,q2) q3 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQREqual(&q4,q3)) {
        fprintf(stderr," CQREqual(&q4,q2) failed\n");
    }
    
    if (CQRNormsq(&normsq,&q4) || normsq !=30.) {
        fprintf(stderr," CQRNormsq(&normsq,&q4) failed\n");
    }
    
    if (CQRInverse(q3,&q4)) {
        fprintf(stderr,"CQRInverse(q3,&q4) failed\n");
    }
    
    if (fabs(q3->w - 1./30.) > DBL_MIN  ||  fabs(q3->x - 2./30.) > DBL_MIN || fabs(q3->y - 3./30.) > DBL_MIN || fabs(q3->z - 4./30.) > DBL_MIN) {
        fprintf(stderr," CQRInverse(q3,&q4) q3 wrong value [ %g, %g, %g, %g ] != [1./30.,2./30.,3./30,4./30.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    
    /* Create quaternions to rotate about the x,y and z-axes by 90 degrees */
    
    if (CQRAxis2Quaternion(&qx,vx,PI/2.)||CQRAxis2Quaternion(&qy,vy,PI/2.)||CQRAxis2Quaternion(&qz,vz,PI/2.)){
        fprintf(stderr,"Axis2Quaternion failed\n");
    }
    
    if (qx.w<0.||fabs(qx.w*qx.w-.5)>1.e-13||fabs(qx.x*qx.x-.5)>1.e-13||fabs(qx.y)>DBL_MIN||fabs(qx.z)>DBL_MIN) {
        fprintf(stderr,"Axis2Quaternion qx wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qx.w,qx.x,qx.y,qx.z);
    }
    
    
    if (qy.w<0.||fabs(qy.w*qy.w-.5)>1.e-13||fabs(qy.y*qy.y-.5)>1.e-13||fabs(qy.x)>DBL_MIN||fabs(qy.z)>DBL_MIN) {
        fprintf(stderr,"Axis2Quaternion qy wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qy.w,qy.x,qy.y,qy.z);
    }
    
    
    if (qz.w<0.||fabs(qz.w*qz.w-.5)>1.e-13||fabs(qz.z*qz.z-.5)>1.e-13||fabs(qz.x)>DBL_MIN||fabs(qz.y)>DBL_MIN) {
        fprintf(stderr,"Axis2Quaternion qz wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qz.w,qz.x,qz.y,qz.z);
    }
    
    if (CQRQuaternion2Matrix(Matx,&qx)||CQRQuaternion2Matrix(Maty,&qy)||CQRQuaternion2Matrix(Matz,&qz)){
        fprintf(stderr," CQRQuaternion2Matrix failed\n");
    }
    
    if (fabs(Matx[0][0]-1.)>1.e-13 ||fabs(Matx[1][2]+1.)>1.e-13||fabs(Matx[2][1]-1.)>1.e-13
        ||fabs(Matx[0][1])>1.e-13 ||fabs(Matx[0][2])>1.e-13
        ||fabs(Matx[1][0])>1.e-13 ||fabs(Matx[1][1])>1.e-13
        ||fabs(Matx[2][0])>1.e-13 ||fabs(Matx[2][2])>1.e-13) {
        fprintf(stderr," CQRQuaternion2Matrix Matx wrong value \n  [ %g, %g, %g ]\n  [ %g, %g, %g ]\n  [ %g, %g, %g ]\n"
                "!=  [1, 0, 0]\n    [0, 0, -1]\n    [0, 1, 0]\n",
                Matx[0][0],Matx[0][1],Matx[0][2],
                Matx[1][0],Matx[1][1],Matx[1][2],
                Matx[2][0],Matx[2][1],Matx[2][2]);
    }
    
    if (fabs(Maty[0][2]-1.)>1.e-13 ||fabs(Maty[1][1]-1.)>1.e-13||fabs(Maty[2][0]+1.)>1.e-13
        ||fabs(Maty[0][0])>1.e-13 ||fabs(Maty[0][1])>1.e-13
        ||fabs(Maty[1][0])>1.e-13 ||fabs(Maty[1][2])>1.e-13
        ||fabs(Maty[2][1])>1.e-13 ||fabs(Maty[2][2])>1.e-13) {
        fprintf(stderr," CQRQuaternion2Matrix Maty wrong value \n  [ %g, %g, %g]\n  [ %g, %g, %g]\n  [ %g, %g, %g ]\n"
                "!=  [0, 0, 1]\n    [0, 1, 0]\n    [-1, 0, 0]\n",
                Maty[0][0],Maty[0][1],Maty[0][2],
                Maty[1][0],Maty[1][1],Maty[1][2],
                Maty[2][0],Maty[2][1],Maty[2][2]);
    }
    
    if (fabs(Matz[0][1]+1.)>1.e-13 ||fabs(Matz[1][0]-1.)>1.e-13||fabs(Matz[2][2]-1.)>1.e-13
        ||fabs(Matz[0][0])>1.e-13 ||fabs(Matz[0][2])>1.e-13
        ||fabs(Matz[1][1])>1.e-13 ||fabs(Matz[1][2])>1.e-13
        ||fabs(Matz[2][0])>1.e-13 ||fabs(Matz[2][1])>1.e-13) {
        fprintf(stderr," CQRQuaternion2Matrix Matz wrong value \n  [ %g, %g, %g]\n  [ %g, %g, %g]\n  [ %g, %g, %g ]\n"
                "!=  [0, -1, 0]\n    [1, 0, 0]\n    [0, 0, 1]\n",
                Matz[0][0],Matz[0][1],Matz[0][2],
                Matz[1][0],Matz[1][1],Matz[1][2],
                Matz[2][0],Matz[2][1],Matz[2][2]);
    }
    
    
    EXX = EXY = EXZ = 0.;
    EYX = EYY = EYZ = 0.;
    EZX = EZY = EZZ = 0.;
    if (CQRQuaternion2Angles(&EXX,&EXY,&EXZ,&qx)
        ||CQRQuaternion2Angles(&EYX,&EYY,&EYZ,&qy)
        ||CQRQuaternion2Angles(&EZX,&EZY,&EZZ,&qz) ){
        fprintf(stderr," CQRQuaternion2Angles failed\n");
    }
    
    if (CQRAngles2Quaternion(q1,EXX,EXY,EXZ)||CQRAngles2Quaternion(q2,EYX,EYY,EYZ)||CQRAngles2Quaternion(q3,EZX,EZY,EZZ)){
        fprintf(stderr," CQRAngles2Quaternion failed\n");
    }
    
    if (CQRDivide(&q4,&qx,q1) || CQRNormsq(&normsq,&q4) || fabs(normsq-1.) > 1.e-13 || fabs(q4.w*q4.w-1.) > 1.e-13) {
        fprintf(stderr,"  CQRAngles2Quaternion q1 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qx.w, qx.x, qx.y, qx.z, q1->w, q1->x, q1->y, q1->z );  
    }
    
    if (CQRDivide(&q4,&qy,q2) || CQRNormsq(&normsq,&q4) || fabs(normsq-1.) > 1.e-13 || fabs(q4.w*q4.w-1.) > 1.e-13) {
        fprintf(stderr,"  CQRAngles2Quaternion q2 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qy.w, qy.x, qy.y, qy.z, q2->w, q2->x, q2->y, q2->z );  
    }
    
    
    if (CQRDivide(&q4,&qz,q3) || CQRNormsq(&normsq,&q4) || fabs(normsq-1.) > 1.e-13 || fabs(q4.w*q4.w-1.) > 1.e-13) {
        fprintf(stderr,"  CQRAngles2Quaternion q3 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qz.w, qz.x, qz.y, qz.z, q3->w, q3->x, q3->y, q3->z );  
    }
    
    
    if(CQRFreeQuaternion(&q1)) {
        fprintf(stderr," CQRFreeQuaternion(&q1) failed\n");
    }
    
    if(CQRFreeQuaternion(&q2)) {
        fprintf(stderr," CQRFreeQuaternion(&q2) failed\n");
    }
    if(CQRFreeQuaternion(&q3)) {
        fprintf(stderr," CQRFreeQuaternion(&q3) failed\n");
    }
    
    
    return 0;
}


