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
    int errorcount;
    int i, j;
    
    PI = 4.*atan2(1.,1.);
    errorcount = 0;
    
    if (CQRCreateEmptyQuaternion(&q1)) {
        errorcount++;
        fprintf(stderr," CQRCreatEmptyQuaternion for q1 failed\n");
    }
    
    if (q1->w !=0.  ||  q1->x !=0. || q1->y !=0. || q1->z !=0.) {
        errorcount++;
        fprintf(stderr," CQRCreateEmptyQuaternion for q1 non-zero [ %g, %g, %g, %g ]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRCreateQuaternion(&q2, 1.,2.,3.,4.)) {
        errorcount++;
        fprintf(stderr," CQRCreateQuaternion for q2 failed\n");
    }
    
    if (q2->w !=1.  ||  q2->x !=2. || q2->y !=3. || q2->z !=4.) {
        errorcount++;
        fprintf(stderr," CQRCreateQuaternion for q2 wrong value [ %g, %g, %g, %g ] != [1.,2,3.,4.]\n",q2->w,q2->x,q2->y,q2->z);
    }
    
    if (CQRCreateEmptyQuaternion(&q3)) {
        errorcount++;
        fprintf(stderr," CQRCreateEmptyQuaternion for q3 failed\n");
    }
    
    if (CQRAdd(q3,q1,q2)) {
        errorcount++;
        fprintf(stderr," CQRAdd(q3,q1,q2) failed\n");
    }
    
    if (q3->w !=1.  ||  q3->x !=2. || q3->y !=3. || q3->z !=4.) {
        errorcount++;
        fprintf(stderr," CQRAdd(q3,q1,q2) q3 wrong value [ %g, %g, %g, %g ] != [1.,2.,3.,4.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQRSetQuaternion(q1,-9999.,-9998.,-9997.,-9996.)) {
        errorcount++;
        fprintf(stderr," CQRSetQuaternion(q1,-9999.,-9998.,-9997.,-9996.)\n");
    }
    
    if (q1->w !=-9999.  ||  q1->x !=-9998. || q1->y !=-9997. || q1->z !=-9996.) {
        errorcount++;
        fprintf(stderr," CQRSetQuaternion q1 wrong value [ %g, %g, %g, %g ] != [-9999.,-9998.,-9997.,-9996.]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRSubtract(q1,q3,q2)) {
        errorcount++;
        fprintf(stderr," CQRSubtract(q1,q3,q2) failed\n");
    }
    
    if (q1->w !=0.  ||  q1->x !=0. || q1->y !=0. || q1->z !=0.) {
        errorcount++;
        fprintf(stderr," CQR Subtract(q1,q3,q2) for q1 non-zero [ %g, %g, %g, %g ]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRSetQuaternion(&q4,1.,-2.,-3.,-4.)) {
        errorcount++;
        fprintf(stderr," CQRSetQuaternion(&q4,1.,-2.,-3.,-4.)\n");
    }
    
    if (q4.w != 1.  ||  q4.x !=-2. || q4.y !=-3. || q4.z !=-4.) {
        errorcount++;
        fprintf(stderr," CQRSetQuaternion &q4 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q4.w,q4.x,q4.y,q4.z);
    }
    
    if (CQRMultiply(q1,&q4,q2)) {
        errorcount++;
        fprintf(stderr," CQRMultiply(q1,&q4,q2) failed\n");
    }
    
    if (fabs(q1->w-30.)>300.*DBL_EPSILON || fabs(q1->x)>30.*DBL_EPSILON || fabs(q1->y)>30.*DBL_EPSILON ||fabs(q1->z)>30.*DBL_EPSILON )  {
        errorcount++;
        fprintf(stderr,"  CQRMultiply(q1,&q4,q2)  q1 wrong value [ %g, %g, %g, %g ] != [30.,0.,0.,0.]\n",q1->w,q1->x,q1->y,q1->z);
    }
    
    if (CQRDivide(q3,q1,q2)) {
        errorcount++;
        fprintf(stderr," CQRDivide(q3,q1,q2) failed\n");
    }
    
    if (fabs(q3->w-1.)>60.*DBL_EPSILON || fabs(q3->x+2.)>60.*DBL_EPSILON || fabs(q3->y+3.)>60.*DBL_EPSILON ||fabs(q3->z+4.)>60.*DBL_EPSILON )  {
        errorcount++;
        fprintf(stderr,"  CQRDivide(q3,q1,q2)  q3 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQRScalarMultiply(q3,q2,3.)) {
        errorcount++;
        fprintf(stderr," CQRScalarMultiply(q3,q1,3.) failed\n");
    }
    
    if (fabs(q3->w-3.)>180.*DBL_EPSILON || fabs(q3->x-6.)>180.*DBL_EPSILON || fabs(q3->y-9.)>180.*DBL_EPSILON ||fabs(q3->z-12.)>180.*DBL_EPSILON )  {
        errorcount++;
        fprintf(stderr,"  CQRScalarMultiply(q3,q2,3.)  q3 wrong value [ %g, %g, %g, %g ] != [3.,6.,9.,12.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQRConjugate(q3,q2)) {
        errorcount++;
        fprintf(stderr," CQRConjugate(q3,q2) failed\n");
    }
    
    if (q3->w != 1.  ||  q3->x !=-2. || q3->y !=-3. || q3->z !=-4.) {
        errorcount++;
        fprintf(stderr," CQRConjugate(q3,q2) q3 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    if (CQREqual(&q4,q3)) {
        errorcount++;
        fprintf(stderr," CQREqual(&q4,q2) failed\n");
    }
    
    if (CQRNormsq(&normsq,&q4) || fabs(normsq-30.) > 300.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr," CQRNormsq(&normsq,&q4) failed\n");
    }
    
    if (CQRInverse(q3,&q4)) {
        errorcount++;
        fprintf(stderr,"CQRInverse(q3,&q4) failed\n");
    }
    
    if (fabs(q3->w - 1./30.) > 2.*DBL_EPSILON  ||  fabs(q3->x - 2./30.) > 2.*DBL_EPSILON || fabs(q3->y - 3./30.) > 2.*DBL_EPSILON || fabs(q3->z - 4./30.) > 2.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr," CQRInverse(q3,&q4) q3 wrong value [ %g, %g, %g, %g ] != [1./30.,2./30.,3./30,4./30.]\n",q3->w,q3->x,q3->y,q3->z);
    }
    
    
    /* Create quaternions to rotate about the x,y and z-axes by 90 degrees */
    
    if (CQRAxis2Quaternion(&qx,vx,PI/2.)||CQRAxis2Quaternion(&qy,vy,PI/2.)||CQRAxis2Quaternion(&qz,vz,PI/2.)){
        errorcount++;
        fprintf(stderr,"Axis2Quaternion failed\n");
    }
    
    if (qx.w<0.||fabs(qx.w*qx.w-.5)>10.*DBL_EPSILON||fabs(qx.x*qx.x-.5)>10.*DBL_EPSILON||fabs(qx.y)>10.*DBL_EPSILON||fabs(qx.z)>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr,"Axis2Quaternion qx wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qx.w,qx.x,qx.y,qx.z);
    }
    
    
    if (qy.w<0.||fabs(qy.w*qy.w-.5)>10.*DBL_EPSILON||fabs(qy.y*qy.y-.5)>10.*DBL_EPSILON||fabs(qy.x)>10.*DBL_EPSILON||fabs(qy.z)>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr,"Axis2Quaternion qy wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qy.w,qy.x,qy.y,qy.z);
    }
    
    
    if (qz.w<0.||fabs(qz.w*qz.w-.5)>10.*DBL_EPSILON||fabs(qz.z*qz.z-.5)>10.*DBL_EPSILON||fabs(qz.x)>10.*DBL_EPSILON||fabs(qz.y)>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr,"Axis2Quaternion qz wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qz.w,qz.x,qz.y,qz.z);
    }
    
    if (CQRQuaternion2Matrix(Matx,&qx)||CQRQuaternion2Matrix(Maty,&qy)||CQRQuaternion2Matrix(Matz,&qz)){
        errorcount++;
        fprintf(stderr," CQRQuaternion2Matrix failed\n");
    }
    
    if (fabs(Matx[0][0]-1.)>10.*DBL_EPSILON ||fabs(Matx[1][2]+1.)>10.*DBL_EPSILON||fabs(Matx[2][1]-1.)>10.*DBL_EPSILON
        ||fabs(Matx[0][1])>10.*DBL_EPSILON ||fabs(Matx[0][2])>10.*DBL_EPSILON
        ||fabs(Matx[1][0])>10.*DBL_EPSILON ||fabs(Matx[1][1])>10.*DBL_EPSILON
        ||fabs(Matx[2][0])>10.*DBL_EPSILON ||fabs(Matx[2][2])>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr," CQRQuaternion2Matrix Matx wrong value \n  [ %g, %g, %g ]\n  [ %g, %g, %g ]\n  [ %g, %g, %g ]\n"
                "!=  [1, 0, 0]\n    [0, 0, -1]\n    [0, 1, 0]\n",
                Matx[0][0],Matx[0][1],Matx[0][2],
                Matx[1][0],Matx[1][1],Matx[1][2],
                Matx[2][0],Matx[2][1],Matx[2][2]);
    }
    
    if (fabs(Maty[0][2]-1.)>DBL_EPSILON ||fabs(Maty[1][1]-1.)>DBL_EPSILON||fabs(Maty[2][0]+1.)>DBL_EPSILON
        ||fabs(Maty[0][0])>DBL_EPSILON ||fabs(Maty[0][1])>DBL_EPSILON
        ||fabs(Maty[1][0])>DBL_EPSILON ||fabs(Maty[1][2])>DBL_EPSILON
        ||fabs(Maty[2][1])>DBL_EPSILON ||fabs(Maty[2][2])>DBL_EPSILON) {
        errorcount++;
        fprintf(stderr," CQRQuaternion2Matrix Maty wrong value \n  [ %g, %g, %g]\n  [ %g, %g, %g]\n  [ %g, %g, %g ]\n"
                "!=  [0, 0, 1]\n    [0, 1, 0]\n    [-1, 0, 0]\n",
                Maty[0][0],Maty[0][1],Maty[0][2],
                Maty[1][0],Maty[1][1],Maty[1][2],
                Maty[2][0],Maty[2][1],Maty[2][2]);
    }
    
    if (fabs(Matz[0][1]+1.)>DBL_EPSILON ||fabs(Matz[1][0]-1.)>DBL_EPSILON||fabs(Matz[2][2]-1.)>DBL_EPSILON
        ||fabs(Matz[0][0])>DBL_EPSILON ||fabs(Matz[0][2])>DBL_EPSILON
        ||fabs(Matz[1][1])>DBL_EPSILON ||fabs(Matz[1][2])>DBL_EPSILON
        ||fabs(Matz[2][0])>DBL_EPSILON ||fabs(Matz[2][1])>DBL_EPSILON) {
        errorcount++;
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
        errorcount++;
        fprintf(stderr," CQRQuaternion2Angles failed\n");
    }
    
    if (CQRAngles2Quaternion(q1,EXX,EXY,EXZ)||CQRAngles2Quaternion(q2,EYX,EYY,EYZ)||CQRAngles2Quaternion(q3,EZX,EZY,EZZ)){
        errorcount++;
        fprintf(stderr," CQRAngles2Quaternion failed\n");
    }
    
    if (CQRDivide(&q4,&qx,q1) || CQRNormsq(&normsq,&q4) || fabs(normsq-1.) > 10.*DBL_EPSILON || fabs(q4.w*q4.w-1.) > 10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr,"  CQRAngles2Quaternion q1 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qx.w, qx.x, qx.y, qx.z, q1->w, q1->x, q1->y, q1->z );  
    }
    
    if (CQRDivide(&q4,&qy,q2) || CQRNormsq(&normsq,&q4) || fabs(normsq-1.) > 10.*DBL_EPSILON || fabs(q4.w*q4.w-1.) > 10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr,"  CQRAngles2Quaternion q2 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qy.w, qy.x, qy.y, qy.z, q2->w, q2->x, q2->y, q2->z );  
    }
    
    
    if (CQRDivide(&q4,&qz,q3) || CQRNormsq(&normsq,&q4) || fabs(normsq-1.) > 10.*DBL_EPSILON || fabs(q4.w*q4.w-1.) > 10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stderr,"  CQRAngles2Quaternion q3 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qz.w, qz.x, qz.y, qz.z, q3->w, q3->x, q3->y, q3->z );  
    }
    
    
    if(CQRFreeQuaternion(&q1)) {
        errorcount++;
        fprintf(stderr," CQRFreeQuaternion(&q1) failed\n");
    }
    
    if(CQRFreeQuaternion(&q2)) {
        errorcount++;
        fprintf(stderr," CQRFreeQuaternion(&q2) failed\n");
    }

    if(CQRFreeQuaternion(&q3)) {
        errorcount++;
        fprintf(stderr," CQRFreeQuaternion(&q3) failed\n");
    }
    
    
    {
        CQRQuaternion q1 = {3.,5.,7.,9.};
        double qw,qx,qy,qz;
        
        CQRGetQuaternionW(&qw,&q1);
        CQRGetQuaternionX(&qx,&q1);
        CQRGetQuaternionY(&qy,&q1);
        CQRGetQuaternionZ(&qz,&q1);
        
        if ( qw != q1.w || qx != q1.x || qy != q1.y || qz != q1.z )
        {
            errorcount++;
            fprintf( stdout, "CQRGetQuaternionW/X/Y/Z failed\n" );
        }
        
        
    }
    
    /*  Tests on [-sqrt(7),2,3,4] = 6*[-sqrt(7)/6,1/3,1/2,2/3]  
     = 6*[-cos(1.11412994158827),sin(1.11412994158827)*[.3713906763541037, .5570860145311556, .7427813527082074]]
     = 6*[cos(2.027462712001523),sin(2.027462712001523)*[.3713906763541037, .5570860145311556, .7427813527082074]]
     = 6*exp([0,.3713906763541037, .5570860145311556, .7427813527082074]*2.027462712001523)
     
     so the log should be
     
     [log(6),0,0,0] +[0,.3713906763541037, .5570860145311556, .7427813527082074]*2.027462712001523]
     =[1.791759469228055, 0.752980747892971, 1.129471121839456, 1.505961495785942]
     
     Note that the log is multivalued
     
     */
    {
        CQRQuaternion q1, q1Im, q1axis, q1log, qtemp, q1logexp, q1exp, q1explog, q1explogexp, q1powi, q1powd;
        double q1angle, norm, q1lognorm;
        double norm1, norm2, norm3, norm4;
        
        CQRMSet(q1,-sqrt(7.),2.,3.,4.);
        CQRGetQuaternionIm(&q1Im,&q1);
        
        if (q1Im.w != 0. || q1Im.x != 2. || q1Im.y != 3. || q1Im.z != 4. )
        {
            errorcount++;
            fprintf( stdout, "CQRGetQuaternionIm failed\n" );
        }
        
        CQRGetQuaternionAxis(&q1axis,&q1);
        
        if (q1axis.w != 0. || q1axis.x != 2./sqrt(4.+9.+16.) 
                           || q1axis.y != 3./sqrt(4.+9.+16.) 
                           || q1axis.z != 4./sqrt(4.+9.+16.) )
        {
            errorcount++;
            fprintf( stdout, "CQRGetQuaternionAxis failed\n" );
        }

        CQRGetQuaternionAngle(&q1angle,&q1);
        
        if (fabs(q1angle-2.027462712001523)>10.*DBL_EPSILON*2.027462712001523)
        {
            errorcount++;
            fprintf( stdout, "CQRGetQuaternionAngle failed, got %g, expected %g\n",q1angle,2.027462712001523 );
        }
        
        CQRLog(&q1log,&q1);
        CQRMSet(qtemp,log(6.), 0.752980747892971, 1.129471121839457, 1.505961495785942)
        CQRMSubtract(qtemp,qtemp,q1log);
        CQRMNorm(norm,qtemp);
        CQRMNorm(q1lognorm,q1log)
        if (norm >  10.*DBL_EPSILON*q1lognorm)
        {
            errorcount++;
            fprintf( stdout, "quaternion log failed log([%g,%g,%g,%g]) = [%g,%g,%g,%g] instead of [%g,%g,%g,%g], normdiff = %g\n",
                    q1.w, q1.x, q1.y, q1.z,
                    q1log.w, q1log.x, q1log.y, q1log.z,
                    log(6.), 0.752980747892971, 1.129471121839457, 1.505961495785942,
                    norm);
            
        }
        
        CQRLog(&q1log,&q1);
        CQRExp(&q1logexp,&q1log);
        CQRExp(&q1exp,&q1);
        CQRLog(&q1explog,&q1exp);
        CQRExp(&q1explogexp,&q1explog);
        CQRMSubtract(qtemp,q1logexp,q1); CQRMNorm(norm1,qtemp);
        CQRMSubtract(qtemp,q1explogexp,q1exp); CQRMNorm(norm2,qtemp);
        CQRMNorm(norm3,q1);
        CQRMNorm(norm4,q1exp)
        
        if (norm1>10.*DBL_EPSILON*norm3 || 
            norm2>10.*DBL_EPSILON*norm4)
        {
            errorcount++;
            fprintf( stdout, "log(exp) or exp(log) failed\n," 
                    " q = [%g,%g,%g,%g],"
                    " log = [%g,%g,%g,%g], exp(log) = [%g,%g,%g,%g],"
                    " exp = [%g,%g,%g,%g], log(exp) = [%g,%g,%g,%g]\n",
                    q1.w,q1.x,q1.y,q1.z,
                    q1log.w,q1log.x,q1log.y,q1log.z,
                    q1logexp.w,q1logexp.x,q1logexp.y,q1logexp.z,
                    q1exp.w,q1exp.x,q1exp.y,q1exp.z,
                    q1explogexp.w,q1explogexp.x,q1explogexp.y,q1explogexp.z
                    );
        }
        
        for (i = -5; i < 6; i++) {
             CQRDoublePower(&q1powd,&q1,(double)i);
             CQRIntegerPower(&q1powi,&q1,i);
             CQRMSubtract(qtemp,q1powd,q1powi);
             CQRMNorm(norm,qtemp);
             CQRMNorm(norm1,q1powi);
            if (norm > 10.*DBL_EPSILON*norm1)
            {
                errorcount++;
                fprintf( stdout, "integer power double power comparison failed\n,"); 
                
            }
        }
        
        
    }
    
    {
        CQRQuaternion q1, q2, q3, qout1, qout2, qout3, qtest1, qtest2, qtest3;
        double norm1, norm2, norm3;
        double normq1, normq2, normq3;
        
        CQRMSet (q1, -4.,0.,0.,0. );
        CQRMSet (q2, -4.,1.,1.,1. );
        CQRMSet (q3,  4.,0.,0.,0. );
        
        for (i = 1; i < 9; i++) {
            for (j = 0; j <  i; j++ ) {
                CQRIntegerRoot(&qout1,&q1,i,j);
                CQRIntegerRoot(&qout2,&q2,i,j);
                CQRIntegerRoot(&qout3,&q3,i,j);
                CQRIntegerPower(&qtest1,&qout1,i);
                CQRIntegerPower(&qtest2,&qout2,i);
                CQRIntegerPower(&qtest3,&qout3,i);
                CQRMNorm(normq1,q1);
                CQRMNorm(normq2,q2);
                CQRMNorm(normq3,q3);
                CQRMDist(norm1,q1,qtest1);
                CQRMDist(norm2,q2,qtest2);
                CQRMDist(norm3,q3,qtest3);
                if (norm1 > 100.*DBL_EPSILON*normq1
                    || norm2 > 100.*DBL_EPSILON*normq2
                    || norm3 > 100.*DBL_EPSILON*normq3) {
                    errorcount++;
                    fprintf(stdout," %d'th root of [%g,%g,%g,%g] = [%g,%g,%g,%g], power = [%g,%g,%g,%g], delta %g\n",
                            i, q1.w, q1.x, q1.y, q1.z,
                            qout1.w, qout1.x, qout1.y, qout1.z,
                            qtest1.w, qtest1.x, qtest1.y, qtest1.z,
                            norm1);
                    fprintf(stdout," %d'th root of [%g,%g,%g,%g] = [%g,%g,%g,%g], power = [%g,%g,%g,%g], delta %g\n",
                            i, q2.w, q2.x, q2.y, q2.z,
                            qout2.w, qout2.x, qout2.y, qout2.z,
                            qtest2.w, qtest2.x, qtest2.y, qtest2.z,
                            norm2);
                    fprintf(stdout," %d'th root of [%g,%g,%g,%g] = [%g,%g,%g,%g], power = [%g,%g,%g,%g], delta %g\n",
                            i, q3.w, q3.x, q3.y, q3.z,
                            qout3.w, qout3.x, qout3.y, qout3.z,
                            qtest3.w, qtest3.x, qtest3.y, qtest3.z,
                            norm3);
                    
                }
            }
        }
        
    }
    
    
    
    
    
    return errorcount;
}


