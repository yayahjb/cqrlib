/*
 *  cqrlib.h
 *  
 *
 *  Created by Herbert J. Bernstein on 2/15/09.
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

/* A utility library for quaternion arithmetic and
 quaternion rotation math.  See
 
 "Quaternions and spatial rotation", Wikipedia
 http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
 
 K. Shoemake, "Quaternions", Department of Computer Science,
 University of Pennsylvania, Philadelphia, PA 19104,
 ftp://ftp.cis.upenn.edu/pub/graphics/shoemake/quatut.ps.Z
 
 K. Shoemake, "Animating rotation with quaternion curves",
 ACM SIGGRAPH Computer Graphics, Vol 19, No. 3, pp 245--254,
 1985.
 
 */

#ifndef CQRLIB_H_INCLUDED
#define CQRLIB_H_INCLUDED

#ifdef __cplusplus

extern "C" {
    
#endif
    
#ifdef CNEARTREE_USE_FAR
#include <malloc.h>
#define FAR __far
#define MALLOC _fmalloc
#define FREE _ffree
#define MEMSET _fmemset
#define MEMMOVE _fmemmove
#else
#include <stdlib.h>
#define FAR
#define MALLOC malloc
#define FREE free
#define MEMSET memset
#define MEMMOVE memmove
#endif
#include <math.h>
    
#define CQR_FAILED       4
#define CQR_NO_MEMORY    2
#define CQR_BAD_ARGUMENT 1
#define CQR_SUCCESS      0
    
    typedef struct {
        double w;
        double x;
        double y;
        double z; } CQRQuaternion;
    
    typedef CQRQuaternion * CQRQuaternionHandle;
    
    /* CQR Macros */
    
#define CQRMCopy(copy,orig) \
(copy).w = (orig).w; (copy).x = (orig).x; (copy).y = (orig).y; (copy).z = (orig).z;

#define CQRMSet(q,qw,qx,qy,qz) \
(q).w = (qw); (q).x = (qx); (q).y = (qy); (q).z = (qz);

#define CQRMAdd(sum,q1,q2) \
(sum).w = (q1).w + (q2).w; sum.x = (q1).x + (q2).x; sum.y = (q1).y + (q2).y; sum.z = (q1).z + (q2).z;

#define CQRMSubtract(sum,q1,q2) \
(sum).w = (q1).w - (q2).w; sum.x = (q1).x - (q2).x; sum.y = (q1).y - (q2).y; sum.z = (q1).z - (q2).z;

#define CQRMMultiply(product,q1,q2 ) \
(product).w = -(q1).z*(q2).z - (q1).y*(q2).y - (q1).x*(q2).x + (q1).w*(q2).w; \
(product).x =  (q1).y*(q2).z - (q1).z*(q2).y + (q1).w*(q2).x + (q1).x*(q2).w; \
(product).y = -(q1).x*(q2).z + (q1).w*(q2).y + (q1).z*(q2).x + (q1).y*(q2).w; \
(product).z =  (q1).w*(q2).z + (q1).x*(q2).y - (q1).y*(q2).x + (q1).z*(q2).w;

#define CQRMScalarMultiply(product,q,s ) \
(product).w = (q).w*s; \
(product).x = (q).x*s; \
(product).y = (q).y*s; \
(product).z = (q).z*s;

#define CQRMConjugate(conjugate,q ) \
(conjugate).w = (q).w; \
(conjugate).x = -(q).x; \
(conjugate).y = -(q).y; \
(conjugate).z = -(q).z;

#define CQRMNormsq(normsq,q) \
normsq = (q).w*(q).w + (q).x*(q).x + (q).y*(q).y + (q).z*(q).z;

#define CQRMInverse(inverseq,q) \
{ double normsq; \
CQRMConjugate(inverseq,q); \
CQRMNormsq(normsq,q); \
if (normsq > 0.) { \
CQRMScalarMultiply(inverserq,1./normsq); \
}
    
    
    /* CQRCreateQuaternion -- create a quaternion = w +ix+jy+kz */
    
    int CQRCreateQuaternion(CQRQuaternionHandle FAR * quaternion, double w, double x, double y, double z); 
    
    /* CQRCreateEmptyQuaternion -- create a quaternion = 0 +i0+j0+k0 */
    
    int CQRCreateEmptyQuaternion(CQRQuaternionHandle FAR * quaternion) ;
    
    /* CQRFreeQuaternion -- free a quaternion  */
    
    int CQRFreeQuaternion(CQRQuaternionHandle FAR * quaternion);        
    
    /* CQRSetQuaternion -- create an existing quaternion = w +ix+jy+kz */
    
    int CQRSetQuaternion( CQRQuaternionHandle quaternion, double w, double x, double y, double z);

    /*  CQRAdd -- add a quaternion (q1) to a quaternion (q2) */
    
    int CQRAdd (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 );
    
    /*  CQRSubtract -- subtract a quaternion (q2) from a quaternion (q1)  */
    
    int CQRSubtract (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 );
    
    /*  CQRMultiply -- multiply a quaternion (q1) by quaternion (q2)  */
    
    int CQRMultiply (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 );
    
    /*  CQRDivide -- Divide a quaternion (q1) by quaternion (q2)  */
    
    int CQRDivide (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 );

    /*  CQRScalarMultiply -- multiply a quaternion (q) by scalar (s)  */
    
    int CQRScalarMultiply (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q, double s );

    /*  CQREqual -- return 0 if quaternion q1 == q2  */
    
    int CQREqual (CQRQuaternionHandle q1, CQRQuaternionHandle q2 );
    
    /*  CQRConjugate -- Form the conjugate of a quaternion qconj */

    int CQRConjugate (CQRQuaternionHandle qconjgate, CQRQuaternionHandle quaternion);
    
    /*  CQRNormsq -- Form the normsquared of a quaternion */
    
    int CQRNormsq (double * normsq, CQRQuaternionHandle quaternion ) ;
    
    /*  CQRInverse -- Form the inverse of a quaternion */
    
    int CQRInverse (CQRQuaternionHandle inversequaternion, CQRQuaternionHandle quaternion );
    
    /* CQRRotateByQuaternion -- Rotate a vector by a Quaternion, w = qvq* */
    
    int CQRRotateByQuaternion(double FAR * w, CQRQuaternionHandle rotquaternion, double FAR * v);        
    
    /* CQRAxis2Quaternion -- Form the quaternion for a rotation around axis v  by angle theta */
    
    int CQRAxis2Quaternion (CQRQuaternionHandle rotquaternion, double FAR * v, double theta);
    
    /* CQRMatrix2Quaterion -- Form the quaternion from a 3x3 rotation matrix R */
    
    int CQRMatrix2Quaternion (CQRQuaternionHandle rotquaternion, double R[3][3]);
    
    /* CQRQuaternion2Matrix -- Form the 3x3 rotation matrix from a quaternion */
    
    int CQRQuaternion2Matrix (double R[3][3], CQRQuaternionHandle rotquaternion);
    
    /* CQRQuaternion2Angles -- Convert a Quaternion into Euler Angles for Rz(Ry(Rx))) convention */
    
    int CQRQuaternion2Angles (double FAR * RotX, double FAR * RotY, double FAR * RotZ, CQRQuaternionHandle rotquaternion);
    
    /* CQRAngles2Quaternion -- Convert Euler Angles for Rz(Ry(Rx))) convention into a quaternion */
    
    int CQRAngles2Quaternion (CQRQuaternionHandle rotquaternion, double RotX, double RotY, double RotZ );
    
#ifdef __cplusplus
    
}

#endif


#endif
