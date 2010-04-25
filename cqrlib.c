/*
 *  cqrlib.c
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
 * YOU MAY REDISTRIBUTE THE CQRlib API UNDER THE TERMS OF THE LGPL   *
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


#ifdef __cplusplus

extern "C" {
    
#endif
    
#include <cqrlib.h>
    
    /* CQRCreateQuaternion -- create a quaternion = w +ix+jy+kz */
    
    int CQRCreateQuaternion(CQRQuaternionHandle FAR * quaternion, double w, double x, double y, double z) {
        
        *quaternion = (CQRQuaternionHandle)MALLOC(sizeof(CQRQuaternion));
        
        if (!*quaternion) return CQR_NO_MEMORY;
        
        (*quaternion)->w = w;
        (*quaternion)->x = x;
        (*quaternion)->y = y;
        (*quaternion)->z = z;
        
        return CQR_SUCCESS;
        
    }
    
    /* CQRCreateEmptyQuaternion -- create a quaternion = 0 +i0+j0+k0 */
    
    int CQRCreateEmptyQuaternion(CQRQuaternionHandle FAR * quaternion) {
        
        *quaternion = (CQRQuaternionHandle)MALLOC(sizeof(CQRQuaternion));
        
        if (!*quaternion) return CQR_NO_MEMORY;
        
        
        (*quaternion)->w = 0;
        (*quaternion)->x = 0;
        (*quaternion)->y = 0;
        (*quaternion)->z = 0;
        
        return CQR_SUCCESS;
        
    }
    
    /* CQRFreeQuaternion -- free a quaternion  */
    
    int CQRFreeQuaternion(CQRQuaternionHandle FAR * quaternion) {
        
        if (!quaternion) return CQR_BAD_ARGUMENT;
        
        FREE(*quaternion);
        
        *quaternion = NULL;
        
        return CQR_SUCCESS;
        
    }
    
    /* CQRSetQuaternion -- create an existing quaternion = w +ix+jy+kz */
    
    int CQRSetQuaternion( CQRQuaternionHandle quaternion, double w, double x, double y, double z) {
        
        if (!quaternion) return CQR_BAD_ARGUMENT;
        
        quaternion->w = w;
        quaternion->x = x;
        quaternion->y = y;
        quaternion->z = z;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQRAdd -- add a quaternion (q1) to a quaternion (q2) */
    
    int CQRAdd (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 ) {
        
        if (!quaternion || !q1 || !q2 ) return CQR_BAD_ARGUMENT;
        
        quaternion->w = q1->w + q2->w;
        quaternion->x = q1->x + q2->x;
        quaternion->y = q1->y + q2->y;
        quaternion->z = q1->z + q2->z;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQRSubtract -- subtract a quaternion (q2) from a quaternion (q1)  */
    
    int CQRSubtract (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 ) {
        
        if (!quaternion || !q1 || !q2 ) return CQR_BAD_ARGUMENT;
        
        quaternion->w = q1->w - q2->w;
        quaternion->x = q1->x - q2->x;
        quaternion->y = q1->y - q2->y;
        quaternion->z = q1->z - q2->z;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQRMultiply -- multiply a quaternion (q1) by quaternion (q2)  */
    
    int CQRMultiply (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 ) {
        
        if (!quaternion || !q1 || !q2 ) return CQR_BAD_ARGUMENT;
        
        quaternion->w = -q1->z*q2->z - q1->y*q2->y - q1->x*q2->x + q1->w*q2->w;
        quaternion->x =  q1->y*q2->z - q1->z*q2->y + q1->w*q2->x + q1->x*q2->w;
        quaternion->y = -q1->x*q2->z + q1->w*q2->y + q1->z*q2->x + q1->y*q2->w;
        quaternion->z =  q1->w*q2->z + q1->x*q2->y - q1->y*q2->x + q1->z*q2->w;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQRDivide -- divide a quaternion (q1) by quaternion (q2)  */
    
    int CQRDivide (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, CQRQuaternionHandle q2 ) {
        
        double norm2sq;
        CQRQuaternion q;
        
        if (!quaternion || !q1 || !q2 ) return CQR_BAD_ARGUMENT;
        
        norm2sq = q2->w*q2->w+q2->x*q2->x+q2->y*q2->y+q2->z*q2->z;
        
        if (norm2sq==0.) return CQR_BAD_ARGUMENT;
        
        q.w =  q1->z*q2->z + q1->y*q2->y + q1->x*q2->x + q1->w*q2->w;
        q.x = -q1->y*q2->z + q1->z*q2->y - q1->w*q2->x + q1->x*q2->w;
        q.y =  q1->x*q2->z - q1->w*q2->y - q1->z*q2->x + q1->y*q2->w;
        q.z = -q1->w*q2->z - q1->x*q2->y + q1->y*q2->x + q1->z*q2->w;
        
        quaternion->w=q.w/norm2sq;
        quaternion->x=q.x/norm2sq;
        quaternion->y=q.y/norm2sq;
        quaternion->z=q.z/norm2sq;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQRScalarMultiply -- multiply a quaternion (q) by scalar (s)  */
    
    int CQRScalarMultiply (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q, double s ) {
        
        if (!quaternion || !q ) return CQR_BAD_ARGUMENT;
        
        quaternion->w = q->w *s;
        quaternion->x = q->x *s;
        quaternion->y = q->y *s;
        quaternion->z = q->z *s;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQREqual -- return 0 if quaternion q1 == q2  */
    
    int CQREqual (CQRQuaternionHandle q1, CQRQuaternionHandle q2 ) {

        if ( !q1 || !q2 ) return CQR_BAD_ARGUMENT;
        
        return ((q1->w==q2->w)&&(q1->x==q2->x)&&(q1->y==q2->y)&&(q1->z==q2->z))?CQR_SUCCESS:CQR_FAILED;

    }
    
    
    /*  CQRConjugate -- Form the conjugate of a quaternion qconj */
    
    int CQRConjugate (CQRQuaternionHandle qconjugate, CQRQuaternionHandle quaternion) {
        
        if (!quaternion || !qconjugate ) return CQR_BAD_ARGUMENT;
        
        qconjugate->w = quaternion->w;
        qconjugate->x = -quaternion->x;
        qconjugate->y = -quaternion->y;
        qconjugate->z = -quaternion->z;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQRNormsq -- Form the normsquared of a quaternion */
    
    int CQRNormsq (double * normsq, CQRQuaternionHandle quaternion ) {
        
        if (!quaternion || !normsq ) return CQR_BAD_ARGUMENT;
        
        *normsq = quaternion->w*quaternion->w + quaternion->x*quaternion->x + quaternion->y*quaternion->y + quaternion->z*quaternion->z;
        
        return CQR_SUCCESS;
        
    }
    
    /*  CQRInverse -- Form the inverse of a quaternion */
    
    int CQRInverse (CQRQuaternionHandle inversequaternion, CQRQuaternionHandle quaternion ) {
        
        double normsq;
        
        if (!quaternion || !inversequaternion ) return CQR_BAD_ARGUMENT;
        
        CQRMConjugate(*inversequaternion,*quaternion);
        CQRMNormsq(normsq,*quaternion);
        if (normsq > 0.) {
            CQRMScalarMultiply(*inversequaternion,*inversequaternion,1./normsq);
            return CQR_SUCCESS;
        } else return CQR_BAD_ARGUMENT;
        
    }
    
    /* CQRRotateByQuaternion -- Rotate a vector by a Quaternion, w = qvq* */
    
    int CQRRotateByQuaternion(double FAR * w, CQRQuaternionHandle rotquaternion, double FAR * v) {
        
        CQRQuaternion vquat, wquat, qconj;
        
        if (!w || !rotquaternion || !v ) return CQR_BAD_ARGUMENT;
        
        CQRMSet(vquat,0,v[0],v[1],v[2]);
       
        CQRMMultiply(wquat,*rotquaternion,vquat);
        CQRMConjugate(qconj,*rotquaternion);
        CQRMMultiply(vquat,wquat,qconj);
        
        w[0] = vquat.x;
        w[1] = vquat.y;
        w[2] = vquat.z;
        
        return CQR_SUCCESS;
        
    }
    
    /* CQRAxis2Quaternion -- Form the quaternion for a rotation around axis v  by angle theta */
    
    int CQRAxis2Quaternion (CQRQuaternionHandle rotquaternion, double FAR * v, double theta) {
        
        double normsq, norm;
        
        if (!rotquaternion || !v ) return CQR_BAD_ARGUMENT;
        
        normsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        
        if (normsq == 0.) return CQR_BAD_ARGUMENT;
        
        if (normsq == 1.) {
            
            CQRMSet(*rotquaternion,cos(theta/2),sin(theta/2)*v[0],sin(theta/2)*v[1],sin(theta/2)*v[2]);
            
        } else {
            
            norm = sqrt(normsq);
            
            CQRMSet(*rotquaternion,cos(theta/2),sin(theta/2)*v[0]/norm,sin(theta/2)*v[1]/norm,sin(theta/2)*v[2]/norm);
        }
        
        return CQR_SUCCESS;
    }
    
    /* CQRMatrix2Quaternion -- Form the quaternion from a 3x3 rotation matrix R */
    
    int CQRMatrix2Quaternion (CQRQuaternionHandle rotquaternion, double R[3][3]) {
        
        double trace, recip, fourxsq, fourysq, fourzsq;
        
        trace = R[0][0] + R[1][1] + R[2][2];
        
        if (trace > -.75) {
            
            rotquaternion->w = sqrt((1.+trace))/2.;
            
            recip = .25/rotquaternion->w;
            
            rotquaternion->z = (R[1][0]-R[0][1])*recip;
            
            rotquaternion->y = (R[0][2]-R[2][0])*recip;
            
            rotquaternion->x = (R[2][1]-R[1][2])*recip;
            
            return CQR_SUCCESS;
            
        }
        
        fourxsq = 1.+ R[0][0] - R[1][1] - R[2][2];
        
        if (fourxsq >= .25) {
            
            rotquaternion->x = sqrt(fourxsq)/2.;
            
            recip = .25/rotquaternion->x;
            
            rotquaternion->y = (R[1][0]+R[0][1])*recip;
            
            rotquaternion->z = (R[0][2]+R[2][0])*recip;
            
            rotquaternion->w = (R[2][1]-R[1][2])*recip;
            
            return CQR_SUCCESS;
            
        }
        
        fourysq = 1.+ R[1][1] - R[0][0] - R[2][2];
        
        if (fourysq >= .25) {
            
            rotquaternion->y = sqrt(fourysq)/2.;
            
            recip = .25/rotquaternion->y;
            
            rotquaternion->x = (R[1][0]+R[0][1])*recip;
            
            rotquaternion->w = (R[0][2]-R[2][0])*recip;
            
            rotquaternion->z = (R[2][1]+R[1][2])*recip;
            
            return CQR_SUCCESS;
            
        }
        
        fourzsq = 1.+ R[2][2] - R[0][0] - R[1][1];
        
        if (fourzsq >= .25) {
            
            rotquaternion->z = sqrt(fourzsq)/2.;
            
            recip = .25/rotquaternion->z;
            
            rotquaternion->w = (R[1][0]-R[0][1])*recip;
            
            rotquaternion->x = (R[0][2]+R[2][0])*recip;
            
            rotquaternion->y = (R[2][1]+R[1][2])*recip;
            
            return CQR_SUCCESS;
            
        }
        
        return CQR_BAD_ARGUMENT;
        
    }
    
    /* CQRQuaternion2Matrix -- Form the 3x3 rotation matrix from a quaternion */
    
    int CQRQuaternion2Matrix (double R[3][3], CQRQuaternionHandle rotquaternion) {
        
        double twoxy, twoyz, twoxz, twowx, twowy, twowz;
        double ww, xx, yy, zz;
        
        if (!R || !rotquaternion) return CQR_BAD_ARGUMENT;
        
        ww = (rotquaternion->w)*(rotquaternion->w);
        xx = (rotquaternion->x)*(rotquaternion->x);
        yy = (rotquaternion->y)*(rotquaternion->y);
        zz = (rotquaternion->z)*(rotquaternion->z);
        
        R[0][0] = ww + xx - yy - zz;
        R[1][1] = ww - xx + yy - zz;
        R[2][2] = ww - xx - yy + zz;
        
        twoxy = 2.*(rotquaternion->x)*(rotquaternion->y);
        twoyz = 2.*(rotquaternion->y)*(rotquaternion->z);
        twoxz = 2.*(rotquaternion->x)*(rotquaternion->z);
        twowx = 2.*(rotquaternion->w)*(rotquaternion->x);
        twowy = 2.*(rotquaternion->w)*(rotquaternion->y);
        twowz = 2.*(rotquaternion->w)*(rotquaternion->z);
        
        R[0][1] = twoxy - twowz;
        R[0][2] = twoxz + twowy;
        R[1][0] = twoxy + twowz;
        R[1][2] = twoyz - twowx;
        R[2][0] = twoxz - twowy;  
        R[2][1] = twoyz + twowx;
        
        return CQR_SUCCESS;
    }
    
    /* CQRQuaternion2Angles -- Convert a Quaternion into Euler Angles for Rz(Ry(Rx))) convention */
    
    int CQRQuaternion2Angles (double FAR * RotX, double FAR * RotY, double FAR * RotZ, CQRQuaternionHandle rotquaternion) {
        
        double SRX, SRY, SRZ, TRX, TRY, TRZ;
        double NSum;
        double TSum;
        
        
        double RMX0, RMX1, RMY0, RMY1, RMZ0, RMZ1, RMZ2;
        
        double PI;
        
        PI = 4.*atan2(1.,1.);
        
        if (!rotquaternion || !RotX || !RotY || !RotZ) return CQR_BAD_ARGUMENT;
        
        RMX0 = rotquaternion->w*rotquaternion->w + rotquaternion->x*rotquaternion->x
        - rotquaternion->y*rotquaternion->y - rotquaternion->z*rotquaternion->z;
        RMX1 = 2.*(rotquaternion->x*rotquaternion->y-rotquaternion->w*rotquaternion->z);
        RMY0 = 2.*(rotquaternion->x*rotquaternion->y+rotquaternion->w*rotquaternion->z);
        RMY1 = rotquaternion->w*rotquaternion->w - rotquaternion->x*rotquaternion->x
        + rotquaternion->y*rotquaternion->y - rotquaternion->z*rotquaternion->z;
        RMZ0 = 2.*(rotquaternion->x*rotquaternion->z-rotquaternion->w*rotquaternion->y);
        RMZ1 = 2.*(rotquaternion->w*rotquaternion->x+rotquaternion->y*rotquaternion->z);
        RMZ2 = rotquaternion->w*rotquaternion->w - rotquaternion->x*rotquaternion->x
        - rotquaternion->y*rotquaternion->y + rotquaternion->z*rotquaternion->z;
        
        if (RMZ0 < 1. ) {
            if (RMZ0 > -1.) {
                SRY = asin(-RMZ0);
            } else {
                SRY = .5*PI;
            }
        } else {
            SRY = -.5*PI;
        } 
        if (RMZ0 > .9999995) {
            SRX = atan2(-RMX1,RMY1);
            SRZ = 0;
        } else {
            if (RMZ0 < -.9999995 ) {
                SRX = atan2(RMX1,RMY1);
                SRZ = 0;
            } else {
                SRX = atan2(RMZ1,RMZ2);
                SRZ = atan2(RMY0,RMX0);
            }
        }

        TRX = PI+SRX;
        if ( TRX > 2.*PI ) TRX -= 2.*PI;
        TRY = PI+SRY;
        if ( TRY > 2.*PI ) TRY -= 2.*PI;
        TRZ = PI+SRZ;
        if ( TRZ > 2.*PI ) TRZ -= 2.*PI;

        NSum = 0;
        TSum = 0;
        NSum += fabs(cos(SRX)-cos(*RotX)) + fabs(sin(SRX)-sin(*RotX))
        + fabs(cos(SRY)-cos(*RotY)) + fabs(sin(SRY)-sin(*RotY))
        + fabs(cos(SRZ)-cos(*RotZ)) + fabs(sin(SRZ)-sin(*RotZ));
        TSum += fabs(cos(TRX)-cos(*RotX)) + fabs(sin(TRX)-sin(*RotX))
        + fabs(cos(TRY)-cos(*RotY)) + fabs(sin(TRY)-sin(*RotY))
        + fabs(cos(TRZ)-cos(*RotZ)) + fabs(sin(TRZ)-sin(*RotZ));
        
        if (NSum < TSum) {
            *RotX = SRX; *RotY = SRY; *RotZ = SRZ;
        } else {
            *RotX = TRX; *RotY = TRY; *RotZ = TRZ;
        }
        
        return CQR_SUCCESS;
    }
    
    /* CQRAngles2Quaternion -- Convert Euler Angles for Rz(Ry(Rx))) convention into a quaternion */
    
    int CQRAngles2Quaternion (CQRQuaternionHandle rotquaternion, double RotX, double RotY, double RotZ ) {
        
        double cx, cy, cz, sx, sy, sz;
        
        if (!rotquaternion) return CQR_BAD_ARGUMENT;
        
        cx = cos(RotX/2);
        sx = sin(RotX/2);
        
        cy = cos(RotY/2);
        sy = sin(RotY/2);
        
        cz = cos(RotZ/2);
        sz = sin(RotZ/2);
        
        rotquaternion->w = cx*cy*cz + sx*sy*sz;
        rotquaternion->x = sx*cy*cz - cx*sy*sz;
        rotquaternion->y = cx*sy*cz + sx*cy*sz;
        rotquaternion->z = cx*cy*sz - sx*sy*cz;
        
        return CQR_SUCCESS;
        
    }
    
#ifdef __cplusplus
    
}

#endif
