#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define USE_LOCAL_HEADERS
#ifndef USE_LOCAL_HEADERS
#include <cqrlib.h>
#else
#include "cqrlib.h"
#endif

int errorcount = 0;

int main ( )
{
    CPPQR<double> q1, q4, qx, qy, qz;
    double normsq;
    double PI;
    double vx[3] = {1.,0.,0.};
    double vy[3] = {0.,1.,0.};
    double vz[3] = {0.,0.,1.};
    double EXX,EXY,EXZ;
    double EYX,EYY,EYZ;
    double EZX,EZY,EZZ;
    
    PI = 4.*atan2(1.,1.);
    errorcount=0;
    
    if (q1.GetW() !=0.  ||  q1.GetX( ) !=0. || q1.GetY( ) !=0. || q1.GetZ( ) !=0.) {
        errorcount++;
        fprintf(stdout," CPPQRCreateEmptyQuaternion for q1 non-zero [ %g, %g, %g, %g ]\n",q1.GetW(),q1.GetX(),q1.GetY( ),q1.GetZ());
    }
    
    CPPQR<double> q2(1.,2.,3.,4.);
    
    if (q2.GetW() !=1.  ||  q2.GetX() !=2. || q2.GetY() !=3. || q2.GetZ() !=4.) {
        errorcount++;
        fprintf(stdout," CPPQRCreateQuaternion for q2 wrong value [ %g, %g, %g, %g ] != [1.,2,3.,4.]\n",q2.GetW(),q2.GetX(),q2.GetY(),q2.GetZ());
    }
    
    CPPQR<double> q3;
    
    q3 = q1 + q2;
    
    if (q3.GetW() !=1.  ||  q3.GetX() !=2. || q3.GetY() !=3. || q3.GetZ() !=4.) {
        fprintf(stdout," CPPQRAdd(q3,q1,q2) q3 wrong value [ %g, %g, %g, %g ] != [1.,2.,3.,4.]\n",q3.GetW(),q3.GetX(),q3.GetY(),q3.GetZ());
    }
    
    q1.Set( -9999.,-9998.,-9997.,-9996.);
    
    if (q1.GetW() !=-9999.  ||  q1.GetX() !=-9998. || q1.GetY() !=-9997. || q1.GetZ() !=-9996.) {
        errorcount++;
        fprintf(stdout," CPPQRSetQuaternion q1 wrong value [ %g, %g, %g, %g ] != [-9999.,-9998.,-9997.,-9996.]\n",q1.GetW(),q1.GetX(),q1.GetY(),q1.GetZ());
    }
    
    q1 = q3 - q2;
    
    if (q1.GetW() !=0.  ||  q1.GetX() !=0. || q1.GetY() !=0. || q1.GetZ() !=0.) {
        errorcount++;
        fprintf(stdout," CPPQR Subtract(q1,q3,q2) for q1 non-zero [ %g, %g, %g, %g ]\n",q1.GetW(),q1.GetX(),q1.GetY(),q1.GetZ());
    }
    
    q4.Set(1.,-2.,-3.,-4.);
    
    if (q4.GetW() != 1.  ||  q4.GetX() !=-2. || q4.GetY() !=-3. || q4.GetZ() !=-4.) {
        errorcount++;
        fprintf(stdout," CPPQRSetQuaternion &q4 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q4.GetW(),q4.GetX(),q4.GetY(),q4.GetZ());
    }
    
    q1 = q4 * q2;
    if (fabs(q1.GetW()-30.)>300.*DBL_EPSILON || fabs(q1.GetX())>300.*DBL_EPSILON || fabs(q1.GetY())>300.*DBL_EPSILON ||fabs(q1.GetZ())>300.*DBL_EPSILON )  {
        errorcount++;
        fprintf(stdout,"  CPPQRMultiply(q1,&q4,q2)  q1 wrong value [ %g, %g, %g, %g ] != [30.,0.,0.,0.]\n",q1.GetW(),q1.GetX(),q1.GetY(),q1.GetZ());
    }
    
    q3 = q1 / q2;
    if (fabs(q3.GetW()-1.)>60.*DBL_EPSILON || fabs(q3.GetX()+2.)>60.*DBL_EPSILON || fabs(q3.GetY()+3.)>60.*DBL_EPSILON ||fabs(q3.GetZ()+4.)>60.*DBL_EPSILON )  {
        errorcount++;
        fprintf(stdout,"  CPPQRDivide(q3,q1,q2)  q3 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q3.GetW(),q3.GetX(),q3.GetY(),q3.GetZ());
    }
    
    q3 = q2 * 3.0;
    if (fabs(q3.GetW()-3.)>180.*DBL_EPSILON || fabs(q3.GetX()-6.)>180.*DBL_EPSILON || fabs(q3.GetY()-9.)>180.*DBL_EPSILON ||fabs(q3.GetZ()-12.)>180.*DBL_EPSILON )  {
        errorcount++;
        fprintf(stdout,"  CPPQRScalarMultiply(q3,q2,3.)  q3 wrong value [ %g, %g, %g, %g ] != [3.,6.,9.,12.]\n",q3.GetW(),q3.GetX(),q3.GetY(),q3.GetZ());
    }
    
    q3 = q2.Conjugate( );
    if (q3.GetW() != 1.  ||  q3.GetX() !=-2. || q3.GetY() !=-3. || q3.GetZ() !=-4.) {
        errorcount++;
        fprintf(stdout," CPPQRConjugate(q3,q2) q3 wrong value [ %g, %g, %g, %g ] != [1.,-2.,-3.,-4.]\n",q3.GetW(),q3.GetX(),q3.GetY(),q3.GetZ());
    }
    
    if ( q4 != q3) {
        errorcount++;
        fprintf(stdout," CPPQREqual(&q4,q2) failed\n");
    }
    
    normsq = q4.Normsq();
    if ( fabs(normsq-30.) > 300.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout," CPPQRNormsq(&normsq,&q4) failed normsq=%g != 30\n", normsq);
    }
    
    q3 = q4.Inverse( );
    
    if (fabs(q3.GetW() - 1./30.) > 2.*DBL_EPSILON  ||  fabs(q3.GetX() - 2./30.) > 2.*DBL_EPSILON || fabs(q3.GetY() - 3./30.) > 2.*DBL_EPSILON || fabs(q3.GetZ() - 4./30.) > 2.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout," CPPQRInverse(q3,&q4) q3 wrong value [ %g, %g, %g, %g ] != [1./30.,2./30.,3./30,4./30.]\n",q3.GetW(),q3.GetX(),q3.GetY(),q3.GetZ());
    }
    
    
    /* Create quaternions to rotate about the x,y and z-axes by 90 degrees */
    
    qx = CPPQR<double>::Axis2Quaternion( vx, PI/2.0 );
    qy = CPPQR<double>::Axis2Quaternion( vy, PI/2.0 );
    qz = CPPQR<double>::Axis2Quaternion( vz, PI/2.0 );
    
    if (qx.GetW()<0.||fabs(qx.GetW()*qx.GetW()-.5)>10.*DBL_EPSILON||fabs(qx.GetX()*qx.GetX()-.5)>10.*DBL_EPSILON||fabs(qx.GetY())>10.*DBL_EPSILON||fabs(qx.GetZ())>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout,"Axis2Quaternion qx wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qx.GetW(),qx.GetX(),qx.GetY(),qx.GetZ());
    }
    
    
    if (qy.GetW()<0.||fabs(qy.GetW()*qy.GetW()-.5)>10.*DBL_EPSILON||fabs(qy.GetY()*qy.GetY()-.5)>10.*DBL_EPSILON||fabs(qy.GetX())>10.*DBL_EPSILON||fabs(qy.GetZ())>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout,"Axis2Quaternion qy wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qy.GetW(),qy.GetX(),qy.GetY(),qy.GetZ());
    }
    
    
    if (qz.GetW()<0.||fabs(qz.GetW()*qz.GetW()-.5)>10.*DBL_EPSILON||fabs(qz.GetZ()*qz.GetZ()-.5)>10.*DBL_EPSILON||fabs(qz.GetX())>10.*DBL_EPSILON||fabs(qz.GetY())>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout,"CPPAxis2Quaternion qz wrong value [ %g, %g, %g, %g ] != [sqrt(1./2.),sqrt(1./2.),0,0]\n",qz.GetW(),qz.GetX(),qz.GetY(),qz.GetZ());
    }
    
    double Matx[9];
    double Maty[9];
    double Matz[9];
    
    CPPQR<double>::Quaternion2Matrix(Matx,qx);
    CPPQR<double>::Quaternion2Matrix(Maty,qy);
    CPPQR<double>::Quaternion2Matrix(Matz,qz);
    
    if (  fabs(Matx[0]-1.)>10.*DBL_EPSILON ||fabs(Matx[5]+1.)>10.*DBL_EPSILON  ||fabs(Matx[7]-1.)>10.*DBL_EPSILON
        ||fabs(Matx[1])>10.*DBL_EPSILON    ||fabs(Matx[2])>10.*DBL_EPSILON
        ||fabs(Matx[3])>10.*DBL_EPSILON    ||fabs(Matx[4])>10.*DBL_EPSILON
        ||fabs(Matx[6])>10.*DBL_EPSILON    ||fabs(Matx[8])>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout," Quaternion2Matrix Matx wrong value \n  [ %g, %g, %g ]\n  [ %g, %g, %g ]\n  [ %g, %g, %g ]\n"
                "!=  [1, 0, 0]\n    [0, 0, -1]\n    [0, 1, 0]\n",
                Matx[0]-1.0,Matx[1],Matx[2],
                Matx[3],Matx[4],Matx[5]+1.0,
                Matx[6],Matx[7]-1.0,Matx[8]);
    }
    
    if (fabs(Maty[2]-1.)>10.*DBL_EPSILON ||fabs(Maty[4]-1.)>10.*DBL_EPSILON||fabs(Maty[6]+1.)>10.*DBL_EPSILON
        ||fabs(Maty[0])>10.*DBL_EPSILON ||fabs(Maty[1])>10.*DBL_EPSILON
        ||fabs(Maty[3])>10.*DBL_EPSILON ||fabs(Maty[5])>10.*DBL_EPSILON
        ||fabs(Maty[7])>10.*DBL_EPSILON ||fabs(Maty[8])>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout," Quaternion2Matrix Maty wrong value \n  [ %g, %g, %g]\n  [ %g, %g, %g]\n  [ %g, %g, %g ]\n"
                "!=  [0, 0, 1]\n    [0, 1, 0]\n    [-1, 0, 0]\n",
                Maty[0],Maty[1],Maty[2]-1.0,
                Maty[3],Maty[4],Maty[5]-1.0,
                Maty[6]+1.0,Maty[7],Maty[8]);
    }
    
    if (fabs(Matz[1]+1.)>10.*DBL_EPSILON ||fabs(Matz[3]-1.)>10.*DBL_EPSILON||fabs(Matz[8]-1.)>10.*DBL_EPSILON
        ||fabs(Matz[0])>10.*DBL_EPSILON ||fabs(Matz[2])>10.*DBL_EPSILON
        ||fabs(Matz[4])>10.*DBL_EPSILON ||fabs(Matz[5])>10.*DBL_EPSILON
        ||fabs(Matz[6])>10.*DBL_EPSILON ||fabs(Matz[7])>10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout," Quaternion2Matrix Matz wrong value \n  [ %g, %g, %g]\n  [ %g, %g, %g]\n  [ %g, %g, %g ]\n"
                "!=  [0, -1, 0]\n    [1, 0, 0]\n    [0, 0, 1]\n",
                Matz[0],Matz[1],Matz[2],
                Matz[3],Matz[4],Matz[5],
                Matz[6],Matz[7],Matz[8]);
    }
    
    
    EXX = EXY = EXZ = 0.;
    EYX = EYY = EYZ = 0.;
    EZX = EZY = EZZ = 0.;
    qx.Quaternion2Angles(EXX,EXY,EXZ);
    qy.Quaternion2Angles(EYX,EYY,EYZ);
    qz.Quaternion2Angles(EZX,EZY,EZZ);
    
    q1 = CPPQR<double>::Angles2Quaternion(EXX,EXY,EXZ);
    q2 = CPPQR<double>::Angles2Quaternion(EYX,EYY,EYZ);
    q3 = CPPQR<double>::Angles2Quaternion(EZX,EZY,EZZ);
    
    q4 = qx / q1;
    normsq = q4.Normsq( );
    
    if ( fabs(normsq-1.) > 10.*DBL_EPSILON || fabs(q4.GetW()*q4.GetW()-1.) > 10.*DBL_EPSILON) {
        errorcount++;
        fprintf(stdout," Angles2Quaternion q1 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qx.GetW(), qx.GetX(), qx.GetY(), qx.GetZ(), q1.GetW(), q1.GetX(), q1.GetY(), q1.GetZ() );  
    }
    
    q4 = qy / q2;
    normsq = q4.Normsq( );
    if (fabs(normsq-1.) > 10.*DBL_EPSILON || fabs(q4.GetW()*q4.GetW()-1.) > 10.*DBL_EPSILON) {
        fprintf(stdout," Angles2Quaternion q2 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qy.GetW(), qy.GetX(), qy.GetY(), qy.GetZ(), q2.GetW(), q2.GetX(), q2.GetY(), q2.GetZ() );  
    }
    
    q4 = qz / q3;
    normsq = q4.Normsq( );
    if (abs(normsq-1.) > 10.*DBL_EPSILON || fabs(q4.GetW()*q4.GetW()-1.) > 10.*DBL_EPSILON) {
        fprintf(stdout," Angles2Quaternion q3 wrong value [%g, %g, %g, %g] != +/-[%g, %g, %g, %g]\n",
                qz.GetW(), qz.GetX(), qz.GetY(), qz.GetZ(), q3.GetW(), q3.GetX(), q3.GetY(), q3.GetZ() );  
    }
    
    { /* lca */
        double m[9] = { 0.0,0.0,1.0, 1.0,0.0,0.0, 0.0,1.0,0.0 };
        double mx[9];
        double sum;
        
        CPPQR<double> qq1, qq2, qq3;
        CPPQR<double> qq4, qqx, qqy, qqz;
        
        CPPQR<double>::Matrix2Quaternion(qq4, m );
        CPPQR<double>::Quaternion2Matrix(mx,qq4 );
        
        sum = 0.0;
        for ( int i=0; i<9; ++i )
            sum += fabs(m[i]-mx[i]);
        if ( sum > 1.0E-8 )
        {
            errorcount++;
            fprintf( stdout, " Matrix2Quaternion difference\n" );
            fprintf( stdout, " qq4 = {%g,%g,%g,%g}\n",qq4.GetW(),qq4.GetX(),qq4.GetY(),qq4.GetZ());
            fprintf( stdout, " m = {{%g,%g,%g},{%g,%g,%g},{%g,%g,%g))\n",
                    m[0],m[1],m[2],m[3],m[4],m[5],m[6],m[7],m[8]);
            fprintf( stdout, " mx = {{%g,%g,%g},{%g,%g,%g},{%g,%g,%g))\n",
                    mx[0],mx[1],mx[2],mx[3],mx[4],mx[5],mx[6],mx[7],mx[8]);
        }
        
    }
    
    
    {
        CPPQR<double> q1( 3,5,7,11 );
        CPPQR<double> q2( q1.GetW(), q1.GetX(), q1.GetY(), q1.GetZ() );
        if ( q1 != q2 || !(q1==q2) )
        {
            errorcount++;
            fprintf( stdout, " CPPQR test constructors, gets, ==, != failed\n" );
        }
    }
    
    {
        CPPQR<double> q;
        if( ::fabs( q.GetW()) + ::fabs( q.GetX()) + ::fabs( q.GetY()) + ::fabs( q.GetZ()) != 0.0 )
        {
            errorcount++;
            fprintf( stdout, " CPPQR default constructor not zero\n" );
        }
    }
    
    {
        CPPQR<double> q1( 3,5,7,11 );
        CPPQR<double> q2( q1 );
        CPPQR<double> q3 = q2;
        if ( q1 != q2 || q1 != q3 )
        {
            errorcount++;
            fprintf( stdout, " CPPQR copy constructor or assignment operator failed\n" );
        }
    }
    
    {
        CPPQR<double> q1( 3,5,7,11 );
        CPPQR<double> q2;
        q2.Set( q1.GetW(), q1.GetX(), q1.GetY(), q1.GetZ() );
        if( q1 != q2 )
        {
            errorcount++;
            fprintf( stdout, " CPPQR Set failed\n" );
        }
    }
    
    {
        if( CPPQR<double>(1,3,5,7)+CPPQR<double>(-1,-3,-5,-7) !=CPPQR<double>(0,0,0,0) )
        {
            errorcount++;
            fprintf( stdout, " CPPQR add failed\n" );
        }
        if( CPPQR<double>(1,3,5,7)-CPPQR<double>(1,3,5,7) !=CPPQR<double>(0,0,0,0) )
        {
            errorcount++;
            fprintf( stdout, " CPPQR subtract failed\n" );
        }
    }
    
    {
        if ( CPPQR<double>(1,3,5,7)*2.0 != CPPQR<double>(4,12,20,28)/2.0 )
        {
            errorcount++;
            fprintf( stdout, " CPPQR multiply or divide by scalar failed\n" );
        }
    }
    
    {
        const CPPQR<double> q1( CPPQR<double>( 3,5,7,9 ) );
        const double normsq = q1.Normsq( );
        if ( normsq != 164.0 )
        {
            errorcount++;
            fprintf ( stdout, " CPPQR Normsq failed \n" );
        }
        
        const CPPQR<double> q2 = q1.UnitQ( );
        if( q1/sqrt(normsq) != q2 )
        {
            errorcount++;
            fprintf( stdout, "UnitQ failed\n" );
        }
    }
    
    {
        const CPPQR<double> q1( CPPQR<double>( 3,5,7,9 ) );
        if ( q1.GetW() != q1[0] || q1.GetX() != q1[1] || q1.GetY() != q1[2] || q1.GetZ() != q1[3] )
        {
            errorcount++;
            fprintf( stdout, "component fetches failed\n" );
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
        const CPPQR<double> q1( CPPQR<double>( -sqrt(7.),2,3,4 ));
        if ( q1.GetIm() != CPPQR<double>( 0,2,3,4 ) )
        {
            errorcount++;
            fprintf( stdout, "GetIm failed\n" );
        }
        
        
        if ( q1.GetAxis().Dist(CPPQR<double>( 0,2./sqrt(4.+9.+16.),3./sqrt(4.+9.+16.),4./sqrt(4.+9.+16.) ))>100.*DBL_EPSILON )
        {
            errorcount++;
            fprintf( stdout, "GetAxis failed\n" );
        }
        
        if (fabs(q1.GetAngle()-2.027462712001523)>40.*DBL_EPSILON*2.027462712001523)
        {
            errorcount++;
            fprintf( stdout, "GetAngle failed, got %g, expected %g\n",q1.GetAngle(),2.027462712001523 );
        }
        
        
        if ((q1.log() - CPPQR<double>(log(6.), 0.752980747892971, 1.129471121839457, 1.505961495785942)).Norm() > 
            40.*DBL_EPSILON*q1.log().Norm())
        {
            errorcount++;
            fprintf( stdout, "quaternion log failed log([%g,%g,%g,%g]) = [%g,%g,%g,%g] instead of [%g,%g,%g,%g], normdiff = %g\n",
                    q1.GetW(), q1.GetX(), q1.GetY(), q1.GetZ(),
                    q1.log().GetW(), q1.log().GetX(), q1.log().GetY(), q1.log().GetZ(),
                    log(6.), 0.752980747892971, 1.129471121839457, 1.505961495785942,
                    (q1.log() - CPPQR<double>(log(6.), 0.752980747892971, 1.129471121839457, 1.505961495785942)).Norm()
            );

        }
        
        if (((q1.log()).exp()-q1).Norm()>10.*DBL_EPSILON*q1.Norm() || 
            ((q1.exp()).log().exp()-q1.exp()).Norm()>10.*DBL_EPSILON*q1.exp().Norm())
        {
            errorcount++;
            fprintf( stdout, "log(exp) or exp(log) failed\n," 
                   " q = [%g,%g,%g,%g],"
                   " log = [%g,%g,%g,%g], exp(log) = [%g,%g,%g,%g],"
                   " exp = [%g,%g,%g,%g], log(exp) = [%g,%g,%g,%g]\n",
                    q1.GetW(),q1.GetX(),q1.GetY(),q1.GetZ(),
                    q1.log().GetW(),q1.log().GetX(),q1.log().GetY(),q1.log().GetZ(),
                    q1.log().exp().GetW(),q1.log().exp().GetX(),q1.log().exp().GetY(),q1.log().exp().GetZ(),
                    q1.exp().GetW(),q1.exp().GetX(),q1.exp().GetY(),q1.exp().GetZ(),
                    q1.exp().log().exp().GetW(),q1.exp().log().exp().GetX(),q1.exp().log().exp().GetY(),q1.exp().log().exp().GetZ()
            );
        }
        
        for (int i = -5; i < 6; i++) {
            if ((q1.pow(i) - q1.pow(double(i))).Norm() > 10.*DBL_EPSILON*(q1.pow(i)).Norm())
            {
                errorcount++;
                fprintf( stdout, "integer power double power comparison failed\n,"); 

            }
        }
        
        
    }
    
    {
        const CPPQR<double> q1( CPPQR<double>( -4.,0.,0.,0. ));
        const CPPQR<double> q2( CPPQR<double>( -4.,1.,1.,1. ));
        const CPPQR<double> q3( CPPQR<double>( 4.,0.,0.,0. ));
        CPPQR <double> qout1, qout2, qout3, qtest1, qtest2, qtest3;
        
        for (int i = 1; i < 9; i++) {
            for (int j = 0; j <  i; j++ ) {
                qout1 = q1.root(i,j);
                qout2 = q2.root(i,j);
                qout3 = q3.root(i,j);
                qtest1 = qout1.pow(i);
                qtest2 = qout2.pow(i);
                qtest3 = qout3.pow(i);
                if (q1.Dist(qtest1) > 100.*DBL_EPSILON*q1.Norm()
                    || q2.Dist(qtest2) > 100.*DBL_EPSILON*q2.Norm()
                    || q3.Dist(qtest3) > 100.*DBL_EPSILON*q3.Norm()) {
                    errorcount++;
                    fprintf(stdout," %d'th root of [%g,%g,%g,%g] = [%g,%g,%g,%g], power = [%g,%g,%g,%g], delta %g\n",
                        i, q1.GetW(), q1.GetX(), q1.GetY(), q1.GetZ(),
                        qout1.GetW(), qout1.GetX(), qout1.GetY(), qout1.GetZ(),
                        qtest1.GetW(), qtest1.GetX(), qtest1.GetY(), qtest1.GetZ(),
                        (q1-qtest1).Norm());
                    fprintf(stdout," %d'th root of [%g,%g,%g,%g] = [%g,%g,%g,%g], power = [%g,%g,%g,%g], delta %g\n",
                        i, q2.GetW(), q2.GetX(), q2.GetY(), q2.GetZ(),
                            qout2.GetW(), qout2.GetX(), qout2.GetY(), qout2.GetZ(),
                            qtest2.GetW(), qtest2.GetX(), qtest2.GetY(), qtest2.GetZ(),
                            (q2-qtest2).Norm());
                    fprintf(stdout," %d'th root of [%g,%g,%g,%g] = [%g,%g,%g,%g], power = [%g,%g,%g,%g], delta %g\n",
                            i, q3.GetW(), q3.GetX(), q3.GetY(), q3.GetZ(),
                            qout3.GetW(), qout3.GetX(), qout3.GetY(), qout3.GetZ(),
                            qtest3.GetW(), qtest3.GetX(), qtest3.GetY(), qtest3.GetZ(),
                            (q3-qtest3).Norm());
                    
                }
            }
        }

    }
    
    
    return errorcount;
}
