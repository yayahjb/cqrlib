                 CQRlib -- ANSI C API for Quaternion Rotations

                                 Release 1.0.4
                                 25 April 2010
                   (c) 2008, 2009, 2010 Herbert J. Bernstein
                      yaya at bernstein-plus-sons dot com
               You may distribute the CQRlib API under the LGPL

   The 1.0.4 release added a version of L. Andrews adaptation to a C++
   template. The 1.0.3 release changed from use of a FAR macro to use of a
   CQR_FAR macro to avoid name conflicts. the macros for malloc, free,
   memmove and memset were also changed. The 1.0.2 release of 14 June 2009
   corrected the Makefile for case-sensitive file systems and to include -lm
   in loading. Release 1.0.1 of 23 February 2009 was a minor documentation
   update to the original 1.0 release of 22 February 2009.

   CQRlib is an ANSI C implementation of a utility library for quaternion
   arithmetic and quaternion rotation math. See
     * "Quaternions and spatial rotation", Wikipedia
       http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
     * K. Shoemake, "Quaternions", Department of Computer Science, University
       of Pennsylvania, Philadelphia, PA 19104,
       ftp://ftp.cis.upenn.edu/pub/graphics/shoemake/quatut.ps.Z
     * K. Shoemake, "Animating rotation with quaternion curves", ACM SIGGRAPH
       Computer Graphics, Vol 19, No. 3, pp 245--254, 1985.

   Work supported in part by NIH NIGMS under grant 1R15GM078077-01 and DOE 
   under grant ER63601-1021466-0009501. Any opinions, findings, and
   conclusions or recommendations expressed in this material are those of the 
   author(s) and do not necessarily reflect the views of the funding
   agencies.

    Installation

   The CQRlib package is available at www.sourceforge.net/projects/cqrlib. A
   source tarball is available at
   downloads.sourceforge.net/cqrlib/CQRlib-1.0.3.tar.gz. Later tarballs may
   be available.

   When the source tarball is downloaded and unpacked, you should have a
   directory CQRlib-1.0.4. To see the current settings for a build execute

   make

   which should give the following information:

 PLEASE READ README_CQRlib.txt and lgpl.txt
 
  Before making the CQRlib library and example programs, check
  that the chosen settings are correct
 
  The current C and C++ compile commands are:
 
    /Users/yaya/bin/libtool --mode=compile gcc -g -O2  -Wall -ansi -pedantic \
             -I.  -c
    /Users/yaya/bin/libtool --mode=compile g++ -g -O2  -Wall -ansi -pedantic \
             -DCQR_NOCCODE=1 -I.  -c
 
  The current library C and C++ link commands are:
 
    /Users/yaya/bin/libtool --mode=link gcc -version-info 2:0:0 -rpath \
             /Users/yaya/lib
    /Users/yaya/bin/libtool --mode=link g++ -version-info 2:0:0 -rpath \
             /Users/yaya/lib
 
  The current C library local, dynamic and static build commands are:
 
    /Users/yaya/bin/libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic -I.
    /Users/yaya/bin/libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic \
             -dynamic -I /Users/yaya/include -L/Users/yaya/lib
    /Users/yaya/bin/libtool --mode=link gcc -g -O2  -Wall -ansi -pedantic \
             -static -I /Users/yaya/include -L/Users/yaya/lib
 
  The current C++ template local, dynamic and static build commands are:
 
    /Users/yaya/bin/libtool --mode=link g++ -g -O2  -Wall -ansi -pedantic \
             -DCQR_NOCCODE=1 -I.
    /Users/yaya/bin/libtool --mode=link g++ -g -O2  -Wall -ansi -pedantic \
             -DCQR_NOCCODE=1 -dynamic -I /Users/yaya/include -L/Users/yaya/lib
    /Users/yaya/bin/libtool --mode=link g++ -g -O2  -Wall -ansi -pedantic \
             -DCQR_NOCCODE=1 -static -I /Users/yaya/include -L/Users/yaya/lib
 
  Before installing the CQRlib library and example programs, check
  that the install directory and install commands are correct:
 
  The current values are :
 
    /usr/local
    /usr/local/bin/libtool --mode=install cp
    
 
  To compile the CQRlib library and example programs type:
 
    make clean
    make all
 
  To run a set of tests type:
 
    make tests
 
  To clean up the directories type:
 
    make clean
 
  To install the library and binaries type:
 
    make install

   If these settings need to be changed, edit Makefile. On some systems, e.g.
   Mac OS X, the default libtool is not appropriate. In that case you should
   install a recent version of libtool. The CQRlib kit has been tested with
   libtool versions 1.3.5 and 1.5.4. If the system libtool is not to be used,
   define the variable LIBTOOL to be the path to the libtool executable, e.g.
   in bash

   export LIBTOOL=$HOME/bin/libtool

   of in the Makefie

   LIBTOOL = $(HOME)/bin/libtool

   If you need to include local header files using #include "..." instead of
   #include <...>, define the variable USE_LOCAL_HEADERS

    Synopsis

   #include <cqrlib.h>

     /* CQRCreateQuaternion -- create a quaternion = w +ix+jy+kz */
    
     int CQRCreateQuaternion(CQRQuaternionHandle * quaternion, double w, 
         double x, double y, double z);
    
     /* CQRCreateEmptyQuaternion -- create a quaternion = 0 +i0+j0+k0 */
    
     int CQRCreateEmptyQuaternion(CQRQuaternionHandle * quaternion) ;
    
     /* CQRFreeQuaternion -- free a quaternion  */
    
     int CQRFreeQuaternion(CQRQuaternionHandle * quaternion);       
    
     /* CQRSetQuaternion -- create an existing quaternion = w +ix+jy+kz */
    
     int CQRSetQuaternion( CQRQuaternionHandle quaternion, double w, double x,
          double y, double z);

     /*  CQRAdd -- add a quaternion (q1) to a quaternion (q2) */
    
     int CQRAdd (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, 
         CQRQuaternionHandle q2 );
    
     /*  CQRSubtract -- subtract a quaternion (q2) from a quaternion (q1)  */
    
     int CQRSubtract (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1,
          CQRQuaternionHandle q2 );
    
     /*  CQRMultiply -- multiply a quaternion (q1) by quaternion (q2)  */
    
     int CQRMultiply (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1,
          CQRQuaternionHandle q2 );
    
     /*  CQRDivide -- Divide a quaternion (q1) by quaternion (q2)  */
    
     int CQRDivide (CQRQuaternionHandle quaternion,  CQRQuaternionHandle q1, 
         CQRQuaternionHandle q2 );

     /*  CQRScalarMultiply -- multiply a quaternion (q) by scalar (s)  */
    
     int CQRScalarMultiply (CQRQuaternionHandle quaternion,  
         CQRQuaternionHandle q, double s );

     /*  CQREqual -- return 0 if quaternion q1 == q2  */
    
     int CQREqual (CQRQuaternionHandle q1, CQRQuaternionHandle q2 );
    
     /*  CQRConjugate -- Form the conjugate of a quaternion qconj */

     int CQRConjugate (CQRQuaternionHandle qconjugate, 
         CQRQuaternionHandle quaternion);
    
     /*  CQRNormsq -- Form the normsquared of a quaternion */
    
     int CQRNormsq (double * normsq, CQRQuaternionHandle quaternion ) ;
    
     /*  CQRNorm -- Form the norm of a quaternion */
    
     int CQRNorm (double * norm, CQRQuaternionHandle quaternion ) ;
    
     /*  CQRInverse -- Form the inverse of a quaternion */
    
     int CQRInverse (CQRQuaternionHandle inversequaternion, 
         CQRQuaternionHandle quaternion );
    
     /* CQRRotateByQuaternion -- Rotate a vector by a Quaternion, w = qvq* */
    
     int CQRRotateByQuaternion(double * w, 
         CQRQuaternionHandle rotquaternion, double * v);       
    
     /* CQRAxis2Quaternion -- Form the quaternion for a rotation around 
        axis v  by angle theta */
    
     int CQRAxis2Quaternion (CQRQuaternionHandle rotquaternion, 
        double * v, double theta);
    
     /* CQRMatrix2Quaterion -- Form the quaternion from a 3x3 rotation 
        matrix R */
    
     int CQRMatrix2Quaternion (CQRQuaternionHandle rotquaternion, 
         double R[3][3]);
    
     /* CQRQuaternion2Matrix -- Form the 3x3 rotation matrix from a 
        quaternion */
    
     int CQRQuaternion2Matrix (double R[3][3], 
         CQRQuaternionHandle rotquaternion);
    
     /* CQRQuaternion2Angles -- Convert a Quaternion into Euler Angles for 
        Rz(Ry(Rx))) convention */
    
     int CQRQuaternion2Angles (double * RotX, double * RotY, 
         double * RotZ, CQRQuaternionHandle rotquaternion);
    
     /* CQRAngles2Quaternion -- Convert Euler Angles for Rz(Ry(Rx))) 
        convention into a quaternion */
    
     int CQRAngles2Quaternion (CQRQuaternionHandle rotquaternion, 
         double RotX, double RotY, double RotZ );

   and for C++

 template< typename DistanceType=double, typename VectorType=double[3], 
   typename MatrixType=double[9] >
 class CPPQR
 {

 public:

      /* Constructors  */
          inline CPPQR( );  // default constructor
          inline CPPQR( const CPPQR& q ); // copy constructor
          inline CPPQR( const DistanceType& wi, const DistanceType& xi, 
            const DistanceType& yi, const DistanceType& zi );

      /* Set -- set the values of an existing quaternion = w +ix+jy+kz */
          inline void Set ( const DistanceType& wi, const DistanceType& xi, 
            const DistanceType& yi, const DistanceType& zi );

      /* Accessors */
          inline DistanceType GetW( void ) const;
          inline DistanceType GetX( void ) const;
          inline DistanceType GetY( void ) const;
          inline DistanceType GetZ( void ) const;
         
      /* Operators */
          inline CPPQR operator+ ( const CPPQR& q ) const;
          inline CPPQR& operator+= ( const CPPQR& q );
          inline CPPQR& operator-= ( const CPPQR& q );
          inline CPPQR operator- ( const CPPQR& q ) const;
          inline CPPQR operator* ( const CPPQR& q ) const;
          inline CPPQR operator/ ( const CPPQR& q2 ) const;
          inline CPPQR operator* ( const DistanceType& d ) const;
          inline CPPQR operator/ ( const DistanceType& d ) const;
          inline CPPQR Conjugate ( void ) const;
          inline CPPQR& operator= ( const CPPQR& q );
          inline bool operator== ( const CPPQR& q ) const;
          inline bool operator!= ( const CPPQR& q ) const;
          inline VectorType& operator* ( const VectorType& v );
          DistanceType operator[] ( const int k ) const;


      /* Normsq -- Form the normsquared of a quaternion */
          inline DistanceType Normsq ( void ) const;

      /* Norm -- Form the norm of a quaternion */
          inline DistanceType Norm ( void ) const;

      /* Inverse -- Form the inverse of a quaternion */
          inline CPPQR Inverse ( void ) const;

      /* RotateByQuaternion -- Rotate a vector by a Quaternion, w = qvq* */
          inline void RotateByQuaternion(VectorType &w, const VectorType v );
          inline VectorType& RotateByQuaternion( const VectorType v );

      /* Axis2Quaternion -- Form the quaternion for a rotation around 
         axis v  by angle theta */
          static inline CPPQR Axis2Quaternion ( const DistanceType& angle, 
            const VectorType v );
          static inline CPPQR Axis2Quaternion ( const VectorType v, 
            const DistanceType& angle  );

      /* Matrix2Quaterion -- Form the quaternion from a 3x3 rotation 
         matrix R */
          static inline void Matrix2Quaternion ( CPPQR& rotquaternion, 
            const MatrixType m );
          static inline void Matrix2Quaternion ( CPPQR& rotquaternion, 
            const DistanceType R[3][3] );

      /* Quaternion2Matrix -- Form the 3x3 rotation matrix from a 
         quaternion */   
          static inline void Quaternion2Matrix( MatrixType& m, const CPPQR q );
          static inline void Quaternion2Matrix( DistanceType m[3][3], 
            const CPPQR q );

      /* Get a unit quaternion from a general one */
          inline CPPQR UnitQ( void ) const;

      /* Quaternion2Angles -- Convert a Quaternion into Euler Angles for 
         Rz(Ry(Rx))) convention */ 
          inline bool Quaternion2Angles ( DistanceType& rotX, 
            DistanceType& rotY, DistanceType& rotZ ) const;

      /* Angles2Quaternion -- Convert Euler Angles for Rz(Ry(Rx))) 
         convention into a quaternion */
          static inline CPPQR Angles2Quaternion ( const DistanceType& rotX, 
            const DistanceType& rotY, const DistanceType& rotZ );
          static inline CPPQR Point2Quaternion( const DistanceType v[3] );

 }; // end class CPPQR



    Description

   The cqrlib.h header file defines the CQRQuaternionHandle type as a pointer
   to a struct of the CQRQuaternion type:

     typedef struct {
         double w;
         double x;
         double y;
         double z; } CQRQuaternion;

   representing w + xi +yj + zk. A quaternion may be declared directly using
   the CQRQuaternion type or dynamically allocated by CQRCreateQuaternion or
   CQRCreateEmptyQuaternion, in which case it is a user responsibility to
   eventually free the allocated memory with CQRFreeQuaternion. The
   components of an existing quaternion may be set by CQRSetQuaternion.

   The rules of quaternion arithmetic are applied:

   -1 = i*i = j*j = k*k, i*j=k=-j*i, j*k=i=-j*k, k*i=j=-i*k

   by CQRAdd, CQRSubtract, CQRMultiply and CQRDivide. CQRScalarMultiply
   multiplies a quaternion by a scalar.

   CQREqual returns 0 if quaternion q1 == q2, component by component.
   CQRConjugate computes a quaternion with the same scalar component and the
   negative of the vector component. CQRNormsq computes the sum of the
   squares of the components. CQRInverse computes the inverse of a non-zero
   quaternion.

   In handling rotations, a right-handed system is assumed.
   CQRRotateByQuaternion rotates a vector by a quaternion, w = qvq*.
   CQRAxis2Quaternion forms the quaternion for a rotation around axis v by
   angle theta. CQRMatrix2Quaterion forms the quaternion equivalent a 3x3
   rotation matrix R. CQRQuaternion2Matrix forms a 3x3 rotation matrix from a
   quaternion. CQRQuaternion2Angles converts a quaternion into Euler Angles
   for the Rz(Ry(Rx))) convention. CQRAngles2Quaternion convert Euler angles
   for the Rz(Ry(Rx))) convention into a quaternion.

   If operating with __cplusplus defined, then the CPPQR template is defined
   allowing the creation of CPPQR quaternion objects. The template has three
   typename arguments: DistanceType, VectorType and MatrixType that default
   to double, double[3] and double[9]. Specializations are provided to
   support a double[3][3] MatrixType.

    Returns

   The CQRlib functions return 0 for normal completion, or the sum of one or
   more of the following non-zero error codes:

     Error Return     Numeric Value    Meaning                                
     CQR_BAD_ARGUMENT    1             /* An argument is not valid */         
     CQR_NO_MEMORY       2             /* A call to allocate memory failed */ 
     CQR_FAILED          4             /* Operation failed */                 

    Examples

   To create a quaternion dynamically from memory, initialized as the x
   vector with a zero scalar value, reporting failure to stderr:

         #include <cqrlib.h>
         #include <stdio.h>
         ...
         CQRQuaternionHandle quathandle;
         ...
         if (CQRCreateQuaternion(&quathandle,0.,1.,0.,0.)) fprintf(stderr," CQRCreateQuaternion failed!!\n");

   To create an x vector quaternion, a y vector quaternion, add then together
   and multiply by a z-vector, and print the result :

         #include <cqrlib.h>
         #include <stdio.h>
         ...
         CQRQuaternion qx, qy, qz, qresult1, qresult2;
         ...
         if (CQRSetQuaternion(&qx,0.,1.,0.,0.)
           ||CQRSetQuaternion(&qy,0.,0.,1.,0.)
           ||CQRSetQuaternion(&qz,0.,0.,0.,1.)) fprintf(stderr," CQRSetQuaternion failed!!\n");
         if (CQRAdd(&qresult1,&qx,&qy)||CQRMultiply(&qresult2,&qresult1,&qz))
           fprintf(stderr," CQR Add or Multiply failed!!\n");
         fprintf(stdout,"Result = ((i+j)*k) = %g %+gi %+gj + %+gk\n",
           qresult2.w, qresult2.x, qresult2.y, qresult2.z);

   The output should be "Result = ((i+j)*k) = 0 +1i -1j +0k".

   To rotate the 3D vector [-1.,0.,1.] 90 degrees clockwise around the vector
   [1.,1.,1.]:

         #include <cqrlib.h>
         #include <math.h>
         #include <stdio.h>
         ...
         double axis[3] = {1.,1.,1.};
         double vector[3] = {-1.,0.,1.};
         double result[3];
         CQRQuaternion rotquat;
        
         double PI;
         PI = 4.*atan2(1.,1.);

         CQRAxis2Quaternion(&rotquat,axis,PI/2);
         CQRRotateByQuaternion(result, &rotquat, vector);
         ...
         fprintf(stdout," [-1.,0.,1.] rotated 90 degrees clockwise"
         " around the vector [1.,1.,1.] = [%g, %g, %g]\n",
         result[0], result[1], result[2]);

   The output should be "[-1.,0.,1.] rotated 90 degrees clockwise around the
   vector [1.,1.,1.] = [0.57735, -1.1547, 0.57735]".

   See the test program CQRlibTest.c.

   For examples of the use of the CPPQR template, see the C++ test program
   CPPQRTest.cpp.

     ----------------------------------------------------------------------

   Updated 25 April 2010
   yaya at bernstein-plus-sons dot com
