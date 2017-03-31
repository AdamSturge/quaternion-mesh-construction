// =============================================================================
// SpinXForm -- Quaternion.h
// Keenan Crane
// August 16, 2011
//
// Quaternion represents an element of the quaternions.  Operator overloading
// makes it possible to write expressions involving quaternions, vectors, and
// scalars.  For instance,
//
//    Quaternion q;
//    q = 1.23;
//
// sets the real part of q to 1.23 and the imaginary part to zero.  Similarly,
//
//    Quaternion q;
//    Eigen::Vector3d v;
//    q = v;
//
// sets the imaginary part to v and the real part to zero.
//

#ifndef SPINXFORM_QUATERNION_H
#define SPINXFORM_QUATERNION_H

#include <Eigen/Dense>
#include <ostream>

class Quaternion
{
   public:
      // CONSTRUCTORS ----------------------------------------------------------
      Quaternion( void );                                      // initializes all components to zero
      Quaternion( const Quaternion& q );                       // initializes from existing quaternion
      Quaternion( double s, double vi, double vj, double vk ); // initializes with specified real (s) and imaginary (v) components
      Quaternion( double s, const Eigen::Vector3d& v );                 // initializes with specified real (s) and imaginary (v) components
      Quaternion( double s );                                  // initializes purely real quaternion with specified real (s) component
      Quaternion( const Eigen::Vector3d& v );                           // initializes purely imaginary quaternion with specified imaginary (v) component
      
      // ASSIGNMENT OPERATORS --------------------------------------------------
      Quaternion& operator=( double s ); // assigns a purely real quaternion with real value s
      Quaternion& operator=( const Eigen::Vector3d& v ); // assigns a purely real quaternion with imaginary value v
			Quaternion& Quaternion :: operator=(const Quaternion& _q); // assigns values of one quaterion to this quaterion

      // ACCESSORS -------------------------------------------------------------
            double& operator[]( int index );       // returns reference to the specified component (0-based indexing: r, i, j, k)
      const double& operator[]( int index ) const; // returns const reference to the specified component (0-based indexing: r, i, j, k)
      void toMatrix( double Q[4][4] ) const;       // builds 4x4 matrix Q representing (left) quaternion multiplication
            double& re( void );                    // returns reference to double part
      const double& re( void ) const;              // returns const reference to double part
            Eigen::Vector3d& im( void );                    // returns reference to imaginary part
      const Eigen::Vector3d& im( void ) const;              // returns const reference to imaginary part

      // VECTOR SPACE OPERATIONS -----------------------------------------------
      Quaternion operator+( const Quaternion& q ) const; // addition
      Quaternion operator-( const Quaternion& q ) const; // subtraction
      Quaternion operator-( void ) const;                // negation
      Quaternion operator*( double c ) const;            // scalar multiplication
      Quaternion operator/( double c ) const;            // scalar division
      void       operator+=( const Quaternion& q );      // addition / assignment
      void       operator+=( double c );                 // addition / assignment of pure real
      void       operator-=( const Quaternion& q );      // subtraction / assignment
      void       operator-=( double c );                 // subtraction / assignment of pure real
      void       operator*=( double c );                 // scalar multiplication / assignment
      void       operator/=( double c );                 // scalar division / assignment

      // ALGEBRAIC OPERATIONS --------------------------------------------------
      Quaternion operator*( const Quaternion& q ) const; // Hamilton product
      void       operator*=( const Quaternion& q );      // Hamilton product / assignment
      Quaternion operator~( void ) const;                // conjugation
      Quaternion inv( void ) const;                      // inverse
      
      // NORMS -----------------------------------------------------------------
      double norm( void ) const;     // returns Euclidean length
      double norm2( void ) const;    // returns Euclidean length squared
      Quaternion unit( void ) const; // returns unit quaternion
      void normalize( void );        // divides by Euclidean length

   protected:
      // STORAGE ---------------------------------------------------------------
      double s; // scalar (double) part
      Eigen::Vector3d v; // vector (imaginary) part
};

// VECTOR SPACE OPERATIONS -----------------------------------------------
Quaternion operator*( double c, const Quaternion& q ); // scalar multiplication

// I/O -------------------------------------------------------------------------
std::ostream& operator<<( std::ostream& os, const Quaternion& q ); // prints components

#endif

