/*!
 * \file linear_solvers_structure.hpp
 * \brief Headers for the classes related to linear solvers (CG, FGMRES, etc)
 *        The subroutines and functions are in the <i>linear_solvers_structure.cpp</i> file.
 */

#ifndef ARIES_MATH_LINEARSOLVER_HPP
#define ARIES_MATH_LINEARSOLVER_HPP

#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include <climits>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <string>


#include "MATH_Vector.hpp"
#include "MATH_Matrix.hpp"
#include "../Common/TBOX_Config.hpp"
#include "../Geometry/GEOM_Geometry.hpp"


/*!
 * \class CSysSolve
 * \brief Class for solving linear systems using classical and Krylov-subspace iterative methods
 * \author J. Hicken.
 * \version 3.2.9 "eagle"
 *
 * The individual solvers could be stand-alone subroutines, but by
 * creating CSysSolve objects we can more easily assign different
 * matrix-vector products and preconditioners to different problems
 * that may arise in a hierarchical solver (i.e. multigrid).
 */

namespace ARIES
{
    namespace MATH
    {
        class MATH_LinearSolver
        {
        private:
            /*!
             * \brief sign transfer function
             * \param[in] x - value having sign prescribed
             * \param[in] y - value that defined the sign
             *
             * this may already be defined as a global function somewhere, if
             * so, feel free to delete this and replace it as needed with the
             * appropriate global function
             */
            double Sign(const double & x, const double & y) const;

            /*!
             * \brief applys a Givens rotation to a 2-vector
             * \param[in] s - sine of the Givens rotation angle
             * \param[in] c - cosine of the Givens rotation angle
             * \param[in, out] h1 - first element of 2x1 vector being transformed
             * \param[in, out] h2 - second element of 2x1 vector being transformed
             */
            void ApplyGivens(const double & s, const double & c, double & h1, double & h2);

            /*!
             * \brief generates the Givens rotation matrix for a given 2-vector
             * \param[in, out] dx - element of 2x1 vector being transformed
             * \param[in, out] dy - element of 2x1 vector being set to zero
             * \param[in, out] s - sine of the Givens rotation angle
             * \param[in, out] c - cosine of the Givens rotation angle
             *
             * Based on givens() of SPARSKIT, which is based on p.202 of
             * "Matrix Computations" by Golub and van Loan.
             */
            void GenerateGivens(double & dx, double & dy, double & s, double & c);

            /*!
             * \brief finds the solution of the upper triangular system Hsbg*x = rhs
             *
             * \param[in] n - size of the reduced system
             * \param[in] Hsbg - upper triangular matrix
             * \param[in] rhs - right-hand side of the reduced system
             * \param[out] x - solution of the reduced system
             *
             * \pre the upper Hessenberg matrix has been transformed into a
             * triangular matrix.
             */
            void SolveReduced(const int & n, const std::vector<std::vector<double> > & Hsbg,
                const std::vector<double> & rhs, std::vector<double> & x);

            /*!
             * \brief Modified Gram-Schmidt orthogonalization
             * \author Based on Kesheng John Wu's mgsro subroutine in Saad's SPARSKIT
             *
             * \tparam Vec - a generic vector class
             * \param[in] i - index indicating which vector in w is being orthogonalized
             * \param[in, out] Hsbg - the upper Hessenberg begin updated
             * \param[in, out] w - the (i+1)th vector of w is orthogonalized against the
             *                    previous vectors in w
             *
             * \pre the vectors w[0:i] are orthonormal
             * \post the vectors w[0:i+1] are orthonormal
             *
             * Reothogonalization is performed if the cosine of the angle between
             * w[i+1] and w[k], k < i+1, is greater than 0.98.  The norm of the "new"
             * vector is kept in nrm0 and updated after operating with each vector
             *
             */
            void ModGramSchmidt(int i, std::vector<std::vector<double> > & Hsbg, std::vector<MATH_Vector> & w);

            /*!
             * \brief writes header information for a CSysSolve residual history
             * \param[in, out] os - ostream class object for output
             * \param[in] solver - string describing the solver
             * \param[in] restol - the target tolerance to solve to
             * \param[in] resinit - the initial residual norm (absolute)
             *
             * \pre the ostream object os should be open
             */
            void WriteHeader(const std::string & solver, const double & restol, const double & resinit);

            /*!
             * \brief writes residual convergence data for one iteration to a stream
             * \param[in] iter - current iteration
             * \param[in] res - the (absolute) residual norm value
             * \param[in] resinit - the initial residual norm
             *
             * \pre the ostream object os should be open
             */
            void WriteHistory(const int & iter, const double & res, const double & resinit);

        public:

            /*! \brief Conjugate Gradient method
             * \param[in] b - the right hand size vector
             * \param[in, out] x - on entry the intial guess, on exit the solution
             * \param[in] mat_vec - object that defines matrix-vector product
             * \param[in] precond - object that defines preconditioner
             * \param[in] tol - tolerance with which to solve the system
             * \param[in] m - maximum size of the search subspace
             * \param[in] monitoring - turn on priting residuals from solver to screen.
             */
            unsigned long CG_LinSolver(const MATH_Vector & b, MATH_Vector & x, MATH_MatrixVectorProduct & mat_vec,
                MATH_Preconditioner & precond, double tol,
                unsigned long m, bool monitoring);

            /*!
             * \brief Flexible Generalized Minimal Residual method
             * \param[in] b - the right hand size vector
             * \param[in, out] x - on entry the intial guess, on exit the solution
             * \param[in] mat_vec - object that defines matrix-vector product
             * \param[in] precond - object that defines preconditioner
             * \param[in] tol - tolerance with which to solve the system
             * \param[in] m - maximum size of the search subspace
             * \param[in] monitoring - turn on priting residuals from solver to screen.
             */
            unsigned long FGMRES_LinSolver(const MATH_Vector & b, MATH_Vector & x, MATH_MatrixVectorProduct & mat_vec,
                MATH_Preconditioner & precond, double tol,
                unsigned long m, double *residual, bool monitoring);

            /*!
           * \brief Biconjugate Gradient Stabilized Method (BCGSTAB)
           * \param[in] b - the right hand size vector
           * \param[in, out] x - on entry the intial guess, on exit the solution
           * \param[in] mat_vec - object that defines matrix-vector product
           * \param[in] precond - object that defines preconditioner
           * \param[in] tol - tolerance with which to solve the system
           * \param[in] m - maximum size of the search subspace
           * \param[in] monitoring - turn on priting residuals from solver to screen.
           */
            unsigned long BCGSTAB_LinSolver(const MATH_Vector & b, MATH_Vector & x, MATH_MatrixVectorProduct & mat_vec,
                MATH_Preconditioner & precond, double tol,
                unsigned long m, double *residual, bool monitoring);

            /*!
             * \brief Solve the linear system using a Krylov subspace method
             * \param[in] Jacobian - Jacobian Matrix for the linear system
             * \param[in] LinSysRes - Linear system residual
             * \param[in] LinSysSol - Linear system solution
             * \param[in] geometry -  Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             */
            unsigned long Solve(MATH_Matrix & Jacobian, MATH_Vector & LinSysRes, MATH_Vector & LinSysSol, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

        };
    }
}

#include "MATH_LinearSolver.inl"

#endif