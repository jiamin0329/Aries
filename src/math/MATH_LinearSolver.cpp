/*!
 * \file linear_solvers_structure.cpp
 * \brief Main classes required for solving linear systems of equations
 */

#include "MATH_LinearSolver.hpp"


namespace ARIES
{
    namespace MATH
    {
        void MATH_LinearSolver::ApplyGivens(const double & s, const double & c, double & h1, double & h2) 
        {
            double temp = c*h1 + s*h2;
            h2 = c*h2 - s*h1;
            h1 = temp;
        }

        void MATH_LinearSolver::GenerateGivens(double & dx, double & dy, double & s, double & c) 
        {
            if ((dx == 0.0) && (dy == 0.0)) 
            {
                c = 1.0;
                s = 0.0;
            }
            else if (fabs(dy) > fabs(dx)) 
            {
                double tmp = dx / dy;
                dx = sqrt(1.0 + tmp*tmp);
                s = Sign(1.0 / dx, dy);
                c = tmp*s;
            }
            else if (fabs(dy) <= fabs(dx)) 
            {
                double tmp = dy / dx;
                dy = sqrt(1.0 + tmp*tmp);
                c = Sign(1.0 / dy, dx);
                s = tmp*c;
            }
            else 
            {
                // dx and/or dy must be invalid
                dx = 0.0;
                dy = 0.0;
                c = 1.0;
                s = 0.0;
            }
            dx = fabs(dx*dy);
            dy = 0.0;
        }

        void MATH_LinearSolver::SolveReduced(const int & n, const std::vector<std::vector<double> > & Hsbg, const std::vector<double> & rhs, std::vector<double> & x)
        {
            // initialize...
            for (int i = 0; i < n; i++)
                x[i] = rhs[i];
            // ... and backsolve
            for (int i = n - 1; i >= 0; i--) 
            {
                x[i] /= Hsbg[i][i];
                for (int j = i - 1; j >= 0; j--) 
                {
                    x[j] -= Hsbg[j][i] * x[i];
                }
            }
        }

        void MATH_LinearSolver::ModGramSchmidt(int i, std::vector<std::vector<double> > & Hsbg, std::vector<MATH_Vector> & w) 
        {
            bool Convergence = true;
            int rank = TBOX::MASTER_NODE;

#ifdef HAVE_MPI
            int size;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

            /*--- Parameter for reorthonormalization ---*/
            static const double reorth = 0.98;

            /*--- Get the norm of the vector being orthogonalized, and find the
            threshold for re-orthogonalization ---*/

            double nrm = dotProd(w[i + 1], w[i + 1]);
            double thr = nrm*reorth;

            /*--- The norm of w[i+1] < 0.0 or w[i+1] = NaN ---*/
            if ((nrm <= 0.0) || (nrm != nrm)) Convergence = false;

            /*--- Synchronization point to check the convergence of the solver ---*/

#ifdef HAVE_MPI

            unsigned short *sbuf_conv = NULL, *rbuf_conv = NULL;
            sbuf_conv = new unsigned short[1]; sbuf_conv[0] = 0;
            rbuf_conv = new unsigned short[1]; rbuf_conv[0] = 0;

            /*--- Convergence criteria ---*/

            sbuf_conv[0] = Convergence;
            MPI_Reduce(sbuf_conv, rbuf_conv, 1, MPI_UNSIGNED_SHORT, MPI_SUM, TBOX::MASTER_NODE, MPI_COMM_WORLD);

            /*-- Compute global convergence criteria in the master node --*/

            sbuf_conv[0] = 0;
            if (rank == TBOX::MASTER_NODE) {
                if (rbuf_conv[0] == size) sbuf_conv[0] = 1;
                else sbuf_conv[0] = 0;
            }

            MPI_Bcast(sbuf_conv, 1, MPI_UNSIGNED_SHORT, TBOX::MASTER_NODE, MPI_COMM_WORLD);

            if (sbuf_conv[0] == 1) Convergence = true;
            else Convergence = false;

            delete [] sbuf_conv;
            delete [] rbuf_conv;

#endif

            if (!Convergence) 
            {
                if (rank == TBOX::MASTER_NODE)
                    std::cout << "\n !!! Error: SU2 has diverged. Now exiting... !!! \n" << std::endl;
#ifndef HAVE_MPI
                exit(TBOX::EXIT_DIVERGENCE);
#else
                MPI_Abort(MPI_COMM_WORLD,1);
#endif
            }

            /*--- Begin main Gram-Schmidt loop ---*/
            for (int k = 0; k < i + 1; k++) 
            {
                double prod = dotProd(w[i + 1], w[k]);
                Hsbg[k][i] = prod;
                w[i + 1].Plus_AX(-prod, w[k]);

                /*--- Check if reorthogonalization is necessary ---*/
                if (prod*prod > thr) 
                {
                    prod = dotProd(w[i + 1], w[k]);
                    Hsbg[k][i] += prod;
                    w[i + 1].Plus_AX(-prod, w[k]);
                }

                /*--- Update the norm and check its size ---*/
                nrm -= Hsbg[k][i] * Hsbg[k][i];
                if (nrm < 0.0) nrm = 0.0;
                thr = nrm*reorth;
            }

            /*--- Test the resulting vector ---*/

            nrm = w[i + 1].norm();
            Hsbg[i + 1][i] = nrm;

            //  if (nrm <= 0.0) {
            //    
            //    /*--- w[i+1] is a linear combination of the w[0:i] ---*/
            //    
            //    cerr << "The FGMRES linear solver has diverged" << std::endl;
            //#ifndef HAVE_MPI
            //    exit(EXIT_DIVERGENCE);
            //#else
            //    MPI_Abort(MPI_COMM_WORLD,1);
            //    MPI_Finalize();
            //#endif
            //    
            //  }

            /*--- Scale the resulting vector ---*/
            w[i + 1] /= nrm;
        }

        void MATH_LinearSolver::WriteHeader(const std::string & solver, const double & restol, const double & resinit) {

            std::cout << "\n# " << solver << " residual history" << std::endl;
            std::cout << "# Residual tolerance target = " << restol << std::endl;
            std::cout << "# Initial residual norm     = " << resinit << std::endl;

        }

        void MATH_LinearSolver::WriteHistory(const int & iter, const double & res, const double & resinit) {

            std::cout << "     " << iter << "     " << res / resinit << std::endl;

        }

        unsigned long MATH_LinearSolver::CG_LinSolver(const MATH_Vector & b, MATH_Vector & x, MATH_MatrixVectorProduct & mat_vec,
            MATH_Preconditioner & precond, double tol, unsigned long m, bool monitoring) 
        {
            int rank = 0;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Check the subspace size ---*/
            if (m < 1) 
            {
                if (rank == TBOX::MASTER_NODE) 
                    std::cerr << "MATH_LinearSolver::ConjugateGradient: illegal value for subspace size, m = " << m << std::endl;
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD,1);
                MPI_Finalize();
#endif
            }

            MATH_Vector r(b);
            MATH_Vector A_p(b);

            /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/
            mat_vec(x, A_p);

            r -= A_p; // recall, r holds b initially
            double norm_r = r.norm();
            double norm0 = b.norm();
            if ((norm_r < tol*norm0) || (norm_r < eps))
            {
                if (rank == TBOX::MASTER_NODE) std::cout << "MATH_LinearSolver::ConjugateGradient(): system solved by initial guess." << std::endl;
                return 0;
            }

            double alpha, beta, r_dot_z;
            MATH_Vector z(r);
            precond(r, z);
            MATH_Vector p(z);

            /*--- Set the norm to the initial initial residual value ---*/
            norm0 = norm_r;

            /*--- Output header information including initial residual ---*/
            int i = 0;
            if ((monitoring) && (rank == TBOX::MASTER_NODE)) 
            {
                WriteHeader("CG", tol, norm_r);
                WriteHistory(i, norm_r, norm0);
            }

            /*---  Loop over all search directions ---*/
            for (i = 0; i < m; i++) 
            {
                /*--- Apply matrix to p to build Krylov subspace ---*/
                mat_vec(p, A_p);

                /*--- Calculate step-length alpha ---*/
                r_dot_z = dotProd(r, z);
                alpha = dotProd(A_p, p);
                alpha = r_dot_z / alpha;

                /*--- Update solution and residual: ---*/
                x.Plus_AX(alpha, p);
                r.Plus_AX(-alpha, A_p);

                /*--- Check if solution has converged, else output the relative residual if necessary ---*/
                norm_r = r.norm();
                if (norm_r < tol*norm0) break;
                if (((monitoring) && (rank == TBOX::MASTER_NODE)) && ((i + 1) % 5 == 0)) WriteHistory(i + 1, norm_r, norm0);

                precond(r, z);

                /*--- Calculate Gram-Schmidt coefficient beta,
                     beta = dotProd(r_{i+1}, z_{i+1}) / dotProd(r_{i}, z_{i}) ---*/
                beta = 1.0 / r_dot_z;
                r_dot_z = dotProd(r, z);
                beta *= r_dot_z;

                /*--- Gram-Schmidt orthogonalization; p = beta *p + z ---*/
                p.Equals_AX_Plus_BY(beta, p, 1.0, z);
            }



            if ((monitoring) && (rank == TBOX::MASTER_NODE)) 
            {
                std::cout << "# Conjugate Gradient final (true) residual:" << std::endl;
                std::cout << "# Iteration = " << i << ": |res|/|res0| = " << norm_r / norm0 << ".\n" << std::endl;
            }

            //  /*--- Recalculate final residual (this should be optional) ---*/
            //  mat_vec(x, A_p);
            //  r = b;
            //  r -= A_p;
            //  double true_res = r.norm();
            //  
            //  if (fabs(true_res - norm_r) > tol*10.0) {
            //    if (rank == TBOX::MASTER_NODE) {
            //      std::cout << "# WARNING in MATH_LinearSolver::ConjugateGradient(): " << std::endl;
            //      std::cout << "# true residual norm and calculated residual norm do not agree." << std::endl;
            //      std::cout << "# true_res - calc_res = " << true_res - norm_r << std::endl;
            //    }
            //  }

            return i;

        }

        unsigned long MATH_LinearSolver::FGMRES_LinSolver(const MATH_Vector & b, MATH_Vector & x, MATH_MatrixVectorProduct & mat_vec,
            MATH_Preconditioner & precond, double tol, unsigned long m, double *residual, bool monitoring) 
        {
            int rank = 0;

#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*---  Check the subspace size ---*/
            if (m < 1)
            {
                if (rank == TBOX::MASTER_NODE) std::cerr << "MATH_LinearSolver::FGMRES: illegal value for subspace size, m = " << m << std::endl;
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD,1);
                MPI_Finalize();
#endif
            }

            /*---  Check the subspace size ---*/

            if (m > 1000) 
            {
                if (rank == TBOX::MASTER_NODE) std::cerr << "MATH_LinearSolver::FGMRES: illegal value for subspace size (too high), m = " << m << std::endl;
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD,1);
                MPI_Finalize();
#endif
            }

            /*---  Define various arrays
               Note: elements in w and z are initialized to x to avoid creating
               a temporary MATH_Vector object for the copy constructor ---*/

            std::vector<MATH_Vector> w(m + 1, x);
            std::vector<MATH_Vector> z(m + 1, x);
            std::vector<double> g(m + 1, 0.0);
            std::vector<double> sn(m + 1, 0.0);
            std::vector<double> cs(m + 1, 0.0);
            std::vector<double> y(m, 0.0);
            std::vector<std::vector<double> > H(m + 1, std::vector<double>(m, 0.0));

            /*---  Calculate the norm of the rhs vector ---*/

            double norm0 = b.norm();

            /*---  Calculate the initial residual (actually the negative residual)
               and compute its norm ---*/

            mat_vec(x, w[0]);
            w[0] -= b;

            double beta = w[0].norm();

            if ((beta < tol*norm0) || (beta < eps)) 
            {
                /*---  System is already solved ---*/
                if (rank == TBOX::MASTER_NODE) std::cout << "MATH_LinearSolver::FGMRES(): system solved by initial guess." << std::endl;
                return 0;
            }

            /*---  Normalize residual to get w_{0} (the negative sign is because w[0]
               holds the negative residual, as mentioned above) ---*/
            w[0] /= -beta;

            /*---  Initialize the RHS of the reduced system ---*/
            g[0] = beta;

            /*--- Set the norm to the initial residual value ---*/

            norm0 = beta;

            /*---  Output header information including initial residual ---*/

            int i = 0;
            if ((monitoring) && (rank == TBOX::MASTER_NODE)) 
            {
                WriteHeader("FGMRES", tol, beta);
                WriteHistory(i, beta, norm0);
            }

            /*---  Loop over all search directions ---*/

            for (i = 0; i < m; i++) 
            {

                /*---  Check if solution has converged ---*/

                if (beta < tol*norm0) break;

                /*---  Precondition the MATH_Vector w[i] and store result in z[i] ---*/

                precond(w[i], z[i]);

                /*---  Add to Krylov subspace ---*/

                mat_vec(z[i], w[i + 1]);

                /*---  Modified Gram-Schmidt orthogonalization ---*/

                ModGramSchmidt(i, H, w);

                /*---  Apply old Givens rotations to new column of the Hessenberg matrix
                     then generate the new Givens rotation matrix and apply it to
                     the last two elements of H[:][i] and g ---*/

                for (int k = 0; k < i; k++)
                    ApplyGivens(sn[k], cs[k], H[k][i], H[k + 1][i]);
                GenerateGivens(H[i][i], H[i + 1][i], sn[i], cs[i]);
                ApplyGivens(sn[i], cs[i], g[i], g[i + 1]);

                /*---  Set L2 norm of residual and check if solution has converged ---*/

                beta = fabs(g[i + 1]);

                /*---  Output the relative residual if necessary ---*/

                if ((((monitoring) && (rank == TBOX::MASTER_NODE)) && ((i + 1) % 50 == 0)) && (rank == TBOX::MASTER_NODE)) WriteHistory(i + 1, beta, norm0);

            }

            /*---  Solve the least-squares system and update solution ---*/

            SolveReduced(i, H, g, y);
            for (int k = 0; k < i; k++) {
                x.Plus_AX(y[k], z[k]);
            }

            if ((monitoring) && (rank == TBOX::MASTER_NODE)) {
                std::cout << "# FGMRES final (true) residual:" << std::endl;
                std::cout << "# Iteration = " << i << ": |res|/|res0| = " << beta / norm0 << ".\n" << std::endl;
            }

            //  /*---  Recalculate final (neg.) residual (this should be optional) ---*/
            //  mat_vec(x, w[0]);
            //  w[0] -= b;
            //  double res = w[0].norm();
            //
            //  if (fabs(res - beta) > tol*10) {
            //    if (rank == TBOX::MASTER_NODE) {
            //      std::cout << "# WARNING in MATH_LinearSolver::FGMRES(): " << std::endl;
            //      std::cout << "# true residual norm and calculated residual norm do not agree." << std::endl;
            //      std::cout << "# res - beta = " << res - beta << std::endl;
            //    }
            //  }

            (*residual) = beta;
            return i;

        }

        unsigned long MATH_LinearSolver::BCGSTAB_LinSolver(const MATH_Vector & b, MATH_Vector & x, MATH_MatrixVectorProduct & mat_vec,
            MATH_Preconditioner & precond, double tol, unsigned long m, double *residual, bool monitoring) {

            int rank = 0;
#ifdef HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

            /*--- Check the subspace size ---*/

            if (m < 1) {
                if (rank == TBOX::MASTER_NODE) std::cerr << "MATH_LinearSolver::BCGSTAB: illegal value for subspace size, m = " << m << std::endl;
#ifndef HAVE_MPI
                exit(EXIT_FAILURE);
#else
                MPI_Abort(MPI_COMM_WORLD,1);
                MPI_Finalize();
#endif
            }

            MATH_Vector r(b);
            MATH_Vector r_0(b);
            MATH_Vector p(b);
            MATH_Vector v(b);
            MATH_Vector s(b);
            MATH_Vector t(b);
            MATH_Vector phat(b);
            MATH_Vector shat(b);
            MATH_Vector A_x(b);

            /*--- Calculate the initial residual, compute norm, and check if system is already solved ---*/

            mat_vec(x, A_x);
            r -= A_x; r_0 = r; // recall, r holds b initially
            double norm_r = r.norm();
            double norm0 = b.norm();
            if ((norm_r < tol*norm0) || (norm_r < eps)) {
                if (rank == TBOX::MASTER_NODE) std::cout << "MATH_LinearSolver::BCGSTAB(): system solved by initial guess." << std::endl;
                return 0;
            }

            /*--- Initialization ---*/

            double alpha = 1.0, beta = 1.0, omega = 1.0, rho = 1.0, rho_prime = 1.0;

            /*--- Set the norm to the initial initial residual value ---*/

            norm0 = norm_r;

            /*--- Output header information including initial residual ---*/

            int i = 0;
            if ((monitoring) && (rank == TBOX::MASTER_NODE)) {
                WriteHeader("BCGSTAB", tol, norm_r);
                WriteHistory(i, norm_r, norm0);
            }

            /*---  Loop over all search directions ---*/

            for (i = 0; i < m; i++) {

                /*--- Compute rho_prime ---*/

                rho_prime = rho;

                /*--- Compute rho_i ---*/

                rho = dotProd(r, r_0);

                /*--- Compute beta ---*/

                beta = (rho / rho_prime) * (alpha / omega);

                /*--- p_{i} = r_{i-1} + beta * p_{i-1} - beta * omega * v_{i-1} ---*/

                double beta_omega = -beta*omega;
                p.Equals_AX_Plus_BY(beta, p, beta_omega, v);
                p.Plus_AX(1.0, r);

                /*--- Preconditioning step ---*/

                precond(p, phat);
                mat_vec(phat, v);

                /*--- Calculate step-length alpha ---*/

                double r_0_v = dotProd(r_0, v);
                alpha = rho / r_0_v;

                /*--- s_{i} = r_{i-1} - alpha * v_{i} ---*/

                s.Equals_AX_Plus_BY(1.0, r, -alpha, v);

                /*--- Preconditioning step ---*/

                precond(s, shat);
                mat_vec(shat, t);

                /*--- Calculate step-length omega ---*/

                omega = dotProd(t, s) / dotProd(t, t);

                /*--- Update solution and residual: ---*/

                x.Plus_AX(alpha, phat); x.Plus_AX(omega, shat);
                r.Equals_AX_Plus_BY(1.0, s, -omega, t);

                /*--- Check if solution has converged, else output the relative residual if necessary ---*/

                norm_r = r.norm();
                if (norm_r < tol*norm0) break;
                if (((monitoring) && (rank == TBOX::MASTER_NODE)) && ((i + 1) % 50 == 0) && (rank == TBOX::MASTER_NODE)) WriteHistory(i + 1, norm_r, norm0);

            }

            if ((monitoring) && (rank == TBOX::MASTER_NODE)) {
                std::cout << "# BCGSTAB final (true) residual:" << std::endl;
                std::cout << "# Iteration = " << i << ": |res|/|res0| = " << norm_r / norm0 << ".\n" << std::endl;
            }

            //  /*--- Recalculate final residual (this should be optional) ---*/
            //	mat_vec(x, A_x);
            //  r = b; r -= A_x;
            //  double true_res = r.norm();
            //  
            //  if ((fabs(true_res - norm_r) > tol*10.0) && (rank == TBOX::MASTER_NODE)) {
            //    std::cout << "# WARNING in MATH_LinearSolver::BCGSTAB(): " << std::endl;
            //    std::cout << "# true residual norm and calculated residual norm do not agree." << std::endl;
            //    std::cout << "# true_res - calc_res = " << true_res <<" "<< norm_r << std::endl;
            //  }

            (*residual) = norm_r;
            return i;
        }

        unsigned long MATH_LinearSolver::Solve(MATH_Matrix & Jacobian, MATH_Vector & LinSysRes, MATH_Vector & LinSysSol, GEOM::GEOM_Geometry *geometry, TBOX::TBOX_Config *config) 
        {

            double SolverTol = config->GetLinear_Solver_Error(), Residual;
            unsigned long MaxIter = config->GetLinear_Solver_Iter();
            unsigned long IterLinSol = 0;

            /*--- Solve the linear system using a Krylov subspace method ---*/

            if (config->GetKind_Linear_Solver() == TBOX::BCGSTAB || config->GetKind_Linear_Solver() == TBOX::FGMRES
                || config->GetKind_Linear_Solver() == TBOX::RESTARTED_FGMRES) 
            {

                MATH_MatrixVectorProduct* mat_vec = new MATH_Matrix_MatrixVectorProduct(Jacobian, geometry, config);

                MATH_Preconditioner* precond = NULL;

                switch (config->GetKind_Linear_Solver_Prec()) {
                case TBOX::JACOBI:
                    Jacobian.BuildJacobiPreconditioner();
                    precond = new MATH_JacobiPreconditioner(Jacobian, geometry, config);
                    break;
                case TBOX::ILU:
                    Jacobian.BuildILUPreconditioner();
                    precond = new MATH_ILUPreconditioner(Jacobian, geometry, config);
                    break;
                case TBOX::LU_SGS:
                    precond = new MATH_LUSGSPreconditioner(Jacobian, geometry, config);
                    break;
                case TBOX::LINELET:
                    Jacobian.BuildJacobiPreconditioner();
                    precond = new MATH_LineletPreconditioner(Jacobian, geometry, config);
                    break;
                default:
                    Jacobian.BuildJacobiPreconditioner();
                    precond = new MATH_JacobiPreconditioner(Jacobian, geometry, config);
                    break;
                }

                switch (config->GetKind_Linear_Solver()) {
                case TBOX::BCGSTAB:
                    IterLinSol = BCGSTAB_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
                    break;
                case TBOX::FGMRES:
                    IterLinSol = FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
                    break;
                case TBOX::RESTARTED_FGMRES:
                    IterLinSol = 0;
                    while (IterLinSol < config->GetLinear_Solver_Iter()) {
                        if (IterLinSol + config->GetLinear_Solver_Restart_Frequency() > config->GetLinear_Solver_Iter())
                            MaxIter = config->GetLinear_Solver_Iter() - IterLinSol;
                        IterLinSol += FGMRES_LinSolver(LinSysRes, LinSysSol, *mat_vec, *precond, SolverTol, MaxIter, &Residual, false);
                        if (LinSysRes.norm() < SolverTol) break;
                        SolverTol = SolverTol*(1.0 / LinSysRes.norm());
                    }
                    break;
                }

                /*--- Dealocate memory of the Krylov subspace method ---*/

                delete mat_vec;
                delete precond;

            }

            /*--- Smooth the linear system. ---*/

            else {
                switch (config->GetKind_Linear_Solver()) {
                case TBOX::SMOOTHER_LUSGS:
                    Jacobian.ComputeLU_SGSPreconditioner(LinSysRes, LinSysSol, geometry, config);
                    break;
                case TBOX::SMOOTHER_JACOBI:
                    Jacobian.BuildJacobiPreconditioner();
                    Jacobian.ComputeJacobiPreconditioner(LinSysRes, LinSysSol, geometry, config);
                    break;
                case TBOX::SMOOTHER_ILU:
                    Jacobian.BuildILUPreconditioner();
                    Jacobian.ComputeILUPreconditioner(LinSysRes, LinSysSol, geometry, config);
                    break;
                case TBOX::SMOOTHER_LINELET:
                    Jacobian.BuildJacobiPreconditioner();
                    Jacobian.ComputeLineletPreconditioner(LinSysRes, LinSysSol, geometry, config);
                    break;
                    IterLinSol = 1;
                }
            }

            return IterLinSol;

        }
    }
}