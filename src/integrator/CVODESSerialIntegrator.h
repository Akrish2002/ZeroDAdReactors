/** 
 * @file CVODESSerialIntegrator.h
 * @brief Serial CVODES integrator wrapper for ODE systems
 * @details
 *  @todo Implementation has to be corrected since it takes Temp + NSpecies ins
 *  tead of Temp + NSpecies - 1
 */

#ifndef SRC_INTEGRATOR_CVODES_SERIAL_INTEGRATOR
#define SRC_INTEGRATOR_CVODES_SERIAL_INTEGRATOR

#include <stdio.h>
#include <iostream>
#include <array>
#include <vector>

/** @name SUNDIALS includes
 *  @brief Public headers for CVODE, N_Vector (serial), dense linear solvers and matrices.
 *  @{ */
#include <cvodes/cvodes.h>          						                    /*!< Prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h> 						                    /*!< Access to serial N_Vector            */
#include <sunlinsol/sunlinsol_dense.h> 						                    /*!< Access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> 						                    /*!< Access to dense SUNMatrix            */
/** @} */

#include "Utility.h"                                                            /*!< Model interface to provide setNEQ(), setInitialState(), evalRHS(), etc. */


/**
 * @class CVODESSerialIntegrator
 * @brief Wrapper for a serial CVODES integration session.
 * @details
 *  TODO: 
 *  Expand on:
 *      - Error handling strategy and debug modes.
 *      - Step/output scheduling via @ref timeop_ and @ref TMULT_.
 */
class CVODESSerialIntegrator
{
    public:
        /**
         * @brief Construct an integrator bound to a model.
         * @param[in] model Reference to utility object.
         * @param[in] debug Debug flag (0 = quiet, 1 = verbose).
         * @post Internals are default-initialized; allocate/setup occurs in initializeandsetupsolver().
         */
        explicit CVODESSerialIntegrator (Utility &model, int debug = 0);
        /**
         * @brief Allocate objects and configure solver/integrator.
         */
        void initializeandsetupsolver();
        /**
         * @brief Run time integration with CVODES and print stats/output.
         * @details
         *   TODO: Describe output cadence (timeop_ × TMULT_ for steps_ iterations) and exit conditions.
         * @post Integration performed; statistics printed; file closed.
         */
        void integrate();
        /**
         * @brief Release all allocated SUNDIALS resources and files.
         * @post All owned resources freed safely (idempotent if called once).
         */
        void freeMemory();

        /* Debugger, getter fns */
        /**
         * @brief Get number of ODE equations.
         * @return NEQ_.
         */
        int getNEQ();
        /**
         * @brief Get the first state component y[0] (if available).
         * @return y[0].
         * @note Provided for quick inspection during debugging.
         */
        double getzeroEqn();
        /**
         * @brief Get the second state component y[1] (if available).
         * @return y[1].
         * @note Provided for quick inspection during debugging.
         */
        double getfirstEqn();


    private:

        int    NEQ_;                                                            /*!< Number of equations (model dimension).                     */
        double time0_;                                                          /*!< Initial time [s].                                          */
        double timeop_;                                                         /*!< Output time for first report [s].                          */
        double time_;                                                           /*!< Current integrator time [s].                               */
        double timefinal_;                                                      /*!< Final time [s].                                            */
        double TMULT_;                                                          /*!< Multiplicative factor for successive output times.         */
        double timeInteg_;                                                      /*!< Wall/CPU integration time (if measured).                   */
        int    steps_;                                                          /*!< Number of output steps (loop count).                       */    

        
        void*     cvode_mem_;                                                   /*!< CVODES memory block pointer (session handle). */

        N_Vector  abstol_;                                                      /*!< Absolute tolerance vector (per-equation). */
        double    RTOL_;                                                        /*!< Relative tolerance (scalar). */
        std::vector<double> ATOL_;                                              /*!< Host-side copy for absolute tolerances. */

        N_Vector  y_;                                                           /*!< State vector. */
        SUNMatrix A_;                                                           /*!< Dense Jacobian/SUNMatrix. */
        SUNLinearSolver LS_;                                                    /*!< Dense linear solver.  */
        SUNContext sunctx_;                                                     /*!< SUNDIALS context (logs/errors/profiling). */


        Utility  &model_;                                                       /*!< Reference to user model providing RHS and initial state. */
        FILE*     FID_;                                                         /*!< File handle for printing solver stats (CSV). */

        int debug_;


        /* ---------------------
        * Function declarations 
        * ---------------------- */ 
        /* 1.1. Initialization */
        /* 1.2. Allocate and create memory for solvers, integrators and matrices */
        /* 1.3. Open file for printing statistics */

        /**
         * @brief Allocate N_Vectors and create SUNContext.
         * @return 0 on success, non-zero on failure.
         */
        int allocateMemory();

        /**
         * @brief Populate @ref y_ with the model's initial state.
         * @post @ref y_ contains initial condition from @ref Utility::setInitialState().
         */
        void setInitialState();

        /**
         * @brief Populate absolute tolerance vector from @ref ATOL_.
         * @post @ref abstol_ is filled with element-wise absolute tolerances.
         */
        void setTolerances();

        /**
         * @brief Create a CVODES session and choose method.
         * @param[in] method CV_ADAMS or CV_BDF (integer tag).
         * @post @ref cvode_mem_ non-null on success.
         * @warning Caller must check for errors if extending this method.
         */
        void allocatesolverMemoryandMethod(int method);

        /**
         * @brief Initialize CVODE with RHS, initial time, and state; attach user data.
         * @post RHS set; initial conditions registered.
         */
        void initializeintegratorMemoryandRHS();

        /**
         * @brief Apply scalar-relative, vector-absolute tolerances to session.
         * @post Tolerances active in CVODES run.
         */
        void setrelativeTolerance();

        /**
         * @brief Create dense SUNMatrix sized NEQ_×NEQ_.
         * @post @ref A_ valid on success.
         */
        void createSUNDenseMatrix();

        /**
         * @brief Create dense SUNLinearSolver for (@ref y_, @ref A_).
         * @post @ref LS_ valid on success.
         */
        void createSUNLinSolObject(); 

        /**
         * @brief Attach matrix and linear solver to CVODES session.
         * @post Linear solver configured for the session.
         */
        void attachMatrixandLinSol();

        /**
         * @brief Open CSV file for printing solver statistics.
         * @post @ref FID_ opened for write.
         */
        void openfileforprinting();


        /* 2. Integrator fns */
        /* int cvode_rhs(double t, N_Vector y, 
         * N_Vector ydot, void* user_data); */                                  /* Cannot have it as a member function */

        /* 3. Free memory */
        /**
         * @brief Destroy owned N_Vectors.
         */
        void destroyN_Vectors();

        /**
         * @brief Free CVODES session memory.
         */
        void freeblockmemory();

        /**
         * @brief Free linear solver memory.
         */
        void freesolvermemory();

        /**
         * @brief Destroy SUNMatrix.
         */
        void freematrix();

        /**
         * @brief Free SUNDIALS context.
         */
        void freeSUNDIALScontext();

        /**
         * @brief Close statistics/output file if open.
         */
        void closefile();
                                                                                
};


#endif /* SRC_INTEGRATOR_CVODES_SERIAL_INTEGRATOR */
