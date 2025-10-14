/**
 * @file IdealGasConstPressureAdiabaticReactorAdapter.h
 * @brief Adapter that exposes IdealGasConstPressureAdiabaticReactor via the Utility interface.
 * @details
 *   - Delegates Utility API calls to an underlying IdealGasConstPressureAdiabaticReactor instance.
 *   - Keeps ownership external (stores a reference).
 */

#ifndef SRC_INCLUDE_ADAPTERS_IDEAL_GAS_CONST_PRESSURE_ADIABATIC_REACTOR_H
#define SRC_INCLUDE_ADAPTERS_IDEAL_GAS_CONST_PRESSURE_ADIABATIC_REACTOR_H

#include <iostream>

/* Headers */
#include "Utility.h"
#include "IdealGasConstPressureAdiabaticReactor.h"


/**
 * @class IdealGasConstPressureAdiabaticReactorAdapter
 * @brief Utility-interface adapter for IdealGasConstPressureAdiabaticReactor.
 * @details
 *   TODO: Add usage notes.
 */
class IdealGasConstPressureAdiabaticReactorAdapter : public Utility
{
    public:
        /**
         * @brief Construct the adapter around an existing reactor.
         * @param[in] r Reactor instance to be adapted (non-owning reference).
         * @param[in] (implicit) debug Prints a message when nonzero (see @ref debug_).
         * @pre @p r must outlive this adapter.
         * @post Adapter is ready to forward Utility calls to @p r.
         */
        explicit IdealGasConstPressureAdiabaticReactorAdapter
        (IdealGasConstPressureAdiabaticReactor &r) : r_(r) 
        {
            if(debug_ == 1)
            {
                std::cout<<"--Reactor Adapter constructor implemented!"<<std::endl;
            }
        }
        
        /**
         * @brief Set the number of ODE equations.
         * @return Dimension of the system as provided by the reactor.
         */
        int setNEQ() override
        {
            return r_.setNEQ();
        }

        /**
         * @brief Write initial state into the provided buffer.
         * @param[out] y State vector (length = setNEQ()).
         * @post @p y is filled with the reactor’s initial condition.
         */
        void setInitialState(double *y) override
        {
            r_.setInitialState(y);
        }

        /**
         * @brief Evaluate the ODE right-hand side.
         * @param[in] t [What is this param(?)].
         * @param[in] y State vector at time @p t (length = setNEQ()).
         * @param[out] ydot Derivative vector (length = setNEQ()).
         * @post @p ydot contains dy/dt evaluated by the reactor.
         */
        void evalRHS(double t, double *y, double *ydot) override
        {
            r_.evalRHS(t, y, ydot);
        } 

    private:
        IdealGasConstPressureAdiabaticReactor &r_;                              ///< Non-owning reference to the wrapped reactor.
        int debug_ = 0;                                                         ///< Debug verbosity: 0 = quiet, 1 = prints constructor msg.
};


#endif
