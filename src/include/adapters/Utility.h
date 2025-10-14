/**
 * @file Utility.h
 * @brief Abstract interface for ODE models used by integrators.
 * @details
 *   Threading/ownership notes, units, and error-handling conventions can go here.
 */

#ifndef SRC_INCLUDE_ADAPTERS_UTILITY_H
#define SRC_INCLUDE_ADAPTERS_UTILITY_H


/* Headers? */


/* Virtual fns because I can change their dfns in the derived class */
/**
 * @class Utility
 * @brief Minimal abstract interface an ODE model must implement.
 * @details
 *   Derive concrete models (reactors, mechanical systems, etc.) from this interface.
 */
class Utility
{
    public:
        virtual ~Utility() = default;
        
        /* CVODES fns */

        /**
         * @brief Return the number of ODE equations.
         * @return System dimension (length of y and ydot).
         */
        /* Reactor fns */
        virtual int setNEQ() = 0;

        /**
         * @brief Populate an external pointer(?) with the initial state.
         * @param[out] y Caller-allocated array of length setNEQ().
         * @post @p y contains the initial condition.
         */
        virtual void setInitialState(double *y) = 0;

        /**
         * @brief Evaluate the ODE right-hand side.
         * @param[in] t [What is this param(?)] 
         * @param[in] y State vector at time @p t (length = setNEQ()).
         * @param[out] ydot Derivative vector (length = setNEQ()).
         * @post @p ydot is filled with dy/dt at (t, y).
         * @note Implementations should avoid allocations inside this call for performance.
         */
        virtual void evalRHS(double t, double *y, double *ydot) = 0;            

        /* virtual void setState(double *y, double temperature) = 0; */

        /* CVODES fns */

        /* Other fns */


};


#endif
