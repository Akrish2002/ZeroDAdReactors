/**
 * @file IdealGasConstPressureAdiabaticReactor.h
 * @brief Module for implementing constant-pressure, adiabatic reactor with ideal-gas law.
 * @details
 *   - Constant pressure
 *   - Adiabatic (no heat exchange)
 *   - Ideal-gas mixture
 */

#ifndef SRC_REACTOR_IDEAL_GAS_CONST_PRESSURE_ADIABATIC_REACTOR
#define SRC_REACTOR_IDEAL_GAS_CONST_PRESSURE_ADIABATIC_REACTOR

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <array>
#include <algorithm>
#include <cstring>

/* Chemgen header files */
#include "types_inl.h"  /* For Species */

/**
 * @class IdealGasConstPressureAdiabaticReactor
 * @brief Constant-pressure, adiabatic reactor model for an ideal-gas mixture.
 * @details
 *   @todo Add details on the state vector layout, rate/thermo sources, and integration expectations.
 */
class IdealGasConstPressureAdiabaticReactor
{
    public:
        /**
         * @brief Construct a 0D-CP reactor for a simple compressible substance.
         * @param[in] n Number of chemical species in the mixture.
         * @param[in] temperature Initial temperature [K]. (default: 300)
         * @param[in] pressure Initial pressure [Pa]. (default: 101325.0)
         * @throws (none) TODO: Specify if implementation throws.
         */
        IdealGasConstPressureAdiabaticReactor(int n, double temperature = 300, double pressure = 101325.0);
    
        /**
         * @brief Refresh thermodynamic and kinetic properties for current state.
         * @details
         *   @todo Enumerate which fields are updated (cp_, h_, cp_bar_, h_bar_, omega_, etc.).
         * @post Properties consistent with current @ref T_, @ref Y_, and @ref C_.
         */
        void getProperties();
    
        /**
         * @brief Compute numerator helper for dT/dt.
         * @details
         * @note  Currently unused.
         */
        void getdTtNumerator();
    
        /**
         * @brief Compute denominator term appearing in dT/dt expression.
         * @details
         * @note  Currently unused
         */
        void getdTtDenominator();
    
        /**
         * @brief Return number of ODE equations for the reactor system.
         * @return Number of equations (state dimension).
         * @note Typically 1 (temperature) + number of species for mass fractions.
         */
        int setNEQ();
    
        /**
         * @brief Populate an external state vector with the initial condition.
         * @param[out] y State array of length setNEQ(); layout [T, Y1..Y_N].
         * @pre @p y is valid and has length = setNEQ().
         * @post @p y contains initial states of @ref T_ and @ref Y_.
         */
        void setInitialState(double* y);
    
        /**
         * @brief Replace internal state from an external vector and temperature.
         * @param[in] y State vector [T, Y1..Y_N]; only mass fractions are used here.
         * @param[in] temperature Temperature override [K].
         * @details
         *   @todo Clarify consistency policy between @p y[0] and @p temperature.
         * @post Internal @ref T_, @ref Y_, and @ref C_ updated.
         */
        void setState(double* y, double temperature);
    
        /**
         * @brief Evaluate ODE right-hand side (RHS): dT/dt and dY/dt.
         * @param[in] t What is this param for? (Unused)
         * @param[in] y State vector at time @p t; layout [T, Y1..Y_N].
         * @param[out] ydot Derivative vector; layout [dT/dt, dY1/dt..dY_N/dt].
         * @details
         *   @todo Document formulas, assumptions (ideal gas, constant P).
         * @pre @p y and @p ydot have length setNEQ().
         * @post @p ydot populated for integrator.
         */
        void evalRHS(double t, double* y, double* ydot);
    
    
        /* ---------------- Debug/Misc accessors ---------------- */
    
        /**
         * @brief Get current temperature.
         * @return Temperature [K].
         */
        double getTemperature();
    
        /**
         * @brief Get current pressure.
         * @return Pressure [Pa].
         */
        double getPressure();
    
        /**
         * @brief Get number of species in the mixture.
         * @return Number of species.
         */
        double getNumberofSpecies();
    
        /**
         * @brief Get molecular weights vector.
         * @details
         *   @todo Document units and whether this triggers a recompute or returns cached data.
         */
        void getMW();
    
        /**
         * @brief Get heat capacities of all species.
         */
        void getcp();
    
        /**
         * @brief Get production rates for each species.
         * @details
         *   @todo Specify units and sign convention.
         */
        void getomega();
    
    
    private:
        /* ---------------- Mixture totals / sizes ---------------- */
    
        double Ctot_;        ///< Total molar concentration [kmol/m^3 or mol/m^3]. @todo Confirm units.
        int    n_species_;   ///< Number of species in the mechanism.
    
        /* ---------------- Thermodynamic state ------------------ */
    
        double T_;           ///< Temperature [K].
        double P_;           ///< Pressure [Pa].
    
        double  MWtot_;      ///< Mixture-averaged molecular weight.
        Species MW_;         ///< Species molecular weights (length = n_species_).
    
        Species cp_;         ///< Mass specific constant-pressure specific heats [J/(kg·K)].
        Species cp_bar_;     ///< Molar specifc constant-pressure specific heats [J/(mol·K)]. 
        Species h_;          ///< Mass specific enthalpies [J/kg].
        Species h_bar_;      ///< Molar specific enthalpies [J/kmol].
    
        /* ---------------- Composition / concentrations --------- */
    
        Species X_;          ///< Mole fractions (sum to 1).        Length = n_species_.
        Species Y0_;         ///< Initial mass fractions.           Length = n_species_.
        Species Y_;          ///< Current mass fractions.           Length = n_species_.
        Species C_;          ///< Species concentrations [mol/m^3]. Length = n_species_.
    
        /* ---------------- Kinetics ----------------------------- */
    
        Species omega_;      ///< Net production rates [kg/(m^3·s) or mol/(m^3·s)]. @todo Confirm basis and units.
    
        /* ---------------- Miscellaneous ------------------------ */
    
        double N_;           ///< For miscellaneous use.
        double D_;           ///< For miscellaneous use.
    
        /* ---------------- Internal helpers --------------------- */
    
        /**
         * @brief Update thermodynamic properties for current temperature/composition.
         * @details
         *   @todo Document basis conversions for cp_→cp_bar_ and h_→h_bar_.
         * @post cp_, h_, cp_bar_, h_bar_ consistent with @ref T_ and @ref MW_.
         */
        void computeThermoProperties();
    
        /**
         * @brief Update chemical source terms for current state.
         * @details
         *   @todo Note any pressure-dependence through @ref C_ and libraries used by @ref source().
         * @post omega_ contains net production rates for each species.
         */
        void computeProductionRates();
    
        /**
         * @brief Compute numerator N_ for dT/dt.
         * @note  Currently unused.
         */
        void computedTtNumerator();
    
        /**
         * @brief Compute denominator D_ for dT/dt.
         * @details
         * @note  Currently unused.
         */
        void computedTtDenominator();
};

#endif /* SRC_REACTOR_IDEAL_GAS_CONST_PRESSURE_ADIABATIC_REACTOR */

