/**
 * @file IdealGasConstPressureAdiabaticReactor.h
 * @brief Module for implementing constant pressure reactor with ideal gas law
 */

#ifndef SRC_REACTOR_IDEAL_GAS_CONST_PRESSURE_ADIABATIC_REACTOR
#define SRC_REACTOR_IDEAL_GAS_CONST_PRESSURE_ADIABATIC_REACTOR


#include<iostream>
#include<vector>
#include<cmath>
#include<string>
#include<fstream>
#include<array>
#include<algorithm> 								
#include<cstring>								    


/* Chemgen header files */
#include "types_inl.h"                                                          /* For Species */


/**
 * @brief   Adiabatic reactor + constant pressure + ideal gas law
 * @tparam   
 */
class IdealGasConstPressureAdiabaticReactor 
{
 
    public:
        /* Only dealing with simple compressible substance */
        IdealGasConstPressureAdiabaticReactor(int n, double temperature = 300, double pressure = 101325.0);
        
        void getProperties();
        void getdTtNumerator();
        void getdTtDenominator();
        int  setNEQ();
        void setInitialState(double* y);
        void setState(double* y, double temperature);                           /* This is to update the reactor's variables with the new y vector, where y --> [T, Y1, Y2, ... YN] */
        void evalRHS(double t, double* y, double* ydot);
        
        /* Debugger, miscellaneous fns */
        double getTemperature();
        double getPressure();
        double getNumberofSpecies();
        void getMW();
        void getcp();
        void getomega();


    private:
	    double Ctot_;
        int n_species_;

        /* Required thermo variables */
        double  T_;
        double  P_;

        double MWtot_;
        Species MW_;                                                            /* Molecular weights        */
	    Species cp_, cp_bar_; 
	    Species h_, h_bar_;

	    Species X_;                                                             /* Mole fraction            */
	    Species Y0_;                                                            /* Initial Mass fraction    */
	    Species Y_;                                                             /* Mass fraction            */
	    Species C_;                                                             /* Concentration            */

        Species omega_;                                                         /* Production rates         */

        /* Miscellaneous variables */
        double N_, D_;

        /* Functions */
        void computeThermoProperties();
        void computeProductionRates();
        void computedTtNumerator();
        void computedTtDenominator();
};


#endif
