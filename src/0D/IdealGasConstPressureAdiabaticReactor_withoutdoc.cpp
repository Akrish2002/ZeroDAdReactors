#include "IdealGasConstPressureAdiabaticReactor.h"
#include "ChemConfig.h"


#include "multiply_divide.h"
#include "pow_gen.h"
#include "exp_gen.h"
#include "array_handling.h"
#include "constants.h"
#include "thermally_perfect.h"
#include "arrhenius.h"
#include "third_body.h"
#include "falloff_troe.h"
#include "falloff_lindemann.h"
#include "falloff_sri.h"
#include "pressure_dependent_arrhenius.h"
#include "reactions.h"
#include "source.h"
#include "chemical_state_functions.h"


/* Public Functions */
IdealGasConstPressureAdiabaticReactor::IdealGasConstPressureAdiabaticReactor(int n, double temperature, double pressure)
{
    n_species_ = n;
    T_         = temperature;
    P_         = pressure;
    MW_        = molecular_weights();
    N_         = 0;
    D_         = 0;
    double temp= 0;
    computeThermoProperties();

    for(int i = 0; i < n_species_ ; i++) 
    {
        Y_[i] = 0.1;                                                            /* This is hardcoded, have to change as input */
    }

    for(int i = 0; i < n_species_ ; i++) 
    {
        temp += Y_[i] / MW_[i];
    }
    MWtot_ = 1 / temp;

    for(int i = 0; i < n_species_ ; i++) 
    {
        	//cp_[i]       = 0;
        	//cp_bar_[i]   = 0;

        	//h_[i]        = 0;
        	//h_bar_[i]    = 0;

            //Y_[i]        = 0.0;
            C_[i]        = (P_ * MWtot_ * Y_[i]) / (ChemConfig::Ru * T_ * MW_[i]); /* Concentration */

        	//omega_[i]    = 0;

        	/* MW_[i]       = 0; */
    }

    computeProductionRates();
}


void IdealGasConstPressureAdiabaticReactor::getProperties()
{
    computeThermoProperties();
    computeProductionRates();
}


/* Unused */
void IdealGasConstPressureAdiabaticReactor::getdTtNumerator()
{
    computedTtNumerator();
}


/* Unused */
void IdealGasConstPressureAdiabaticReactor::getdTtDenominator()
{
    computedTtDenominator();
}


int IdealGasConstPressureAdiabaticReactor::setNEQ()
{
    return(n_species_ + 1);
}


void IdealGasConstPressureAdiabaticReactor::setInitialState(double* y)
{
    /* For setting the temperature */
    //y[0] = ChemConfig::TEMP0;
    y[0] = T_;
    
    /* For mole fraction of species */
    for(int i = 1; i <= n_species_ ; i++)
    {
        /* y[i] = Y0_[i - 1]; */
        y[i] = Y_[i - 1];
    }

}


void IdealGasConstPressureAdiabaticReactor::setState(double* y, double temperature)
{
    /* Y_[0] = temperature; */
    T_ = temperature;

    for(int i = 1; i <= n_species_; i++)
    {
        /* y[i] = Y_[i - 1]; */
        Y_[i - 1] = y[i];

        C_[i - 1] = (P_ * MWtot_ * Y_[i - 1]) / (ChemConfig::Ru * T_ * MW_[i-1]);     /* Concentration */

    }

}


void IdealGasConstPressureAdiabaticReactor::evalRHS(double t, double* y, double* ydot)
{
    /* Declaring these vars for easy access */
    double temp       = y[0];
    double *massFracs = &y[1];
    double *dTdt      = &ydot[0];
    double *dYdt      = &ydot[1];
            
    /* Update the state */
    /* setState(y, y[0]); */
    setState(y, temp); 
            
    /* Getting the required properties */   
    getProperties();

    /* Energy eqn */
    N_ = 0.0;
    D_ = 0.0;
    for(int i = 0; i < n_species_; i++)
    {
        N_ += -h_bar_[i] * omega_[i];
        D_ += ((Y_[i] * P_ * MWtot_) * cp_bar_[i]) / (ChemConfig::Ru * T_ * MW_[i]);
    }
    *dTdt = N_ / D_;
    
    /* Species eqn */
    double omega_sum            = 0.0;
    double omega_temp           = 0.0;
    double concentration_sum    = 0.0;
    double concentration        = 0.0;

    for(int i = 0; i < n_species_; i++)
    {
        omega_sum         += omega_[i];
        concentration_sum += (Y_[i] * P_ * MWtot_) / (ChemConfig::Ru * T_ * MW_[i]);
    }
    
    for(int i = 0; i < n_species_; i++)
    {
        /* Change for loops from 1 to <= n_species_ */
        /* omega_temp   = omega_[i - 1] * (MW_[i - 1] * ChemConfig::Ru * T_) / (P_ * MWtot_); */
        /* ydot[i]      = omega_temp - Y_[i - 1] * ((omega_sum / concentration_sum) + (*dTdt / T_)); */
        omega_temp   = omega_[i] * (MW_[i] * ChemConfig::Ru * T_) / (P_ * MWtot_); 
        dYdt[i]      = omega_temp - Y_[i] * ((omega_sum / concentration_sum) + (*dTdt / T_));
    }
    
 }


double IdealGasConstPressureAdiabaticReactor::getTemperature()
{
    return T_;
}


double IdealGasConstPressureAdiabaticReactor::getNumberofSpecies()
{
    return n_species_;
}


double IdealGasConstPressureAdiabaticReactor::getPressure()
{
    return P_;
}


void IdealGasConstPressureAdiabaticReactor::getMW()
{
    std::cout<<"--Molecular weights:"<<std::endl;
    for(int i = 0; i < n_species_; i++)
    {
        std::cout<<"> "<<MW_[i]<<std::endl;
    }
    
    std::cout<<"--Molecular weight of mixture:"<<std::endl<<"> "<<MWtot_<<std::endl;
}


void IdealGasConstPressureAdiabaticReactor::getcp()
{
    std::cout<<"--Molar specific heats (J/kmol-K):"<<std::endl;
    for(int i = 0; i < n_species_; i++)
    {
        std::cout<<"> "<<cp_bar_[i]<<std::endl;
    }
}


void IdealGasConstPressureAdiabaticReactor::getomega()
{
    std::cout<<"--Production rates:"<<std::endl;
    for(int i = 0; i < n_species_; i++)
    {
        std::cout<<"> "<<omega_[i]<<std::endl;
    }
}


/* Private Functions */
void IdealGasConstPressureAdiabaticReactor::computeThermoProperties()
{
     cp_  = species_specific_heat_constant_pressure_mass_specific(T_);
     h_   = species_enthalpy_mass_specific(T_);
     for(int i = 0; i < n_species_ ; i++)
     {
         cp_bar_[i]  = cp_[i] * MW_[i];
         h_bar_[i]   = h_[i]  * MW_[i];
     }
}


void IdealGasConstPressureAdiabaticReactor::computeProductionRates()
{
    omega_ = source(C_, T_);
}


void IdealGasConstPressureAdiabaticReactor::computedTtNumerator()
{
	/* for(int i = 0; i < n_species_; i++)  
    {
         N_ += -h_bar_[i] * omega_[i]; 
    } 
    */
}


void IdealGasConstPressureAdiabaticReactor::computedTtDenominator()
{
	/* for(int i = 0; i < n_species_; i++)  
    {
         D_ += C_[i] * (cp_bar_[i]); 
    }
    */
} 


