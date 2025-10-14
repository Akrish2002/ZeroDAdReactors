#include "CVODESSerialIntegrator.h"
#include "IdealGasConstPressureAdiabaticReactor.h"
#include "Utility.h"
#include "IdealGasConstPressureAdiabaticReactorAdapter.h"
#include "ChemConfig.h"

#include <iostream>


int main()
{

/*------------------------------Chemgen routine------------------------------*/

    /* 1. Load the .yaml file
     * 2. Run to generate the .h files (But then how will it compile?) */

/*---------------------------------------------------------------------------*/


/*------------------------------Reactor--------------------------------------*/

    IdealGasConstPressureAdiabaticReactor reactor(10, 300.9470);

    /* Print stats */
    std::cout<<"--Temperature: "<<reactor.getTemperature()<<std::endl;
    std::cout<<"--Pressure: "<<reactor.getPressure()<<std::endl;
    std::cout<<"--Number of Species: "<<reactor.getNumberofSpecies()<<std::endl;
    reactor.getMW();
    reactor.getcp();
    reactor.getomega();

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/

    IdealGasConstPressureAdiabaticReactorAdapter adapter(reactor);
    CVODESSerialIntegrator integ(adapter);
    std::cout<<"--Number of Eqns: "<<integ.getNEQ()<<std::endl;
    integ.initializeandsetupsolver();
    integ.integrate();

/*---------------------------------------------------------------------------*/
    
    /* Print stats */
    //std::cout<<
    

}

