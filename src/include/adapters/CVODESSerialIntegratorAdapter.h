#ifndef 
#define


#include "Utility.h"
#include "CVODESSerialIntegrator.h"


class CVODESSerialIntegratorAdapter : public Utility
{

    public:
        explicit CVODESSerialIntegratorAdapter
        (CVODESSerialIntegrator &var) : var_(var) {}
        
    private:

        CVODESSerialIntegrator& var_; 

}
