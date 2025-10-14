#ifndef SRC_CHEM_CONFIG
#define SRC_CHEM_CONFIG


namespace ChemConfig
{
    /* Physical constants */
    inline constexpr double TEMP0 = 350.00;                                     /* Kelvin */
    inline constexpr double p0    = 1.0e5;                                      /* Pascal */
    inline constexpr double Ru    = 8.314;                                      /* J/mol-K */
    
    /* Tolerances */
    inline constexpr double RTOL = 1.0e-4;
    inline constexpr double ATOL = 1.0e-8;

}


#endif
