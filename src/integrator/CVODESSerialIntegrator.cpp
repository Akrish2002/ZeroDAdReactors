#include "CVODESSerialIntegrator.h"


static int cvode_rhs(double t, N_Vector y, N_Vector ydot, void* user_data)
{
    Utility *f = static_cast<Utility*>(user_data);

    double *ydata    = N_VGetArrayPointer(y);
    double *ydotdata = N_VGetArrayPointer(ydot);

    f->evalRHS(t, ydata, ydotdata);

    return 0;
}


static int check_retval(void* returnvalue, const char* funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL)
  {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  /* Check if retval < 0 */
  else if (opt == 1)
  {
    retval = (int*)returnvalue;
    if (*retval < 0)
    {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
              funcname, *retval);
      return (1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL)
  {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return (1);
  }

  return (0);
}


static void PrintOutput(sunrealtype t, sunrealtype y1, sunrealtype y2,
                        sunrealtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}


static void PrintAllOutput(sunrealtype t, N_Vector y, int NEQ)
{
    printf("At t = %0.4e"); 
    double *ydata = N_VGetArrayPointer(y);

    for(int i = 0; i < NEQ; i++)
    {
        std::cout<<std::endl;
        std::cout<<">"<<ydata[i];
    }
    std::cout<<std::endl;
  return;
}


/* Made a subtle mistake here, I thought using model_ was valid, but it gave me garbage though */
CVODESSerialIntegrator::CVODESSerialIntegrator(Utility &model, int debug) : model_(model), NEQ_(model.setNEQ())            
{
    cvode_mem_  = nullptr;                        
    abstol_     = nullptr;
    RTOL_       = 1.0e-8;
    ATOL_.resize(NEQ_);
    for(int i = 0; i < NEQ_; i++)
    {
        ATOL_[i] = 1e-8;
    }
    y_          = nullptr;
    A_          = nullptr;
    LS_         = nullptr;

    time0_      = 0.0;
    timeop_     = 1.0;                                                          /* Hardcoded */
    TMULT_      = 10.0;                                                         /* Hardcoded */
    steps_      = 12;

    debug_      = debug;

    if(debug_ == 1)
    {
        std::cout<<"--Constructor of CVODES implemented!"<<std::endl;
    }
} 


int CVODESSerialIntegrator::allocateMemory()
{
    int retval;    

    if(debug_ == 1)
    {
        std::cout<<"--Starting memory allocation!"<<std::endl;
    }

    retval = SUNContext_Create(SUN_COMM_NULL, &sunctx_);

    /* 1. Allocating memory for y_ N_Vector */
    y_ = N_VNew_Serial(NEQ_, sunctx_);
    if(check_retval((void*)y_, "N_VNew_Serial", 0))
    {
        return (1);
    }

    /* 2. Allocating memory for abstol N_Vector */
    abstol_ = N_VNew_Serial(NEQ_, sunctx_);
    if(check_retval((void*)abstol_, "N_VNew_Serial", 0))
    {
        return (1);
    }
    
    if(debug_ == 1)
    {
        std::cout<<"--Memory allocated!"<<std::endl;
    } 

    return(0);
}


void CVODESSerialIntegrator::setInitialState()
{
    double *ydata     = N_VGetArrayPointer(y_);
    model_.setInitialState(ydata);

    if(debug_ == 1)
    {
        PrintAllOutput(0.0, y_, NEQ_);
    }
}


void CVODESSerialIntegrator::setTolerances()
{
    for(int i = 0; i < NEQ_; i++)
    {
        NV_Ith_S(abstol_, i) = ATOL_[i]; 
    }

    if(debug_ == 1)
    {
        std::cout<<"--Tolerances set!"<<std::endl;
    }

}


void CVODESSerialIntegrator::allocatesolverMemoryandMethod(int method)
{
    cvode_mem_ = CVodeCreate(method, sunctx_);
    
    if(debug_ == 1)
    {
        std::cout<<"--Allocating solver memory and method"<<std::endl;
    }
}


void CVODESSerialIntegrator::initializeintegratorMemoryandRHS()
{
    int flag = CVodeInit(cvode_mem_, cvode_rhs, time0_, y_);                    /* Initializing integrator memory and specify the right hand side fn, intial time and initial dependend variable vector */
    int flag1 = CVodeSetUserData(cvode_mem_, &model_);                          /* There is no error checker with flag1 */

    if(flag != CV_SUCCESS)
    {
        /* throw */
        if(flag == CV_MEM_FAIL)
        {
            std::cout<<"Memory allocation failed!"<<std::endl;
        }

        else if(flag == CV_ILL_INPUT)
        {
            std::cout<<"Illegal value for CVodeInit input argument"<<std::endl;
        }

        else
        {
            std::cout<<"CVodeInit Failed"<<std::endl;
        }
    }

    if(debug_ == 1)
    {
        std::cout<<"--Initialized integrator memory and RHS"<<std::endl;
    }

}

void CVODESSerialIntegrator::setrelativeTolerance()
{
    int flag = CVodeSVtolerances(cvode_mem_, RTOL_, abstol_);
    if (check_retval(&flag, "CVodeSVtolerances", 1)) 
    { 
        std::cout<<"Success in setting relative tolerance!"<<std::endl;
    }
}

void CVODESSerialIntegrator::createSUNDenseMatrix()
{
    A_ = SUNDenseMatrix(NEQ_, NEQ_, sunctx_);
    /* Write error checks */
}


void CVODESSerialIntegrator::createSUNLinSolObject()
{
    LS_ = SUNLinSol_Dense(y_, A_, sunctx_);
    /* Write error checks */
}


void CVODESSerialIntegrator::attachMatrixandLinSol()
{
    int flag = CVodeSetLinearSolver(cvode_mem_, LS_, A_);
    /* Write error checks using the flag value */
}
    

void CVODESSerialIntegrator::openfileforprinting()
{
    FID_ = fopen("0DCPAdReactor_Stats.csv", "w");
}


void CVODESSerialIntegrator::destroyN_Vectors()
{
    N_VDestroy(y_);       
    N_VDestroy(abstol_);       
}


void CVODESSerialIntegrator::freeblockmemory()
{
    CVodeFree(&cvode_mem_);
}


void CVODESSerialIntegrator::freesolvermemory()
{
    SUNLinSolFree(LS_);
}


void CVODESSerialIntegrator::freematrix()
{
    SUNMatDestroy(A_);
}


void CVODESSerialIntegrator::freeSUNDIALScontext()
{
    SUNContext_Free(&sunctx_);
}


void CVODESSerialIntegrator::closefile()
{
    fclose(FID_);
}


/* Public functions */
void CVODESSerialIntegrator::initializeandsetupsolver()
{
    allocateMemory();
    setInitialState();
    setTolerances();
    allocatesolverMemoryandMethod(1);
    initializeintegratorMemoryandRHS();
    setrelativeTolerance();
    createSUNDenseMatrix();
    createSUNLinSolObject(); 
    attachMatrixandLinSol();
    openfileforprinting();

    if(debug_ == 1)
    {
        std::cout<<"--Initialized and setup solver!"<<std::endl;
    }

}


void CVODESSerialIntegrator::integrate()
{
    int iout;
    int tout; 
    int flag; 

    iout = 0;
    tout = timeop_;

    while(1)
    {
        flag = CVode(cvode_mem_, tout, y_, &time_, CV_NORMAL);
        /* flag = CVode(cvode_mem_, tout, y_, &time_, CV_ONE_STEP); */
        /* PrintOutput(time_, NV_Ith_S(y_, 0), NV_Ith_S(y_, 1), NV_Ith_S(y_, 2)); */
        //PrintAllOutput(time_, y_, NEQ_);

        if(check_retval(&flag, "CVode", 1)) 
        {
            std::cout<<"--Error in integration!"<<std::endl;
            break;
        }

        if(flag == CV_SUCCESS)
        {
            iout++;
            tout *= TMULT_;
        }

        flag = CVodePrintAllStats(cvode_mem_, FID_, SUN_OUTPUTFORMAT_CSV);
        
        if(iout == steps_)
        {
            std::cout<<"Integration done!"<<std::endl;
            break;
        }
    }
    fclose(FID_);

    if(debug_ == 1)
    {
        printf("\nFinal Statistics:\n");
        flag = CVodePrintAllStats(cvode_mem_, stdout, SUN_OUTPUTFORMAT_TABLE);
    }

}


void CVODESSerialIntegrator::freeMemory()
{
    destroyN_Vectors();   
    freeblockmemory();    
    freesolvermemory();   
    freematrix();         
    freeSUNDIALScontext();
    closefile();
}


/* Debugger, getter fns */
int CVODESSerialIntegrator::getNEQ()
{
    return NEQ_;
}


double CVODESSerialIntegrator::getzeroEqn()
{
    return NV_Ith_S(y_, 0);
}


double CVODESSerialIntegrator::getfirstEqn()
{
    return NV_Ith_S(y_, 1);
}
