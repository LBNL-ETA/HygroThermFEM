#include <iostream>
#include <string>
#include <memory>

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
//#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <sunmatrix/sunmatrix_sparse.h>    /* access to sparse SUNMatrix           */
#include <sundials/sundials_types.h>   /* definition of type realtype          */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to spgmr SUNLinearSolver      */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver       */

#include "HygroThermFEM2D.hxx"

using HygroThermFEM::NodePool;
using HygroThermFEM::MaterialPool;
using HygroThermFEM::State;

typedef struct {
  HygroThermFEM::MoistureDomain* domain;
  //realtype dTime;
  size_t timestepIndex{0};
  N_Vector pp;
} *UserData;


/* 
 * This is the central piece for interfacing HygroThermFEM and IDA from SUNDIALS
 * It computes the residual for F(u_dot,u,t) = 0 for a given u_dot and u
 * For us it looks like F = C(u)*u_dot + K(u)*u -RHS = 0
 *    where C is the (possibly) nonlinear capacitance matrix 
 *    and K is the (possibly) nonlinear stiffness matrix
 */
int residual(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{

    realtype *yval, *ypval, *rval;

    /* Initialize rr to uu, to take care of boundary equations. 
     * ... this should only matter for dirichlet conditions*/
    N_VScale(1.0, yy, rr);

    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);

    sunindextype neq = N_VGetLength(yy);
    UserData data;
    data = (UserData) user_data;

    auto timestepIndex = data->timestepIndex;

    // Get C matrix 
    // Simon divides by dTime when getting mass matrix... 
    // so we defined a new function to get raw (lumped) mass matrix
    auto C_eig = data->domain->m_Elements.getCMatrix();
    // Get stiffness matrix
    auto K_eig = data->domain->m_Elements.conductanceMatrix();
    // apply bcs no stiffness matrix
    K_eig += data->domain->m_BCs.HMatrix(timestepIndex);
    // RHS vector
    auto RHS = data->domain->m_BCs.RVector(timestepIndex);
    // apply bcs on RHS (funkyness to add vectors together... thanks, StackOverflow)
    auto RHSbc = data->domain->m_BCs.RVector(timestepIndex);
    std::transform (RHS.begin(), RHS.end(), RHSbc.begin(), RHS.begin(), std::plus<double>());

    double LHS;
    for (int i=0; i<neq; i++) {
        LHS = 0.0;
        for (int j=0; j<neq; j++) { 
            LHS += C_eig.m_Matrix.coeff(i,j) * ypval[j] + K_eig.m_Matrix.coeff(i,j) * yval[j]; 
        }
        rval[i] = LHS - RHS[i];
    }

    // convert solution vector to format HygroThermFEM can understand and update nodal solutions
    std::vector<double> yvec;
    yvec.clear();
    for (int i=0; i<neq; i++) {
      yvec.push_back(yval[i]);
    }
    updateNodeValues(yvec, false);
    NodePool::Instance().updateNodeValues(yvec, HygroThermFEM::BaseVariable::humidity, false);

  return(0);
}

/*
 * This fella finds an initial time derivative of the solution vector using eigen
 * In other words, solves 
 *  C(u)*u_dot + K(u)*u = RHS 
 *  --> C(u)*u_dot = fstar = RHS - K(u)*u
 * for u_dot given an initial u
 * */
std::vector<double> getInitialUdot(N_Vector uu, void *user_data)
{
  std::vector<double> udot0;
  std::vector<double> u0;
  realtype *uval;
  sunindextype neq = N_VGetLength(uu);
  UserData data;
  data = (UserData) user_data;
  auto timestepIndex = data->timestepIndex;

  /*Grab system information from HygroThermFEM*/
  auto C_eig = data->domain->m_Elements.getCMatrix();
  auto K_eig = data->domain->m_Elements.conductanceMatrix();
  K_eig += data->domain->m_BCs.HMatrix(timestepIndex);
  auto RHS = data->domain->m_BCs.RVector(timestepIndex);
  auto RHSbc = data->domain->m_BCs.RVector(timestepIndex);
  std::transform (RHS.begin(), RHS.end(), RHSbc.begin(), RHS.begin(), std::plus<double>());
  
  // turn into vector that Eigen can understand
  uval = N_VGetArrayPointer(uu);
  for (int i=0; i<neq; i++) {
    u0.push_back(uval[i]);
  }

  auto fstar = (-1.0)*K_eig*u0;
  std::transform (fstar.begin(), fstar.end(), RHS.begin(), fstar.begin(), std::plus<double>());

  udot0 = HygroThermFEM::CLinearSolver::solveEigen(C_eig, fstar);

  return(udot0);
}

int main(int, char**){
    std::cout << "This is a demo to show sundials can help with unstable cases." << std::endl;
    std::cout << "It uses Moisture_2D_TwoElements_1.TestExample_1.unit.cxx as the test case with high initial humidity where existing solver fails.\n" << std::endl;

    const auto initialTemperature = 20;
    const auto initialHumidity = 0.999;
    const auto initialPressure = 101325.0;
    const auto liquidPercent = 1.0;

    const State state(initialTemperature, initialHumidity, initialPressure, liquidPercent);
    NodePool::Instance().createNode(1, 0.15, 0.05, state);
    NodePool::Instance().createNode(2, 0.15, 0, state);
    NodePool::Instance().createNode(3, 0.05, 0.05, state);
    NodePool::Instance().createNode(4, 0.05, 0, state);
    NodePool::Instance().createNode(5, 0, 0.05, state);
    NodePool::Instance().createNode(6, 0, 0, state);

    // Material Properties (Cottaer Sandstone)
    const double thermalConductivityDry{1.8};
    const double density{2050.0};
    const double porosity{0.22};
    const double specificHeatCapacityDry{850.0};
    const double diffusionResistanceFactor{15.0};
    const std::vector<FenestrationCommon::point> thermalConductivityMoistureDependent = {
      {0.0, 1.8}, {180, 1.8}};
    const double thermalConductivityMeasuredAtTemperature{0};
    const std::vector<FenestrationCommon::point> thermalConductivityTemperatureDependent = {
      {0.0, 1.8}, {1, 1.8}};
    const double thermalConductivityMeasuredAtHumidity{0};
    const std::vector<FenestrationCommon::point> liquidTransportationCurve = {{0, 0},
                                                                              {27, 1E-8},
                                                                              {45, 1.1E-8},
                                                                              {90, 2E-8},
                                                                              {126, 3.5E-8},
                                                                              {144, 5E-8},
                                                                              {162, 1E-7},
                                                                              {171, 2E-7},
                                                                              {180, 7E-7}};

    const std::vector<FenestrationCommon::point> moistureStorageFunction = {{0, 0},
                                                                            {0.5, 5.3},
                                                                            {0.65, 8.4},
                                                                            {0.8, 12},
                                                                            {0.93, 17},
                                                                            {0.95, 25},
                                                                            {0.99, 63},
                                                                            {0.995, 83},
                                                                            {0.999, 120},
                                                                            {1, 180}};

    auto & material =
      MaterialPool::Instance().createSolidMaterial("Cottaer Sandstone",
                                                   thermalConductivityDry,
                                                   density,
                                                   porosity,
                                                   specificHeatCapacityDry,
                                                   diffusionResistanceFactor,
                                                   thermalConductivityMoistureDependent,
                                                   thermalConductivityMeasuredAtTemperature,
                                                   thermalConductivityTemperatureDependent,
                                                   thermalConductivityMeasuredAtHumidity,
                                                   liquidTransportationCurve,
                                                   moistureStorageFunction);

    HygroThermFEM::MoistureDomain domain;

    /// Create elements
    domain.createElement(3, 4, 2, 1, material.name());
    domain.createElement(6, 4, 3, 5, material.name());

    // Create Boundary Conditions
    const auto airTemperature = 20.0;
    const auto airHumidity = 0.0;
    const auto hc = 10.0;

    const HygroThermFEM::FixedBCHCCoefficients bcCoeff{airTemperature, hc, airHumidity};

    domain.createBC_FixedHc(5, 6, bcCoeff);

    const auto dTime = 3600;
    const auto nSteps = 30;

    auto humidities = NodePool::Instance().properties(HygroThermFEM::Variable::humidity);
    std::vector<double> timesteps;
    std::vector<std::vector<double>> waterContentSolution;
    std::vector<std::vector<HygroThermFEM::NodeFlux>> fluxSolution;
   /*
   * Replace the time integration and the nonlinear solve with IDA from SUNDIALS... 
   * HygroThermFEM is essentially unchanged.
   */

    size_t timestepIndex{0};

    int retval; // IDA functions return 0 (good) or something less than zero (Bad)

// Initalize Solver
    // SUNContext object is the orchestra conductor
    SUNContext ctx;
    retval = SUNContext_Create(NULL, &ctx);

    // make some SUNDIALS-native vectors
    N_Vector ud, uu, rr, vatol;
    const auto neq = NodePool::Instance().maxIndex();
    uu = N_VNew_Serial(neq, ctx);
    ud = N_VNew_Serial(neq, ctx);
    rr = N_VNew_Serial(neq, ctx);
    vatol = N_VNew_Serial(neq, ctx);

    // set initial condition on humidity
    N_VConst(initialHumidity,uu);

    N_VConst(29.0,rr);

    // get access to SUNDIALS arrays
    realtype *uuvals, *udvals, *rvals;
    uuvals = N_VGetArrayPointer(uu);
    udvals = N_VGetArrayPointer(ud);
    rvals = N_VGetArrayPointer(rr);

    // initialize solution with IDA
    void *mem = nullptr;
    mem = IDACreate(ctx);

    // tell sundials how to get to domain object and stuff so it can construct a residual
    UserData data;
    data = (UserData) malloc(sizeof *data);
    data->domain = &domain;
    //data->dTime = dTime;
    data->timestepIndex = timestepIndex;   
    // this pp vector is a carryover from making a user defined preconditioner... see **_messy
    data->pp = nullptr;
    data->pp = N_VClone(uu);
    retval = IDASetUserData(mem, data);  

    // This gets the consistent IC for udot
    std::vector<double> ud0;
    ud0 = getInitialUdot(uu, data);
    N_VConst(0.0,ud);
    for (int j = 0; j < neq; j++) {
      udvals[j] = ud0[j];
    }

    const realtype t0 = 0.0;
    retval = IDAInit(mem, residual, t0, uu, ud);
    
    realtype reltol = RCONST(1.0e-5);
    realtype abstol = RCONST(1.0e-4);
    realtype tret, tsd, t;
    t = RCONST(0.0);
    retval = IDASStolerances(mem, reltol, abstol);
    N_VConst(abstol,vatol);
    //retval = IDASVtolerances(mem, reltol, vatol);

    // Chose the solver (Still in initialize)
    SUNLinearSolver LS;
    
     /* Create banded SUNMatrix for use in linear solves */
    SUNMatrix A;
    A = SUNBandMatrix(neq, neq, neq, ctx);

    /* Create banded SUNLinearSolver object */
    LS = SUNLinSol_Band(uu, A, ctx);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);

    retval = IDASetMaxNumSteps(mem, 10000);
    //retval = IDASetMinStep(mem, dTime/1000.);

    std::vector<double> waterContent = {0.0,0.0,0.0,0.0,0.0,0.0};
    

    // This is transient solver
    for(unsigned i = 0; i < nSteps; ++i)
    {
        t = (i+1)*dTime;
        std::cout << "#################" << std::endl;
        std::cout << "# t_intended = " << t << std::endl;
        retval = IDASolve(mem, t, &tret, uu, ud, IDA_NORMAL);
        std::cout << "# t_IDA: " << tret << std::endl;
        std::cout << "#################" << std::endl;
        //std::cout << "IDASolve error code: " << retval << std::endl;
        waterContent = NodePool::Instance().properties(HygroThermFEM::Variable::water);
        waterContentSolution.push_back(waterContent);
        

    }



  std::vector<double> correctH20at86400 = {15.8460316, 15.8460316, 15.7388886, 15.7388886, 0.0420899289, 0.0420899289};
  std::vector<double> correctH20at3600 = {15.8461537, 15.8461537, 15.8448713, 15.8448713, 8.68767655, 8.68767655};

  for (int i = 0; i < humidities.size(); ++i)
  
    {
      std::cout << "correct: " << correctH20at86400[i] << ", with IDA: " << waterContent[i] << " (humidity: "  << uuvals[i] << ")" << std::endl;
    }

  // Leave no trace
  IDAFree(&mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(uu);
  N_VDestroy(ud);
  N_VDestroy(rr);
  N_VDestroy(data->pp);
  free(data);
  SUNContext_Free(&ctx);

  return (0);
}