#include "SolverSundials.hxx"

#include <ida/ida.h>                /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h> /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_sparse.h> /* access to sparse SUNMatrix           */
#include <sundials/sundials_types.h>    /* definition of type realtype          */
#include <sunlinsol/sunlinsol_spgmr.h>  /* access to spgmr SUNLinearSolver      */
#include <sunmatrix/sunmatrix_band.h>   /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>   /* access to band SUNLinearSolver       */

#include "SingleDomain.hxx"
#include "NodePool.hxx"
#include "LinearSolver.hxx"
#include "Exceptions.hxx"

namespace Sundials
{
    SunInitialization::~SunInitialization()
    {
        IDAFree(&mem);
        N_VDestroy(uu);
        N_VDestroy(ud);
        N_VDestroy(rr);
        N_VDestroy(data.data->pp);
        free(data.data);
        SUNContext_Free(&ctx);
    }

    int residual(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void * user_data)
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
        data = (UserData)user_data;

        auto timestepIndex = data->timestepIndex;

        // Get C matrix
        // Simon divides by dTime when getting mass matrix...
        // so we defined a new function to get raw (lumped) mass matrix
        auto C_eig = data->domain->elements.getCMatrix();
        // Get stiffness matrix
        auto K_eig = data->domain->elements.conductanceMatrix();
        // apply bcs no stiffness matrix
        K_eig += data->domain->boundaryConditions.HMatrix(timestepIndex);
        // RHS vector
        auto RHS = data->domain->boundaryConditions.RVector(timestepIndex);
        // apply bcs on RHS (funkyness to add vectors together... thanks, StackOverflow)
        auto RHSbc = data->domain->boundaryConditions.RVector(timestepIndex);
        std::transform(RHS.begin(), RHS.end(), RHSbc.begin(), RHS.begin(), std::plus<>());

        double LHS;
        for(int i = 0; i < neq; i++)
        {
            LHS = 0.0;
            for(int j = 0; j < neq; j++)
            {
                LHS += C_eig(i, j) * ypval[j] + K_eig(i, j) * yval[j];
            }
            rval[i] = LHS - RHS[i];
        }

        // convert solution vector to format HygroThermFEM can understand and update nodal solutions

        data->solution.clear();
        for(int i = 0; i < neq; i++)
        {
            data->solution.push_back(yval[i]);
        }

        return (0);
    }

    SunUserData getInitialUdot(N_Vector uu, void * user_data)
    {
        std::vector<double> udot0;
        std::vector<double> u0;
        realtype * uval;
        sunindextype neq = N_VGetLength(uu);
        UserData data;
        data = (UserData)user_data;
        auto timestepIndex = data->timestepIndex;

        /*Grab system information from HygroThermFEM*/
        auto C_eig = data->domain->elements.getCMatrix();
        auto K_eig = data->domain->elements.conductanceMatrix();
        K_eig += data->domain->boundaryConditions.HMatrix(timestepIndex);
        auto RHS = data->domain->boundaryConditions.RVector(timestepIndex);
        auto RHSbc = data->domain->boundaryConditions.RVector(timestepIndex);
        std::transform(RHS.begin(), RHS.end(), RHSbc.begin(), RHS.begin(), std::plus<>());

        // turn into vector that Eigen can understand
        uval = N_VGetArrayPointer(uu);
        for(int i = 0; i < neq; i++)
        {
            u0.push_back(uval[i]);
        }

        auto fstar = (-1.0) * K_eig * u0;
        std::transform(fstar.begin(), fstar.end(), RHS.begin(), fstar.begin(), std::plus<>());

        udot0 = HygroThermFEM::CLinearSolver::solveEigen(C_eig, fstar);

        return {udot0, data};
    }

    SunInitialization initializeSolver(double initialValue, HygroThermFEM::SingleDomain & domain)
    {
        int retval;   // IDA functions return 0 (good) or something less than zero (Bad)

        // Initalize Solver
        // SUNContext object is the orchestra conductor
        SUNContext ctx;
        retval = SUNContext_Create(nullptr, &ctx);

        // make some SUNDIALS-native vectors
        N_Vector ud, uu, rr, vatol;
        const auto neq = HygroThermFEM::maxNodeIndex();
        uu = N_VNew_Serial(neq, ctx);
        ud = N_VNew_Serial(neq, ctx);
        rr = N_VNew_Serial(neq, ctx);
        vatol = N_VNew_Serial(neq, ctx);

        // set initial condition
        N_VConst(initialValue, uu);

        N_VConst(29.0, rr);

        // get access to SUNDIALS arrays
        realtype *uuvals, *udvals, *rvals;
        uuvals = N_VGetArrayPointer(uu);
        udvals = N_VGetArrayPointer(ud);
        rvals = N_VGetArrayPointer(rr);

        // initialize solution with IDA
        void * mem = nullptr;
        mem = IDACreate(ctx);

        // tell sundials how to get to domain object and stuff so it can construct a residual
        UserData data;
        data = (UserData)malloc(sizeof *data);
        data->domain = &domain;
        // data->dTime = dTime;
        data->timestepIndex = 0;
        // this pp vector is a carryover from making a user defined preconditioner... see **_messy
        data->pp = nullptr;
        data->pp = N_VClone(uu);
        retval = IDASetUserData(mem, data);

        // This gets the consistent IC for udot
        SunUserData userData;

        userData = getInitialUdot(uu, data);
        N_VConst(0.0, ud);
        for(int j = 0; j < neq; j++)
        {
            udvals[j] = userData.uDot0[j];
        }

        const realtype t0 = 0.0;
        retval = IDAInit(mem, residual, t0, uu, ud);

        realtype reltol = RCONST(1.0e-5);
        realtype abstol = RCONST(1.0e-4);
        realtype tret, tsd, t;
        t = RCONST(0.0);
        retval = IDASStolerances(mem, reltol, abstol);
        N_VConst(abstol, vatol);
        // retval = IDASVtolerances(mem, reltol, vatol);

        // Chose the solver (Still in initialize)
        SUNLinearSolver LS;

        /* Create banded SUNMatrix for use in linear solves */
        SUNMatrix A;
        A = SUNBandMatrix(neq, neq, neq, ctx);

        /* Create banded SUNLinearSolver object */
        LS = SUNLinSol_Band(uu, A, ctx);

        /* Attach the matrix and linear solver */
        retval = IDASetLinearSolver(mem, LS, A);

        constexpr auto maxSteps = 10000;
        retval = IDASetMaxNumSteps(mem, maxSteps);
        // retval = IDASetMinStep(mem, dTime/1000.);

        return {retval, ctx, mem, uu, ud, rr, userData};
    }

    std::vector<HygroThermFEM::SingleTimestepSolution>
      transient(HygroThermFEM::SingleDomain & domain,
                const std::vector<double> & previousTimestepValues,
                double t_DTime)
    {
        std::vector<HygroThermFEM::SingleTimestepSolution> solution;
        auto sunInit{initializeSolver(0.99, domain)};
        auto currentTime{0.0};
        const auto nSteps{10u};
        for(auto i = 0u; i < nSteps; ++i)
        {
            auto retval =
              IDASolve(sunInit.mem, t_DTime * (i + 1), &currentTime, sunInit.uu, sunInit.ud, IDA_NORMAL);
                if(retval != IDA_SUCCESS)
                {
                    throw HygroThermFEM::SolutionFailedToConvergeException();;
                }
            solution.emplace_back(sunInit.data.data->solution, currentTime);
        }

        return solution;
    }
}   // namespace Sundials