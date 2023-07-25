#include "SolverSundials.hxx"

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* definition of type realtype          */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to spgmr SUNLinearSolver      */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver       */

#include "MultiDomain.hxx"
#include "SingleDomain.hxx"
#include "NodePool.hxx"
#include "LinearSolver.hxx"
#include "Exceptions.hxx"
#include "VectorOperators.hxx"
#include "SimulationProperties.hxx"

namespace Sundials
{
    struct UserData
    {
        explicit UserData(HygroThermFEM::SingleDomain & domain) : domain(domain)
        {}
        HygroThermFEM::SingleDomain & domain;
        size_t timestepIndex{0};
        std::vector<double> solution;
    };

    struct SunInitialization
    {
        SunInitialization(int error,
                          SUNContext ctx,
                          void * mem,
                          N_Vector uu,
                          N_Vector ud,
                          N_Vector rr,
                          N_Vector pp,
                          const std::vector<double> & uDot0,
                          const std::shared_ptr<UserData> & data) :
            error(error),
            ctx(ctx),
            mem(mem),
            uu(uu),
            ud(ud),
            rr(rr),
            pp(pp),
            uDot0(uDot0),
            data(data)
        {}

        ~SunInitialization()
        {
            IDAFree(&mem);
            N_VDestroy(uu);
            N_VDestroy(ud);
            N_VDestroy(rr);
            N_VDestroy(pp);
            SUNContext_Free(&ctx);
        }

        // Copy constructor
        SunInitialization(const SunInitialization & other) :
            error(other.error),
            ctx(other.ctx),           // assuming SUNContext can be copied directly
            uu(N_VClone(other.uu)),   // assuming this is the correct way to copy a N_Vector
            ud(N_VClone(other.ud)),
            rr(N_VClone(other.rr)),
            pp(N_VClone(other.pp)),
            uDot0(other.uDot0),
            data(new UserData(*other.data))   // deep copy UserData
        {
            mem = IDACreate(ctx);   // create new mem with the copied context
        }

        // Copy assignment operator
        SunInitialization & operator=(const SunInitialization & other)
        {
            if(this != &other)
            {   // self-assignment check expected
                // clean up the destination object
                IDAFree(&mem);
                N_VDestroy(uu);
                N_VDestroy(ud);
                N_VDestroy(rr);
                N_VDestroy(pp);
                SUNContext_Free(&ctx);

                // copy from source to destination
                error = other.error;
                ctx = other.ctx;           // assuming SUNContext can be copied directly
                uu = N_VClone(other.uu);   // assuming this is the correct way to copy a N_Vector
                ud = N_VClone(other.ud);
                rr = N_VClone(other.rr);
                pp = N_VClone(other.pp);
                uDot0 = other.uDot0;
                data = other.data;

                mem = IDACreate(ctx);   // create new mem with the copied context
            }
            return *this;
        }

        int error{0};
        SUNContext ctx{nullptr};
        void * mem{nullptr};
        N_Vector uu{nullptr};
        N_Vector ud{nullptr};
        N_Vector rr{nullptr};
        N_Vector pp{nullptr};
        std::vector<double> uDot0;
        std::shared_ptr<UserData> data;
    };

    int residual(realtype, N_Vector yy, N_Vector yp, N_Vector rr, void * user_data)
    {
        /* Initialize rr to uu, to take care of boundary equations.
         * ... this should only matter for dirichlet conditions*/
        N_VScale(1.0, yy, rr);

        const realtype * yval = N_VGetArrayPointer(yy);
        const realtype * ypval = N_VGetArrayPointer(yp);
        realtype * rval = N_VGetArrayPointer(rr);

        sunindextype neq = N_VGetLength(yy);
        UserData * data;
        data = (UserData *)user_data;

        auto timestepIndex = data->timestepIndex;

        // Get C matrix
        // Simon divides by dTime when getting mass matrix...
        // so we defined a new function to get raw (lumped) mass matrix
        auto C_eig = data->domain.elements.getCMatrix();
        // Get stiffness matrix
        auto K_eig = data->domain.elements.conductanceMatrix();
        // apply bcs no stiffness matrix
        K_eig += data->domain.boundaryConditions.HMatrix(timestepIndex);
        // RHS vector
        auto RHS = data->domain.boundaryConditions.RVector(timestepIndex);
        auto RHSe = data->domain.elements.RVector();

        // Sum up RHS and RHSe using trasnform
        std::transform(RHS.begin(), RHS.end(), RHSe.begin(), RHS.begin(), std::plus<>());

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

        // convert solution vector to format HygroThermFEM can understand and update nodal
        // solutions
        data->solution.clear();
        for(int i = 0; i < neq; i++)
        {
            data->solution.push_back(yval[i]);
        }

        // This does not seem necessary but keeping it here for now to check
        // HygroThermFEM::updateNodeValues(*data->solution, baseVariableOf(*data->domain));

        return 0;
    }

    std::vector<double> getInitialUdot(N_Vector uu, const UserData & user_data)
    {
        sunindextype neq = N_VGetLength(uu);

        auto timestepIndex = user_data.timestepIndex;

        auto C_eig = user_data.domain.elements.getCMatrix();
        auto K_eig = user_data.domain.elements.conductanceMatrix();
        K_eig += user_data.domain.boundaryConditions.HMatrix(timestepIndex);
        auto RHS = user_data.domain.boundaryConditions.RVector(timestepIndex);

        // turn into vector that Eigen can understand
        const realtype * uval = N_VGetArrayPointer(uu);

        std::vector<double> u0;
        u0.reserve(static_cast<size_t>(neq));
        for(int i = 0; i < neq; i++)
        {
            u0.push_back(uval[i]);
        }

        auto fstar = (-1.0) * K_eig * u0;
        std::transform(fstar.begin(), fstar.end(), RHS.begin(), fstar.begin(), std::plus<>());

        return HygroThermFEM::CLinearSolver::solveEigen(C_eig, fstar);
    }

    N_Vector CreateVector(int neq, SUNContext ctx)
    {
        return N_VNew_Serial(neq, ctx);
    }

    N_Vector CreateVectorFromData(const std::vector<double> & data, SUNContext ctx)
    {
        return N_VMake_Serial(data.size(), const_cast<realtype *>(data.data()), ctx);
    }

    void InitUd(const std::vector<double> & uDot0, N_Vector ud)
    {
        realtype * udvals = N_VGetArrayPointer(ud);
        for(size_t j = 0u; j < uDot0.size(); j++)
        {
            udvals[j] = uDot0[j];
        }
    }

    std::unique_ptr<UserData> CreateUserData(HygroThermFEM::SingleDomain & domain)
    {
        auto data{std::make_unique<UserData>(domain)};
        data->timestepIndex = 0;
        return data;
    }

    SUNLinearSolver CreateSolver(N_Vector uu, SUNMatrix A, SUNContext ctx)
    {
        return SUNLinSol_Band(uu, A, ctx);
    }

    std::map<HygroThermFEM::DomainType, double> errorForSundials{
      {HygroThermFEM::DomainType::Thermal,
       HygroThermFEM::SimulationProperties::Instance().errorToleranceTemperature()},
      {HygroThermFEM::DomainType::Moisture,
       HygroThermFEM::SimulationProperties::Instance().errorToleranceMoisture()}};

    SunInitialization initializeSolver(const std::vector<double> & initialValues,
                                       HygroThermFEM::SingleDomain & domain)
    {
        SUNContext ctx;
        int retval{SUNContext_Create(nullptr, &ctx)};
        const auto neq = HygroThermFEM::maxNodeIndex();

        N_Vector ud = CreateVector(neq, ctx);
        N_Vector rr = CreateVector(neq, ctx);
        N_Vector vatol = CreateVector(neq, ctx);
        N_Vector uu = CreateVectorFromData(initialValues, ctx);

        void * mem = IDACreate(ctx);
        auto data = CreateUserData(domain);
        retval = IDASetUserData(mem, data.get());

        auto uDot0 = getInitialUdot(uu, *data);
        InitUd(uDot0, ud);

        const realtype t0 = 0.0;
        retval = IDAInit(mem, residual, t0, uu, ud);

        realtype reltol{RCONST(errorForSundials.at(domain.domainType))};
        realtype abstol{RCONST(errorForSundials.at(domain.domainType) * 1000)};
        retval = IDASStolerances(mem, reltol, abstol);
        N_VConst(abstol, vatol);

        SUNMatrix A = SUNBandMatrix(neq, neq, neq, ctx);
        SUNLinearSolver LS = CreateSolver(uu, A, ctx);

        retval = IDASetLinearSolver(mem, LS, A);
        const auto maxSteps{HygroThermFEM::SimulationProperties::Instance().maxNumberOfIterations()
                            * 100};
        retval = IDASetMaxNumSteps(mem, maxSteps);

        N_Vector pp = N_VClone(uu);

        return {retval, ctx, mem, uu, ud, rr, pp, uDot0, std::move(data)};
    }

    void resetValues(const std::vector<double> & values, SunInitialization & sunInit)
    {
        realtype * yval = N_VGetArrayPointer(sunInit.uu);
        for(size_t i = 0u; i < values.size(); ++i)
        {
            yval[i] = values[i];
        }
    }

    struct SolverPimpl
    {
        SolverPimpl(HygroThermFEM::SingleDomain & domain,
                    const std::vector<double> & previousTimestepValue) :
            sunInit(Sundials::initializeSolver(previousTimestepValue, domain)),
            previousTimestepValue(previousTimestepValue)
        {}

        int solve(double t_DTime)
        {
            resetValues(previousTimestepValue, sunInit);
            return IDASolve(sunInit.mem, t_DTime, &currentTime, sunInit.uu, sunInit.ud, IDA_NORMAL);
        }

        SunInitialization sunInit;
        std::vector<double> previousTimestepValue;
        double currentTime{0.0};
    };

    TransientSingleDomainSundials::TransientSingleDomainSundials(
      HygroThermFEM::SingleDomain & domain) :
        SingleDomainTransientSolver(domain)
    {}

    HygroThermFEM::SingleTimestepSolution TransientSingleDomainSundials::transient(
      const std::vector<double> & previousTimestepValues, double t_DTime, size_t timestepIndex)
    {
        auto & pimpl = getSolverPimpl(m_Domain, previousTimestepValues);

        auto retval{pimpl->solve(t_DTime * (timestepIndex + 1))};
        if(retval != IDA_SUCCESS)
        {
            throw HygroThermFEM::SolutionFailedToConvergeException();
        }
        HygroThermFEM::updateNodeValues(pimpl->sunInit.data->solution,
                                        HygroThermFEM::baseVariableOf(pimpl->sunInit.data->domain),
                                        false);
        return {pimpl->sunInit.data->solution, pimpl->currentTime};
    }

    HygroThermFEM::SingleTimestepSolution
      TransientSingleDomainSundials::transient(const std::vector<double> & previousTimestepValues,
                                               double t_DTime)
    {
        return transient(previousTimestepValues, t_DTime, 0);
    }

    std::shared_ptr<SolverPimpl> & TransientSingleDomainSundials::getSolverPimpl(
      HygroThermFEM::SingleDomain & domain, const std::vector<double> & previousTimestepValues)
    {
        if(!m_pimpl)
        {
            m_pimpl = std::make_shared<SolverPimpl>(domain, previousTimestepValues);
            m_pimpl->currentTime = 0.0;
        }

        return m_pimpl;
    }

    SolverIDA::SolverIDA(HygroThermFEM::MultiDomain & domain) : TransientSolver(domain)
    {}

    std::unique_ptr<HygroThermFEM::SingleDomainTransientSolver>
      SolverIDA::createSolver(HygroThermFEM::SingleDomain & domain)
    {
        return std::make_unique<Sundials::TransientSingleDomainSundials>(domain);
    }

    std::vector<HygroThermFEM::SingleTimestepSolution>
      transient(HygroThermFEM::SingleDomain & domain,
                const std::vector<double> & previousTimestepValues,
                double t_DTime,
                size_t nTimesteps)
    {
        std::vector<HygroThermFEM::SingleTimestepSolution> solution;
        auto sunInit{initializeSolver(previousTimestepValues, domain)};
        auto currentTime{0.0};
        for(auto i = 0u; i < nTimesteps; ++i)
        {
            auto retval = IDASolve(
              sunInit.mem, t_DTime * (i + 1), &currentTime, sunInit.uu, sunInit.ud, IDA_NORMAL);
            if(retval != IDA_SUCCESS)
            {
                throw HygroThermFEM::SolutionFailedToConvergeException();
            }
            solution.emplace_back(sunInit.data->solution, currentTime);
            HygroThermFEM::updateNodeValues(
              sunInit.data->solution, HygroThermFEM::baseVariableOf(sunInit.data->domain), true);
        }

        return solution;
    }
}   // namespace Sundials
