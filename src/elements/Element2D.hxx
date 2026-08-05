#pragma once

#include <array>
#include <memory>
#include <span>
#include <type_traits>
#include <utility>

#include "Functions.hxx"
#include "Material.hxx"
#include "Materials.hxx"
#include "Node2D.hxx"
#include "Nodes.hxx"
#include "Quadrilateral2D.hxx"
#include "SquareMatrix.hxx"

namespace HygroThermFEM
{
    //! Constant that holds number of nodes in certain elementsCreator
    const std::size_t numOfQuadrilateralNodes = 4;

    //! \brief Dense, inline 4x4 matrix for element-level integration.
    //!
    //! The integration matrices are read sixteen times per integration point on every call to
    //! integrate(), which happens once per coefficient term per element per assembly. Holding
    //! them in a sparse matrix meant every one of those reads was a binary search through an
    //! index structure, to reach one of sixteen values that are almost all non-zero anyway.
    struct ElementMatrix2D
    {
        double & operator()(const std::size_t row, const std::size_t col)
        {
            return values[row * numOfQuadrilateralNodes + col];
        }

        double operator()(const std::size_t row, const std::size_t col) const
        {
            return values[row * numOfQuadrilateralNodes + col];
        }

        ElementMatrix2D & operator+=(const ElementMatrix2D & other)
        {
            for(std::size_t idx = 0; idx < values.size(); ++idx)
            {
                values[idx] += other.values[idx];
            }
            return *this;
        }

        //! Matrix-vector product on the stack; no heap allocation.
        [[nodiscard]] std::array<double, numOfQuadrilateralNodes>
          operator*(std::span<const double> vec) const
        {
            std::array<double, numOfQuadrilateralNodes> product{};
            for(std::size_t row = 0; row < numOfQuadrilateralNodes; ++row)
            {
                for(std::size_t col = 0; col < numOfQuadrilateralNodes; ++col)
                {
                    product[row] += (*this)(row, col) * vec[col];
                }
            }
            return product;
        }

        std::array<double, numOfQuadrilateralNodes * numOfQuadrilateralNodes> values{};
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  IQLEMatrix2D
    //////////////////////////////////////////////////////////////////////////////

    //! \brief General class for quadratic linear element integrator.
    //!
    //! Hold basic calculations for matrix integration in global coordinate system. Integration is
    //! performed over different equations and therefore, different integration matrices will be
    //! used. To create new integrator type, user need to inherit from this class.
    class IQLEIntegrator2D
    {
    public:
        virtual ~IQLEIntegrator2D() = default;
        IQLEIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);

        //! Integrate matrix over all points of integration.
        //!
        //! Not virtual: no integrator overrides it. The derived types differ only in the
        //! integration matrix their constructors build.
        ElementMatrix2D
          integrate(std::span<const double>
                      t_Values   //!< Nodal values for which integration will be performed
                    ) const;

        //! \brief Integrate with the coefficient interpolated at each Gauss point via the
        //! shape functions (coeff_g = sum_l psi_l(g) * t_Values[l]) instead of the pairwise
        //! nodal average 0.5*(k_row + k_col) used by integrate().
        //!
        //! The pairwise average breaks discrete conservation when the coefficient varies
        //! inside an element: the column sums sum_row M(row,col) * 0.5*(k_row + k_col) do
        //! not vanish, which acts as a spurious interior moisture source wherever the
        //! coefficient has a gradient (e.g. c_sat(T) under a temperature gradient). With
        //! the Gauss-point-interpolated coefficient the column sums telescope to exactly
        //! zero for ANY spatially varying coefficient, because sum_row grad(psi_row) = 0
        //! pointwise. Used by the moisture element; thermal assembly keeps integrate().
        ElementMatrix2D integrateInterpolated(std::span<const double> t_Values) const;

    protected:

        const QuadrilateralLinearGlobal2D & m_Global2D;

        //! This matrix will hold different forms of shape functions operations. This will mainly
        //! depend on what integrator will be used to integrate matrices.
        std::vector<ElementMatrix2D> m_IntegrationMatrix;
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDDuIntegrator2D
    //////////////////////////////////////////////////////////////////////////////

    //! \brief Class to handle conductance matrix in global coordinate system
    //! Conductance equation is D/Dx(Du/Dx)+D/Dy(Du/Dy) where u is state variable
    //!
    //! Not much is going on in this class except that during construction, adequate integration
    //! matrix is created. This integrator should be used in creation of conductance matrix.
    class QLEDDuIntegrator2D : public IQLEIntegrator2D
    {
    public:
        virtual ~QLEDDuIntegrator2D() = default;

        QLEDDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDpDuIntegrator2D
    //////////////////////////////////////////////////////////////////////////////

    //! \brief Class to handle double derivative matrix in global coordinate system
    //! Equation is (Dp/Dx)(Du/Dx)+(Dp/Dy)(Du/Dy) where u is state variable and
    //! p is independent variable that can have different values in different nodes.
    //!
    //! In this case, creation of integration matrix is little bit more complicated since it needs
    //! to include additional independent variable.
    class QLEDpDuIntegrator2D : public IQLEIntegrator2D
    {
    public:
        virtual ~QLEDpDuIntegrator2D() = default;

        QLEDpDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);

        //! Update independent variables (corresponds to variable p in above equation) and creates
        //! new integration matrix.
        void setIndependentVariables(std::span<const double> t_Values);
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  QLECapacitance2D
    //////////////////////////////////////////////////////////////////////////////

    //! \brief Class to handle capacitance matrix in global coordinate system. Equation is
    //! derived from C * Du/Dt where u is state variable.
    //!
    //! Integrator is acutally calculating simple integration matrix Ksi,i * Ksi,j which is used
    //! to calculate capacitance (or mass) matrix. It is more important to put this matrix at
    //! correct place in matrix system of equations.
    class QLECapacitanceIntegrator2D : public IQLEIntegrator2D
    {
    public:
        virtual ~QLECapacitanceIntegrator2D() = default;

        QLECapacitanceIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);
    };

    //! Simple structure to hold (x, y) direction of heat flux in the node.
    struct NodeFlux
    {
        NodeFlux(double t_x, double t_y) : x(t_x), y(t_y)
        {}
        double x;
        double y;
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D
    //////////////////////////////////////////////////////////////////////////////

    //! \brief Class that handles base quadratic linear 2D element functionality.
    //!
    //! First important functionality of this class is that it handles integration of various
    //! matrix types. In child classes, user is suppose to provide concrete integration matrices
    //! that will be integrated here. Beside integrating matrices, this class provides functionality
    //! to create right hand side vectors. It also hold information if element will work in linear
    //! solution.
    class IElementLinear2D
    {
    public:
        virtual ~IElementLinear2D() = default;
        //! Construction of base element. Constructor will make all necessary shape functions in
        //! background (since they are dependent on element geometry). In constructor will
        //! calculate weighting factors that will be used to calculate some local node properties.
        //! This is important since one node can belong to several different elements and when
        //! doing averaging, program will have information about how much of weight one node
        //! this element will have.
        IElementLinear2D(
          Nodes & nodePool,                //!< Reference to NodePool for node lookup
          Materials & materialPool,           //!< Reference to MaterialPool for material lookup
          size_t index1,                      //!< Node 1 index
          size_t index2,                      //!< Node 2 index
          size_t index3,                      //!< Node 3 index
          size_t index4,                      //!< Node 4 index
          const std::string & materialName,   //!< SolidMaterial assigned to the element
          Variable
            variable,   //! Variable is used to determine for which property flux will be calculated
          bool isLinear = true   //!< States if element equations can be solved by linear
                                 //!< approach or non-linear iterations needs to be used.
        );

        //! Integrates all matrices that are part of K * (D/Dx(Du/Dx) + D/Dy(Du/Dy)) equation.
        ElementMatrix2D DDuMatrices() const;

        //! Integrates all matrices that are part of K * ((Dp/Dx)(Du/Dx) + (Dp/Dy)(Du/Dy)) equation.
        ElementMatrix2D DpDuMatrices() const;

        //! Integrates all matrices that are part of K * Du/Dt equation.
        ElementMatrix2D capacitanceMatrices() const;

        //! Integrates right hand-side vector.
        std::array<double, numOfQuadrilateralNodes> rightSideVector() const;

        //! \brief Sets a constant volumetric source for this element [W/m^3 for the
        //! thermal domain]. Contributes the consistent load q * integral(psi_i dA) via
        //! volumetricSourceVector(); zero (the default) adds nothing.
        void setVolumetricSource(double value);

        //! True when a nonzero volumetric source is set; lets the assembly skip the
        //! element entirely in the common no-source case.
        [[nodiscard]] bool hasVolumetricSource() const;

        //! \brief Consistent load vector of the element's volumetric source,
        //! q * integral(psi_i dA); all zeros when no source is set.
        [[nodiscard]] std::array<double, numOfQuadrilateralNodes> volumetricSourceVector() const;

        std::vector<NodeFlux> flux() const;

        //! Heat flux (sigma = -k * grad(state)) evaluated at each 2x2 Gauss point, in the
        //! integration-point order. Unlike flux(), the values are not extrapolated to the corner
        //! nodes; the recovery-based error estimator consumes the raw Gauss-point flux.
        [[nodiscard]] std::vector<NodeFlux> fluxAtGaussPoints() const;

        //! Global (x, y) coordinates of the element's 2x2 Gauss points, in the same order as
        //! fluxAtGaussPoints().
        [[nodiscard]] std::vector<std::array<double, 2>> gaussPointGlobalCoordinates() const;

        //! Representative scalar conductivity of the element: the mean of the nodal conductivity
        //! values used in the flux calculation. For the constant-conductivity solid materials that
        //! dominate steady-state thermal models this equals the material's dry conductivity. Its
        //! reciprocal is the inverse constitutive scalar that defines the estimator's energy norm.
        [[nodiscard]] double meanConductivity() const;

        INode2D & getNode(std::size_t index) const;

        //! Returns true of element have both nodes. This is used to test if certain boundary
        //! condition is part of the element.
        bool haveBothNodes(size_t index1,   //!< Node 1 index
                           size_t index2    //!< Node 2 index
                           ) const;

        //! Returns indexes of all four nodes in the element
        std::vector<std::size_t> nodeIndexes() const;

        //! Returns material that is assigned to the element.
        const IMaterial & getMaterial() const;

        //! Returns if element satisfies linear problem. Not virtual: neither element type
        //! overrides it, they differ only in the terms their constructors register.
        bool isLinear() const;

    protected:
        //! Template function that will add K*(D/Dx(Du/Dx) + D/Dy(Du/Dy)) matrix into the system,
        //! i.e. transport of the state variable itself. The special case of the two-argument
        //! DDu with a unit nodal factor, registered without one so the assembly skips the
        //! scaling.
        template<typename T>
        void DDu(T && coefficient)
        {
            using CoefficientType = typename std::decay<T>::type;
            static_assert(std::is_base_of<IValue, CoefficientType>::value,
                          "DDu coefficient must derive from IValue");
            m_DDuFunctions.emplace_back(
              std::make_unique<CoefficientType>(std::forward<T>(coefficient)), nullptr);
        }

        //! Template function that will add K * Du/Dt matrix into the system.
        template<typename T>
        void Cap(T & t,
                 typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_CapacitanceFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        template<typename T>
        void Cap(T && t,
                 typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_CapacitanceFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        //! Template function that will add K * Du/Dt matrix into the system with the
        //! capacity lumped nodally (v_i * c_i on the diagonal) regardless of
        //! m_LumpCapacityNodally. Required for near-singular capacities such as the
        //! latent-heat-of-fusion secant: the consistent pairwise-average integration
        //! smears the spike onto neighbouring rows, which then book phantom storage
        //! energy against their own temperature change.
        template<typename T>
        void CapNodal(T && t,
                      typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * =
                        nullptr)
        {
            m_NodalCapacitanceFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        //! Template function that will add K*((Dp/Dx)(Du/Dx) + (Dp/Dy)(Du/Dy)) matrix into the
        //! system.
        template<typename T, typename U>
        void DpDu(T & t,
                  U & u,
                  typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_DpDuFunctions.emplace_back(std::unique_ptr<T>(new T(t)),
                                         std::unique_ptr<U>(new U(u)));
        }

        template<typename T, typename U>
        void DpDu(T && t,
                  U & u,
                  typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_DpDuFunctions.emplace_back(std::unique_ptr<T>(new T(t)),
                                         std::unique_ptr<U>(new U(u)));
        }

        template<typename T, typename U>
        void DpDu(T & t,
                  U && u,
                  typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_DpDuFunctions.emplace_back(std::unique_ptr<T>(new T(t)),
                                         std::unique_ptr<U>(new U(u)));
        }

        template<typename T, typename U>
        void DpDu(T && t,
                  U && u,
                  typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_DpDuFunctions.emplace_back(std::unique_ptr<T>(new T(t)),
                                         std::unique_ptr<U>(new U(u)));
        }

        //! Transport of a PRODUCT potential, div(coefficient grad(factor * u)), assembled as
        //! the coefficient's stiffness matrix with column j scaled by the nodal factor at
        //! node j:
        //!
        //!     K_ij = S_ij(coefficient) * factor_j
        //!
        //! The general form of the transport term; the one-argument DDu above is this with
        //! factor = 1.
        //!
        //! Used for the moisture vapour term, where the transported potential is the vapour
        //! content c_sat(T) * phi while the state variable is phi. Assembling it this way,
        //! rather than as a stiffness in phi plus a separately integrated coupling term in
        //! c_sat, gives two properties that the split form does not have together:
        //!
        //!   - equilibrium is exact: when the nodal products factor_j * u_j are all equal the
        //!     result is factor * sum_j S_ij = 0, since a stiffness matrix has vanishing row
        //!     sums (sum_j psi_j = 1). A sealed strip settles ON its closed-form profile on
        //!     ANY mesh instead of converging to it as the mesh is refined;
        //!   - moisture is conserved: sum_i K_ij = factor_j * sum_i S_ij = 0 by the same
        //!     argument applied to the columns, for any spatially varying coefficient.
        //!
        //! See doc/Moisture Governing Equations.md (D1) and
        //! tst/units/validation/SealedStrip_SteadyGradient.unit.cxx.
        template<typename T, typename U>
        void DDu(T && coefficient, U && nodalFactor)
        {
            using CoefficientType = typename std::decay<T>::type;
            using FactorType = typename std::decay<U>::type;
            static_assert(std::is_base_of<IValue, CoefficientType>::value,
                          "DDu coefficient must derive from IValue");
            static_assert(std::is_base_of<IValue, FactorType>::value,
                          "DDu nodal factor must derive from IValue");
            m_DDuFunctions.emplace_back(
              std::make_unique<CoefficientType>(std::forward<T>(coefficient)),
              std::make_unique<FactorType>(std::forward<U>(nodalFactor)));
        }

        //! Template function that will create functions used in equivalent material conductivity
        //! and therefore used in flux calculations.
        template<typename T>
        void CondFlux(T & t,
                  typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_ConductanceFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        template<typename T>
        void CondFlux(T && t,
                  typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_ConductanceFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        //! \brief Supports multiplication of between matrix and vector in differential equations.
        //!
        //! Differential equations sometimes have mutliplication between matrix of that is function
        //! of certain node properties and vector of some node property as well. For example,
        //! heat transfer equation have [H]*{vapor content} that needs to be evaluated. This
        //! structure is used to support matrix-vector pair that will be used in differential
        //! equations creation.
        struct MatrixVector
        {
            //! Simple constructor for matrix-vector pair. It contains function that will be used to
            //! create matrix and vector that will be created from some of node properties.
            MatrixVector(
              iValue && MatrixFunction,   //!< Function that will be used to create matrix
              Variable PropertyVector     //!< Node property that will be used to create vector
            );

            //! Matrix-vector pair whose vector is a computed per-node function instead of a raw
            //! node property (e.g. the vapour flux potential c_sat * phi).
            MatrixVector(
              iValue && MatrixFunction,   //!< Function that will be used to create matrix
              iValue && VectorFunction    //!< Function evaluated per node to create the vector
            );
            iValue MatrixFunction;
            Variable PropertyVector{Variable::temperature};
            iValue VectorFunction;   //!< When set, takes precedence over PropertyVector.
        };

        //! Template function that will create necessary multiplier between matrix and vector.
        template<typename T>
        void multiplies(
          T & t,
          Variable property,
          const typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_Matrix_x_Vector.emplace_back(std::unique_ptr<T>(new T(t)), property);
        }

        //! Multiplier whose vector side is a computed per-node function (any IValue).
        template<typename T, typename U>
        void multiplies(
          T && t,
          U && u,
          const typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr,
          const typename std::enable_if<std::is_base_of<IValue, U>::value, U>::type * = nullptr)
        {
            m_Matrix_x_Vector.emplace_back(iValue(new T(t)), iValue(new U(u)));
        }

        //! Template function that will create necessary multiplier between matrix and vector.
        template<typename T>
        void multiplies(
          T && t,
          Variable property,
          const typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = nullptr)
        {
            m_Matrix_x_Vector.emplace_back(std::unique_ptr<T>(new T(t)), property);
        }

        const IMaterial & m_Material;
        Variable m_FluxVariable;   // State variable used to calculate flux

        //! When true, the capacitance is lumped nodally as diag(v_i * xi_i) with v_i the
        //! nodal volume integral(psi_i) and xi_i the nodal capacity value, instead of the
        //! coefficient-averaged consistent mass integral(xi psi_i psi_j) later row-lumped.
        //! With the secant moisture capacity this makes C*(phi - phi_prev) telescope exactly
        //! to the change in stored moisture, giving machine-precision conservation (D6). Set
        //! by the moisture element; thermal keeps the default (false) so its mass is unchanged.
        bool m_LumpCapacityNodally{false};

        //! When true, DDu/DpDu coefficients are interpolated at each Gauss point via the
        //! shape functions (IQLEIntegrator2D::integrateInterpolated) instead of the pairwise
        //! nodal average. The pairwise average is non-conservative for coefficients that vary
        //! inside an element (e.g. delta*c_sat(T) under a temperature gradient) and was the
        //! source of the ~1% closed-strip moisture creation in the D1 diagnostic. Set by the
        //! moisture element; thermal keeps the default (false) so its assembly is unchanged.
        bool m_InterpolateCoefficientsAtGaussPoints{false};

        //! Constant volumetric source over the element; zero contributes nothing.
        double m_VolumetricSource{0.0};

    private:
        //! \brief Used to evaluate complex differential equation with two variable functions. One
        //! standing as multiplicator and other one as part of derivative. Equation k*(Dp/Dx)(Du/Dx)
        //! + k*(Dp/Dy)(Du/Dy) have fixed function (k) and independent function (p) that is part of
        //! derivative. Note that u is state variable for which matrix equations will be formed.
        struct DerivativeFunction
        {
            DerivativeFunction(
              std::unique_ptr<IValue>
                fixedValue,   //!< Stands as free function (equal to k from above comment)
              std::unique_ptr<IValue>
                derivativeValue   //!< Stands as derivative function (equal to p from above comment)
            );

            std::unique_ptr<IValue> fixedValue;
            std::unique_ptr<IValue> derivativeValue;
        };

        //! \brief Coefficient and nodal factor of a product-form transport term, see the
        //! two-argument DDu. The coefficient is integrated into a stiffness matrix; the nodal
        //! factor scales its columns.
        struct NodalProductFunction
        {
            NodalProductFunction(
              std::unique_ptr<IValue> coefficient,   //!< Multiplies the gradient (e.g. delta)
              std::unique_ptr<IValue> nodalFactor    //!< Turns the state into the transported
                                                     //!< potential (e.g. c_sat)
            );

            std::unique_ptr<IValue> coefficient;
            std::unique_ptr<IValue> nodalFactor;
        };

        std::vector<NodalProductFunction> m_DDuFunctions;
        std::vector<iValue> m_ConductanceFunctions;
        std::vector<iValue> m_CapacitanceFunctions;
        //! Capacitance functions that are always lumped nodally (see CapNodal).
        std::vector<iValue> m_NodalCapacitanceFunctions;
        std::vector<DerivativeFunction> m_DpDuFunctions;

        //! Vector of values that will simply be evaluated on right hand side.
        //! This is in form [M]*{V} (Matrix * vector). First property is simply set of functions
        //! that form matrix and second is simply property of vectors.
        std::vector<MatrixVector> m_Matrix_x_Vector;

        //! Calculates angle between two vectors made of (node1-node2) and (node1-node3)
        double angleBetweenNodes(const Node2D & node1, const Node2D & node2, const Node2D & node3) const;

        //! State-variable gradients (d/dx, d/dy) at each Gauss point. Shared building block of
        //! flux() and fluxAtGaussPoints().
        std::vector<std::array<double, 2>> stateGradientsAtGaussPoints() const;

        //! Exact nodal-lumped capacity diag(v_i * nodalCapacity_i), v_i = integral(psi_i).
        //! Used when m_LumpCapacityNodally is set (moisture). See m_LumpCapacityNodally.
        [[nodiscard]] ElementMatrix2D
          nodalLumpedCapacity(std::span<const double> nodalCapacity) const;


        //! Circular vector connects first and last node so that program can easily iterate
        //! through nodes. Vector iterator will skip identical nodes (in trianglular element
        //! for example). Angle between nodes algorithm is using this vector.
        class CircularNodesVector : public INodes
        {
        public:
            explicit CircularNodesVector(
              const std::initializer_list<std::reference_wrapper<INode2D>> & __l) :
                INodes(__l)
            {}

            explicit CircularNodesVector(Node2D & node1, Node2D & node2) : INodes(node1, node2)
            {}

            explicit CircularNodesVector(Node2D & node1,
                                         Node2D & node2,
                                         Node2D & node3,
                                         Node2D & node4) :
                INodes(node1, node2, node3, node4)
            {}

            [[nodiscard]] std::vector<std::reference_wrapper<INode2D>>::const_iterator begin() const
            {
                return m_Nodes.begin();
            }

            [[nodiscard]] std::vector<std::reference_wrapper<INode2D>>::const_iterator end() const
            {
                return m_Nodes.end();
            }

            /// Keeps iterating over unique elements of the vector
            const INode2D & current()
            {
                return m_Nodes[currentIndex];
            }

            bool last() const
            {
                return passedLast;
            }

            const INode2D & previous()
            {
                auto validIndex = checkPrevIndex(currentIndex);

                while(m_Nodes[validIndex] == m_Nodes[currentIndex])
                {
                    validIndex = checkPrevIndex(validIndex);
                }

                return m_Nodes[validIndex];
            }

            const INode2D & next()
            {
                auto validIndex = checkNextIndex(currentIndex);
                while(m_Nodes[validIndex] == m_Nodes[currentIndex])
                {
                    validIndex = checkNextIndex(validIndex);
                }

                return m_Nodes[validIndex];
            }

            void moveToNext()
            {
                auto nextIndex = checkNextIndex(currentIndex);
                while(m_Nodes[nextIndex] == m_Nodes[currentIndex])
                {
                    nextIndex = checkNextIndex(nextIndex);
                }
                passedLast = nextIndex < currentIndex;
                currentIndex = nextIndex;
            }

        private:
            std::size_t checkNextIndex(const std::size_t index) const
            {
                auto validIndex = index;

                if(validIndex != m_Nodes.size() - 1)
                {
                    validIndex = index + 1;
                }
                else
                {
                    validIndex = 0;
                }

                return validIndex;
            }

            std::size_t checkPrevIndex(const std::size_t index) const
            {
                auto validIndex = index;

                if(validIndex != 0)
                {
                    validIndex = index - 1;
                }
                else
                {
                    validIndex = m_Nodes.size() - 1;
                }

                return validIndex;
            }
            std::size_t currentIndex{0};
            bool passedLast{false};
        };

        CircularNodesVector m_Nodes;

        QuadrilateralLinearGlobal2D m_Global2D;
        QLECapacitanceIntegrator2D m_QLECapacitance2D;

        //! Integration matrix is derived purely from the (fixed) element geometry, so it is built
        //! once here rather than on every assembly.
        QLEDDuIntegrator2D m_QLEDDu2D;

        //! setIndependentVariables() overwrites every entry before each use, so the integrator is
        //! reused instead of reconstructed. Mutable because the assembly accessors are const.
        mutable QLEDpDuIntegrator2D m_QLEDpDu2D;

        //! Nodal volumes integral(psi_i). Geometry only, so computed once at construction --
        //! nodalLumpedCapacity() previously rebuilt them on every call, which made it the single
        //! most expensive thing in the assembly.
        const std::array<double, numOfQuadrilateralNodes> m_NodalVolumes;

        const bool m_Linear;
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementThermalLinear2D
    //////////////////////////////////////////////////////////////////////////////

    //! Creates element that will be used in thermal domain to find temperature distribution
    //! solution.
    class ElementThermalLinear2D : public IElementLinear2D
    {
    public:
        ElementThermalLinear2D(
          Nodes & nodePool,               //!< Reference to NodePool for node lookup
          Materials & materialPool,          //!< Reference to MaterialPool for material lookup
          size_t index1,                     //!< Node 1 index
          size_t index2,                     //!< Node 2 index
          size_t index3,                     //!< Node 3 index
          size_t index4,                     //!< Node 4 index
          const std::string & materialName   //!< SolidMaterial name assigned to the element
        );
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementMoistureLinear2D
    //////////////////////////////////////////////////////////////////////////////

    //! Creates element that will be used in moisture domain to find water distribution solution
    class ElementMoistureLinear2D : public IElementLinear2D
    {
    public:
        ElementMoistureLinear2D(
          Nodes & nodePool,               //!< Reference to NodePool for node lookup
          Materials & materialPool,          //!< Reference to MaterialPool for material lookup
          size_t index1,                     //!< Node 1 index
          size_t index2,                     //!< Node 2 index
          size_t index3,                     //!< Node 3 index
          size_t index4,                     //!< Node 4 index
          const std::string & materialName   //!< SolidMaterial name assigned to the element
        );
    };

}   // namespace HygroThermFEM
