#pragma once

#include <memory>

#include "Functions.hxx"
#include "Material.hxx"
#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM
{
    //! Constant that holds number of nodes in certain elementsCreator
    const std::size_t numOfQuadrilateralNodes = 4;

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
        IQLEIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);

        //! Integrate matrix over all points of integration
        virtual FenestrationCommon::SquareMatrix
          integrate(const std::vector<double> &
                      t_Values   //!< Nodal values for which integration will be performed
                    ) const final;

    protected:
        //! Helper function that integrates matrix in given integration point. Integration is
        //! identical with no matter what integration matrix looks like.
        virtual void calculateMatrixInIntegrationPoint(
          FenestrationCommon::SquareMatrix &
            matrix,   //!< Matrix that integration results will be added to
          const std::vector<double> &
            t_Values,   //!< Nodal values for which integration will be performed
          const std::size_t t_IntegrationPointIndex   //!< Integration point
          ) const final;

        const QuadrilateralLinearGlobal2D & m_Global2D;

        //! This matrix will hold different forms of shape functions operations. This will mainly
        //! depend on what integrator will be used to integrate matrices.
        std::vector<FenestrationCommon::SquareMatrix> m_IntegrationMatrix;
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
        void setIndependentVariables(const std::vector<double> & t_Values);
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
        //! Construction of base element. Constructor will make all necessary shape functions in
        //! background (since they are dependent on element geometry). In constructor will
        //! calculate weighting factors that will be used to calculate some local node properties.
        //! This is important since one node can belong to several different elements and when
        //! doing averaging, program will have information about how much of weight one node
        //! this element will have.
        IElementLinear2D(
          const size_t index1,                //!< Node 1 index
          const size_t index2,                //!< Node 2 index
          const size_t index3,                //!< Node 3 index
          const size_t index4,                //!< Node 4 index
          const std::string & materialName,   //!< Material assigned to the element
          const bool isLinear = true   //!< States if element equations can be solved by linear
                                       //!< approach or non-linear iterations needs to be used.
        );

        //! Integrates all matrices that are part of K * (D/Dx(Du/Dx) + D/Dy(Du/Dy)) equation.
        FenestrationCommon::SquareMatrix DDuMatrices() const;

        //! Integrates all matrices that are part of K * ((Dp/Dx)(Du/Dx) + (Dp/Dy)(Du/Dy)) equation.
        FenestrationCommon::SquareMatrix DpDuMatrices();

        //! Integrates all matrices that are part of K * Du/Dt equation.
        FenestrationCommon::SquareMatrix capacitanceMatrices() const;

        //! Integrates right hand-side vector.
        std::vector<double> rightSideVector() const;

        std::vector<double> flux() const;

        Node2D & getNode(std::size_t index);

        //! Returns true of element have both nodes. This is used to test if certain boundary
        //! condition is part of the element.
        bool haveBothNodes(const size_t index1,   //!< Node 1 index
                           const size_t index2    //!< Node 2 index
                           ) const;

        //! Returns indexes of all four nodes in the element
        std::vector<std::size_t> nodeIndexes() const;

        //! Returns material that is assigned to the element.
        const Material & getMaterial() const;

        //! Returns if element satisfies linear problem.
        virtual bool isLinear() const final;

    protected:
        //! Template function that will add DDu matrix into the system.
        template<typename T>
        void DDu(T & t,
                 const typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = 0)
        {
            m_DDuFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        //! Template function that will add K*(D/Dx(Du/Dx) + D/Dy(Du/Dy)) matrix into the system.
        template<typename T>
        void DDu(T && t,
                 const typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = 0)
        {
            m_DDuFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        //! Template function that will add K * Du/Dt matrix into the system.
        template<typename T>
        void Cap(T & t, typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = 0)
        {
            m_CapacitanceFunctions.emplace_back(std::unique_ptr<T>(new T(t)));
        }

        //! Template function that will add K*((Dp/Dx)(Du/Dx) + (Dp/Dy)(Du/Dy)) matrix into the
        //! system.
        template<typename T, typename U>
        void DpDu(T & t,
                  U & u,
                  typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = 0)
        {
            m_DpDuFunctions.emplace_back(std::unique_ptr<T>(new T(t)),
                                         std::unique_ptr<U>(new U(u)));
        }

        //! Template function that will create functions used in equivalent material conductivity
        //! and therefore used in flux calculations.
        template<typename T>
        void Cond(T & t, typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = 0)
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
              iValue && MatrixFunction,       //!< Function that will be used to create matrix
              const Variable PropertyVector   //!< Node property that will be used to create vector
            );
            iValue MatrixFunction;
            Variable PropertyVector;
        };

        //! Template function that will create necessary multiplier between matrix and vector.
        template<typename T>
        void multiplies(
          T & t,
          Variable property,
          const typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = 0)
        {
            m_Matrix_x_Vector.emplace_back(std::unique_ptr<T>(new T(t)), property);
        }

        //! Template function that will create necessary multiplier between matrix and vector.
        template<typename T>
        void multiplies(
          T && t,
          Variable property,
          const typename std::enable_if<std::is_base_of<IValue, T>::value, T>::type * = 0)
        {
            m_Matrix_x_Vector.emplace_back(std::unique_ptr<T>(new T(t)), property);
        }

        const Material & m_Material;

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

        std::vector<iValue> m_DDuFunctions;
        std::vector<iValue> m_ConductanceFunctions;
        Variable m_FluxVariable;   // State variable used to calculate flux
        std::vector<iValue> m_CapacitanceFunctions;
        std::vector<DerivativeFunction> m_DpDuFunctions;

        //! Vector of values that will simply be evaluated on right hand side.
        //! This is in form [M]*{V} (Matrix * vector). First property is simply set of functions
        //! that form matrix and second is simply property of vectors.
        std::vector<MatrixVector> m_Matrix_x_Vector;

        //! Calculates angle between two vectors made of (node1-node2) and (node1-node3)
        double angleBetweenNodes(const Node2D & node1, const Node2D & node2, const Node2D & node3);


        //! Circular vector connects first and last node so that program can easily iterate
        //! through nodes. Vector iterator will skip identical nodes (in trianglular element
        //! for example). Angle between nodes algorithm is using this vector.
        class CircularNodesVector : public INodes
        {
        public:
            explicit CircularNodesVector(
              const std::initializer_list<std::reference_wrapper<Node2D>> & __l) :
                INodes(__l),
                currentIndex(0),
                passedLast(false)
            {}

            explicit CircularNodesVector(Node2D & node1, Node2D & node2) :
                INodes(node1, node2),
                currentIndex(0),
                passedLast(false)
            {}

            explicit CircularNodesVector(Node2D & node1,
                                         Node2D & node2,
                                         Node2D & node3,
                                         Node2D & node4) :
                INodes(node1, node2, node3, node4),
                currentIndex(0),
                passedLast(false)
            {}

            typename std::vector<std::reference_wrapper<Node2D>>::const_iterator begin() const
            {
                return m_Nodes.begin();
            }

            typename std::vector<std::reference_wrapper<Node2D>>::const_iterator end() const
            {
                return m_Nodes.end();
            }

            /// Keeps iterating over unique elements of the vector
            const Node2D & current()
            {
                return m_Nodes[currentIndex];
            }

            bool last()
            {
                return passedLast;
            }

            const Node2D & previous()
            {
                auto validIndex = checkPrevIndex(currentIndex);

                while(m_Nodes[validIndex] == m_Nodes[currentIndex])
                {
                    validIndex = checkPrevIndex(validIndex);
                }

                return m_Nodes[validIndex];
            }

            const Node2D & next()
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
        QLEDDuIntegrator2D m_DDuIntegrator;

        /// This one depends on independent variables and should be stored in vector
        std::vector<QLEDpDuIntegrator2D> m_QLEDpDuIntegrator2D;

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
          const size_t index1,               //!< Node 1 index
          const size_t index2,               //!< Node 2 index
          const size_t index3,               //!< Node 3 index
          const size_t index4,               //!< Node 4 index
          const std::string & materialName   //!< Material name assigned to the element
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
          const size_t index1,               //!< Node 1 index
          const size_t index2,               //!< Node 2 index
          const size_t index3,               //!< Node 3 index
          const size_t index4,               //!< Node 4 index
          const std::string & materialName   //!< Material name assigned to the element
        );
    };

}   // namespace MoisThermFEM
