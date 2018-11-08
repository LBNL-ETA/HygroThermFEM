#pragma once

#include <memory>

#include "Functions.hxx"
#include "Material.hxx"
#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM
{
    // Constant that holds number of nodes in certain elementsCreator
    const std::size_t numOfQuadrilateralNodes = 4;
    // const std::size_t numOfIntegrationPoints = 4;

    //////////////////////////////////////////////////////////////////////////////
    ///  IQLEMatrix2D
    //////////////////////////////////////////////////////////////////////////////

    // Abstract class that forces users to perform matrix calculation in inherited
    // version Depending on equations, matrices will have different calculation
    // methods (for example capacitance and conductance matrices have different
    // form)
    class IQLEIntegrator2D
    {
    public:
        IQLEIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);

        // Integrate matrix over all points of integration
        virtual FenestrationCommon::SquareMatrix
          integrate(const std::vector<double> & t_Values) const final;

    protected:
        virtual void
          calculateMatrixInIntegrationPoint(FenestrationCommon::SquareMatrix & matrix,
                                            const std::vector<double> & t_Values,
                                            const std::size_t t_IntegrationPointIndex) const final;

        const QuadrilateralLinearGlobal2D & m_Global2D;

        std::vector<FenestrationCommon::SquareMatrix> m_DPsiDxDyMatrix;
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDDuIntegrator2D
    //////////////////////////////////////////////////////////////////////////////

    // Class to handle conductance matrix in global coordinate system
    // Conductance equation is D/Dx(Du/Dx)+D/Dy(Du/Dy) where u is state variable
    class QLEDDuIntegrator2D : public IQLEIntegrator2D
    {
    public:
        virtual ~QLEDDuIntegrator2D() = default;

        QLEDDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDpDuIntegrator2D
    //////////////////////////////////////////////////////////////////////////////

    // Class to handle double derivative matrix in global coordinate system
    // Equation is (Dp/Dx)(Du/Dx)+(Dp/Dy)(Du/Dy) where u is state variable and
    // p is independent variable that can have different values in different nodes.
    class QLEDpDuIntegrator2D : public IQLEIntegrator2D
    {
    public:
        virtual ~QLEDpDuIntegrator2D() = default;

        QLEDpDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element);

        // Update independent variables (corresponds to variable p in above equation)
        void setIndependentVariables(const std::vector<double> & t_Values);

        void clearIntegrationMatrix();
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  QLECapacitance2D
    //////////////////////////////////////////////////////////////////////////////

    // Class to handle capacitance matrix in global coordinate system
    class QLECapacitance2D : public IQLEIntegrator2D
    {
    public:
        virtual ~QLECapacitance2D() = default;

        QLECapacitance2D(const QuadrilateralLinearGlobal2D & t_Element);
    };

    /// Equation k*(Dp/Dx)(Du/Dx) + k*(Dp/Dy)(Du/Dy) have fixed function (k) and
    /// independent function (p) that is part of derivative. Note that u is state
    /// variable for which matrix equations will be formed.
    struct DerivativeFunction
    {
        DerivativeFunction(std::unique_ptr<IValue> fixedValue,
                           std::unique_ptr<IValue> derivativeValue);

        std::unique_ptr<IValue> fixedValue;
        std::unique_ptr<IValue> derivativeValue;
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D
    //////////////////////////////////////////////////////////////////////////////

    /// Class that handles creation of conductance and capacitance matrices in
    /// linear 2D world. This class will be used by multiple governing equations
    /// since basis of matrix creation are identical with only difference in what
    /// coefficients are passed
    class IElementLinear2D
    {
    public:
        IElementLinear2D(const Node2D & t_Node1,
                         const Node2D & t_Node2,
                         const Node2D & t_Node3,
                         const Node2D & t_Node4,
                         const Material & t_Material,
                         const bool isLinear = true);

        FenestrationCommon::SquareMatrix DDuMatrices() const;

        FenestrationCommon::SquareMatrix DpDuMatrices();

        FenestrationCommon::SquareMatrix capacitanceMatrix() const;

        Node2D & getNode(std::size_t index);

        bool haveBothNodes(const Node2D & t_Node1, const Node2D & t_Node2) const;

        std::vector<std::size_t> nodeIndexes() const;

        const Material & getMaterial() const;

        virtual bool isLinear() const final;

    protected:
        std::vector<iValue> m_DDuFunctions;
        std::vector<iValue> m_Capacitance;
        std::vector<DerivativeFunction> m_DuDpFunctions;

        const Material & m_Material;

    private:
        /// calculates angle between two vectors made of (node1-node2) and (node1-node3)
        double angleBetweenNodes(const Node2D & node1, const Node2D & node2, const Node2D & node3);


        /// Circular vector connects first and last node so that program can easily iterate
        /// through nodes. Vector iterator will skip identical nodes (in trianglular element
        /// for example). Angle between nodes algorithm is using this vector.
        class CircularNodesVector : public INodes
        {
        public:
            explicit CircularNodesVector(const std::initializer_list<Node2D> & __l) :
                INodes(__l),
                currentIndex(0),
                passedLast(false)
            {}

            typename std::vector<Node2D>::const_iterator begin() const
            {
                return m_Nodes.begin();
            }

            typename std::vector<Node2D>::const_iterator end() const
            {
                return m_Nodes.end();
            }

            /// Keeps iterating over unique elements of the vector
            Node2D & current()
            {
                return m_Nodes[currentIndex];
            }

            bool last()
            {
                return passedLast;
            }

            Node2D & previous()
            {
                auto validIndex = checkPrevIndex(currentIndex);

                while(m_Nodes[validIndex] == m_Nodes[currentIndex])
                {
                    validIndex = checkPrevIndex(validIndex);
                }

                return m_Nodes[validIndex];
            }

            Node2D & next()
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
        QLECapacitance2D m_QLECapacitance2D;
        QLEDDuIntegrator2D m_DDuIntegrator;
        /// This one depends on functions and must be stored for every submatrix
        std::vector<QLEDpDuIntegrator2D> m_QLEDpDuIntegrator2D;

        const bool m_Linear;
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementThermalLinear2D
    //////////////////////////////////////////////////////////////////////////////

    class ElementThermalLinear2D : public IElementLinear2D
    {
    public:
        ElementThermalLinear2D(const Node2D & t_Node1,
                               const Node2D & t_Node2,
                               const Node2D & t_Node3,
                               const Node2D & t_Node4,
                               const Material & mat);
    };

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementMoistureLinear2D
    //////////////////////////////////////////////////////////////////////////////

    class ElementMoistureLinear2D : public IElementLinear2D
    {
    public:
        ElementMoistureLinear2D(const Node2D & t_Node1,
                                const Node2D & t_Node2,
                                const Node2D & t_Node3,
                                const Node2D & t_Node4,
                                const Material & mat);
    };

}   // namespace MoisThermFEM
