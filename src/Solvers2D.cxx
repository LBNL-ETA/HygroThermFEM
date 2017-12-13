#include "SquareMatrix.hxx"
#include "LinearSolver.hxx"
#include "Solvers2D.hxx"
#include "NodePool.hxx"
#include "Elements2D.hxx"
#include "BoundaryConditions2D.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

    /////////////////////////////////////////////////////////////////////////
    // ISolver2D
    /////////////////////////////////////////////////////////////////////////

    ISolver2D::ISolver2D(Elements2DLinear const &t_Elements, BoundaryConditions2D const &t_BCs) :
            m_Elements(t_Elements), m_BCs(t_BCs), m_Solved(false) {

    }

    /////////////////////////////////////////////////////////////////////////
    // SteadyStateSolver2D
    /////////////////////////////////////////////////////////////////////////

    SteadyStateSolver2D::SteadyStateSolver2D(
            Elements2DLinear const &t_Elements,
            BoundaryConditions2D const &t_BCs) :
            ISolver2D(t_Elements, t_BCs) {

    }

    std::vector<double> SteadyStateSolver2D::solution() {
        if (!m_Solved) {
            solve();
            m_Solved = true;
        }

        return m_Solution;
    }

    void SteadyStateSolver2D::solve() {
        CSquareMatrix condMat = m_Elements.thermalConductivity();
        CSquareMatrix H = m_BCs.matrixA();

        std::vector< double > B = m_BCs.vectorR();
        condMat = condMat.add( H );

        CLinearSolver aSolver;
        m_Solution = aSolver.solveSystem(condMat, B);

    }

    /////////////////////////////////////////////////////////////////////////
    // TransientSolver2D
    /////////////////////////////////////////////////////////////////////////

    TransientSolver2D::TransientSolver2D(
            Elements2DLinear const &t_Elements,
            BoundaryConditions2D const &t_BCs,
            double const DTemp,
            size_t const NTimeSteps) :
            ISolver2D(t_Elements, t_BCs), m_DTemp(DTemp), m_NSteps(NTimeSteps) {

    }

    std::vector<std::vector<double> > TransientSolver2D::solution() {
        if (!m_Solved) {
            solve();
            m_Solved = true;
        }

        return m_Solution;
    }

    void TransientSolver2D::solve() {
        auto size = NodePool::Instance().maxIndex();

        // ublas::compressed_matrix< double, ublas::column_major, 0 > A( size, size, 3 * size );
        // ublas::permutation_matrix< size_t > pm( A.size1() );
        // ublas::vector< double > y( size );

        // vector< vector< double > > condMat = m_Elements.thermalConductivity();
        // vector< vector< double > > H = m_BCs.matrixA();
        auto RhoCp = m_Elements.rhoCp();
        std::vector<double> M(size);

        // Creates lump matrix
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                M[i] += RhoCp[i][j];
            }
            M[i] /= m_DTemp;
        }

        auto rBC = m_BCs.vectorR();

        auto condMat = m_Elements.thermalConductivity();
				CSquareMatrix H { m_BCs.matrixA() };

        condMat = condMat.add( H );
        condMat = condMat.addDiagonal(M);

        auto temperatures = NodePool::Instance().nodeTemperatures();
        std::vector<double> B(temperatures.size());


        for (unsigned i = 0; i < m_NSteps; ++i) {
            for (unsigned j = 0; j < size; ++j) {
                B[j] = temperatures[j] * M[j] + rBC[j];
            }

            CLinearSolver aSolver;
            B = aSolver.solveSystem( condMat, B );

            std::vector<double> aSolution(size);
            for (unsigned j = 0; j < size; ++j) {
                aSolution[j] = B[j];
                temperatures[j] = B[j];
            }
            m_Solution.push_back(aSolution);
        }
    }

}