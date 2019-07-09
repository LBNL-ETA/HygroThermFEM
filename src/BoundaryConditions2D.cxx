#include "BoundaryConditions2D.hxx"
#include "NodePool.hxx"

namespace HygroThermFEM
{
    /// BoundaryConditions2D::BoundaryConditions2D(
    /// 		std::vector< std::reference_wrapper< IBCLinear2D > > & t_BCs ) :
    /// 		m_BCs( t_BCs ), m_Linear( true ) {
    /// 	for ( const auto & aBc : m_BCs ) {
    /// 		m_Linear = m_Linear && aBc.get().isLinear();
    /// 	}
    /// }

    SquareMatrix BoundaryConditions2D::HMatrix() const
    {
        std::vector<std::vector<double>> result{
          NodePool::Instance().maxIndex(), std::vector<double>(NodePool::Instance().maxIndex(), 0)};
        for(const std::unique_ptr<IBCLinear2D> & aBc : m_BCs)
        {
            auto indexes = aBc->getNodeIndexes();
            auto matH = aBc->H_Matrix();
            for(size_t i = 0; i < 2; ++i)
            {
                for(size_t j = 0; j < 2; ++j)
                {
                    result[indexes[i] - 1][indexes[j] - 1] += matH(i, j);
                }
            }
        }

        //TODO: Get matrices from vector
        return SquareMatrix{result};
    }

    std::vector<double> BoundaryConditions2D::RVector() const
    {
        std::vector<double> result(NodePool::Instance().maxIndex(), 0);
        // Create full size matrices
        for(const std::unique_ptr<IBCLinear2D> & aBc : m_BCs)
        {
            auto indexes = aBc->getNodeIndexes();
            auto vecR = aBc->R_Vector();
            for(size_t i = 0; i < 2; ++i)
            {
                result[indexes[i] - 1] += vecR[i];
            }
        }

        //TODO: Get vectors from timestep
        return result;
    }

    bool BoundaryConditions2D::isLinear() const
    {
        return m_Linear;
    }

    BoundaryConditions2D::BoundaryConditions2D() : m_Linear(true)
    {}

    void BoundaryConditions2D::assignBC(std::unique_ptr<IBCLinear2D> && bc)
    {
        m_Linear = m_Linear && bc->isLinear();
        m_BCs.push_back(std::move(bc));
    }

    void BoundaryConditions2D::assignTimestepBCs(std::vector<std::unique_ptr<IBCLinear2D>> && bc)
    {
        std::for_each(
          bc.begin(), bc.end(), [&](const auto & aBC) { m_Linear = m_Linear && aBC->isLinear(); });

        m_TransientBCs.push_back(bc);
    }
}   // namespace HygroThermFEM
