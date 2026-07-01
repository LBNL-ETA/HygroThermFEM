#include "BoundaryConditions2D.hxx"

namespace HygroThermFEM
{
    /// BoundaryConditions2D::BoundaryConditions2D(
    /// 		std::vector< std::reference_wrapper< IBCLinear2D > > & t_BCs ) :
    /// 		m_BCs( t_BCs ), m_Linear( true ) {
    /// 	for ( const auto & aBc : m_BCs ) {
    /// 		m_Linear = m_Linear && aBc.get().isLinear();
    /// 	}
    /// }

    SquareMatrix BoundaryConditions2D::HMatrix(const size_t maxNodeIndex,
                                               const size_t timestepIndex) const
    {
        // Assemble directly into a triplet list. The H-matrix has only the boundary-node 2x2 blocks
        // populated, so the previous dense maxNodeIndex x maxNodeIndex intermediate allocated and
        // scanned O(n^2) storage for a handful of nonzeros -- setFromTriplets sums the duplicates.
        std::vector<Eigen::Triplet<double>> tripletList;
        forEachActiveBC(timestepIndex, [&](const IBCLinear2D & bcItem)
        {
            const auto indexes = bcItem.getNodeIndexes();
            const auto matH = bcItem.H_Matrix();
            for(size_t row = 0; row < 2; ++row)
            {
                for(size_t col = 0; col < 2; ++col)
                {
                    tripletList.emplace_back(static_cast<int>(indexes[row] - 1),
                                             static_cast<int>(indexes[col] - 1),
                                             matH(row, col));
                }
            }
        });
        return SquareMatrix{maxNodeIndex, tripletList};
    }

    std::vector<double> BoundaryConditions2D::RVector(const size_t maxNodeIndex,
                                                      const size_t timestepIndex) const
    {
        std::vector<double> result(maxNodeIndex, 0);
        forEachActiveBC(timestepIndex, [&](const IBCLinear2D & bcItem)
        {
            const auto indexes = bcItem.getNodeIndexes();
            const auto vecR = bcItem.R_Vector();
            for(size_t idx = 0; idx < 2; ++idx)
            {
                result[indexes[idx] - 1] += vecR[idx];
            }
        });
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

        m_TransientBCs.push_back(std::move(bc));
    }

    void BoundaryConditions2D::clear()
    {
        m_BCs.clear();
        m_TransientBCs.clear();
    }
}   // namespace HygroThermFEM
