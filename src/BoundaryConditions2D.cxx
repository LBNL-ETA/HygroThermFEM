#include "BoundaryConditions2D.hxx"
#include "NodePool.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

/// BoundaryConditions2D::BoundaryConditions2D(
/// 		std::vector< std::reference_wrapper< IBCLinear2D > > & t_BCs ) :
/// 		m_BCs( t_BCs ), m_Linear( true ) {
/// 	for ( const auto & aBc : m_BCs ) {
/// 		m_Linear = m_Linear && aBc.get().isLinear();
/// 	}
/// }

SquareMatrix BoundaryConditions2D::HMatrix() const {
  std::vector<std::vector<double>> result{
      NodePool::Instance().maxIndex(),
      std::vector<double>(NodePool::Instance().maxIndex(), 0)};
  for (const std::unique_ptr<IBCLinear2D> &aBc : m_BCs) {
    auto indexes = aBc->getNodeIndexes();
    auto matH = aBc->H_Matrix();
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        result[indexes[i] - 1][indexes[j] - 1] += matH(i, j);
      }
    }
  }
  return SquareMatrix{result};
}

std::vector<double> BoundaryConditions2D::RVector() const {
  std::vector<double> result(NodePool::Instance().maxIndex(), 0);
  // Create full size matrices
  for (const std::unique_ptr<IBCLinear2D> &aBc : m_BCs) {
    auto indexes = aBc->getNodeIndexes();
    auto vecR = aBc->R_Vector();
    for (size_t i = 0; i < 2; ++i) {
      result[indexes[i] - 1] += vecR[i];
    }
  }
  return result;
}

bool BoundaryConditions2D::isLinear() const { return m_Linear; }

/// FenestrationCommon::SquareMatrix< double > BoundaryConditions2D::DHMatrix()
/// const { 	FenestrationCommon::SquareMatrix< double > HDMatrix{
/// NodePool::Instance().maxIndex() }; 	for ( const auto & aBc : m_BCs ) { 		const
/// auto & bc = aBc.get(); 		auto indexes = bc.getNodeIndexes(); 		auto matDH =
/// bc.D_HMatrix(); 		for ( size_t i = 0; i < numOfBCNodes; ++i ) { 			for ( size_t j
/// = 0; j < numOfBCNodes; ++j ) { 				HDMatrix[ indexes[ i ] - 1 ][ indexes[ j ] -
/// 1 ] += matDH[ i ][ j ];
/// 			}
/// 		}
/// 	}
/// 	return HDMatrix;
/// }

void BoundaryConditions2D::updateNodeValues(const std::vector<double> &values,
                                            const Property property) {
  for (auto &aBc : m_BCs) {
    for (auto i = 0u; i < numOfBCNodes; ++i) {
      auto &node = aBc->getNode(i);
      auto index = node.getNodeNumber();
      node.setProperty(property, values[index - 1]);
    }
  }
}

} // namespace MoisThermFEM