#include "Elements2D.hxx"
#include "NodePool.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	ElementsThermalLinear2D::ElementsThermalLinear2D( const std::vector< ElementThermalLinear2D > & t_Elements )
			: m_Conductivity( NodePool::Instance().maxIndex() ),
			  m_RhoCp( NodePool::Instance().maxIndex() ) {

		// now copy matrices into global matrix
		for( const auto & aElement : t_Elements ) {
			auto indexes = aElement.nodeIndexes();
			auto conductivity = aElement.conductivity();
			auto rhoCp = aElement.rhoCp();
			for( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					m_Conductivity[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += conductivity[ i ][ j ];
					m_RhoCp[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += rhoCp[ i ][ j ];
				}
			}
		}

	}

	SquareMatrix< double > & ElementsThermalLinear2D::thermalConductivity() {
		return m_Conductivity;
	}

	SquareMatrix< double > & ElementsThermalLinear2D::rhoCp() {
		return m_RhoCp;
	}

}