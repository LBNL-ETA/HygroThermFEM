#include "Elements2D.hxx"
#include "NodePool.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	ElementsLinear2D::ElementsLinear2D( const std::vector< ElementLinear2D > & t_Elements )
			: m_Conductance( NodePool::Instance().maxIndex() ),
			  m_Capacitance( NodePool::Instance().maxIndex() ) {

		// now integrate element matrices into global matrix
		for( const auto aElement : t_Elements ) {
			auto indexes = aElement.nodeIndexes();
			auto conductance = aElement.thermalConductanceMatrix();
			auto capacitance = aElement.thermalCapacitanceMatrix();
			for( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					m_Conductance[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += conductance[ i ][ j ];
					m_Capacitance[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += capacitance[ i ][ j ];
				}
			}
		}

	}

	SquareMatrix< double > & ElementsLinear2D::thermalConductanceMatrix() {
		return m_Conductance;
	}

	SquareMatrix< double > & ElementsLinear2D::thermalCapacitanceMatrix() {
		return m_Capacitance;
	}

}