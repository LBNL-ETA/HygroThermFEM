#include "Elements2D.hxx"
#include "NodePool.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	ElementsLinear2D::ElementsLinear2D(
			const std::vector< std::reference_wrapper< const IElementLinear2D > > & t_Elements )
			: m_Conductance( NodePool::Instance().maxIndex() ),
				m_Capacitance( NodePool::Instance().maxIndex() ) {

		// now integrate element matrices into global matrix
		for ( const IElementLinear2D & aElement : t_Elements ) {
			auto indexes = aElement.nodeIndexes();
			auto conductance = aElement.conductanceMatrix();
			auto capacitance = aElement.capacitanceMatrix();
			for ( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for ( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					m_Conductance[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += conductance[ i ][ j ];
					m_Capacitance[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += capacitance[ i ][ j ];
				}
			}
		}

	}

	SquareMatrix< double > & ElementsLinear2D::conductanceMatrix() {
		return m_Conductance;
	}

//	SquareMatrix< double > & ElementsLinear2D::thermalCapacitanceMatrix() {
//		return m_Capacitance;
//	}

	std::vector< double > ElementsLinear2D::getLumpedMass( const double DTime ) {
		auto size = m_Capacitance.size();

		std::vector< double > M( size );

		// Creates lump matrix
		for ( size_t i = 0; i < size; ++i ) {
			for ( size_t j = 0; j < size; ++j ) {
				M[ i ] += m_Capacitance[ i ][ j ];
			}
			M[ i ] /= DTime;
		}

		return M;
	}


}