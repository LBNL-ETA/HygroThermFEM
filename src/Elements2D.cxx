#include "Elements2D.hxx"
#include "NodePool.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	SquareMatrix< double > ElementsLinear2D::conductanceMatrix() {
		SquareMatrix< double > result{ NodePool::Instance().maxIndex() };
		// now integrate element matrices into global matrix
		for ( auto & aElement : m_Elements ) {
			auto indexes = aElement->nodeIndexes();
			auto conductance = aElement->conductanceMatrix();
			for ( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for ( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					result[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += conductance[ i ][ j ];
				}
			}
		}
		return result;
	}

//	SquareMatrix< double > & ElementsLinear2D::thermalCapacitanceMatrix() {
//		return m_Capacitance;
//	}

	Vector< double > ElementsLinear2D::getLumpedMass( const double DTime ) {
		FenestrationCommon::SquareMatrix< double > Capacitance { NodePool::Instance().maxIndex() };

		// now integrate element matrices into global matrix
		for ( auto & aElement : m_Elements ) {
			auto indexes = aElement->nodeIndexes();
			auto capacitance = aElement->capacitanceMatrix();
			for ( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for ( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					Capacitance[ indexes[ i ] - 1 ][ indexes[ j ] - 1 ] += capacitance[ i ][ j ];
				}
			}
		}

		auto size = Capacitance.size();

		FenestrationCommon::Vector< double > M( size, 0 );

		// Creates lump matrix
		for ( size_t i = 0; i < size; ++i ) {
			for ( size_t j = 0; j < size; ++j ) {
				M[ i ] += Capacitance[ i ][ j ];
			}
			M[ i ] /= DTime;
		}

		return M;
	}

	ElementsLinear2D::ElementsLinear2D() : m_Linear( true ) {

	}

	bool ElementsLinear2D::isLinear() const {
		return m_Linear;
	}

	void ElementsLinear2D::updateNodeValues( const std::vector< double > & values,
																					 const Property property ) {
		for ( auto & aBc : m_Elements ) {
			for( auto i = 0u; i < numOfQuadrilateralNodes; ++i ) {
				auto & node = aBc->getNode( i );
				auto index = node.getNodeNumber();
				node.setProperty( property, values[ index - 1 ] );
			}

		}
	}


}