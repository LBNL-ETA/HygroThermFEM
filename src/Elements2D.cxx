#include "Elements2D.hxx"
#include "NodePool.hxx"

using namespace FenestrationCommon;

namespace MoisThermFEM {

	SparceSquareMatrix< double > ElementsLinear2D::conductanceMatrix() {
		SparceSquareMatrix< double > result{ NodePool::Instance().maxIndex() };
		// now integrate element matrices into global matrix
		for ( auto & aElement : m_Elements ) {
			auto indexes = aElement->nodeIndexes();
			auto conductance = aElement->conductanceMatrix();
			auto condDer = aElement->conductanceDerivativeMatrix();
			for ( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for ( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					result( indexes[ i ] - 1, indexes[ j ] - 1 ) +=
							( conductance( i, j ) + condDer( i, j) );
				}
			}
		}
		return result;
	}

//	SquareMatrix< double > & ElementsLinear2D::thermalCapacitanceMatrix() {
//		return m_Capacitance;
//	}

	std::vector< double > ElementsLinear2D::getLumpedMass( const double DTime ) {
		FenestrationCommon::SparceSquareMatrix< double > Capacitance{ NodePool::Instance().maxIndex() };

		// now integrate element matrices into global matrix
		for ( auto & aElement : m_Elements ) {
			auto indexes = aElement->nodeIndexes();
			auto capacitance = aElement->capacitanceMatrix();
			for ( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for ( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					Capacitance( indexes[ i ] - 1, indexes[ j ] - 1 ) += capacitance( i, j );
				}
			}
		}

		auto size = Capacitance.size();

		std::vector< double > M( size, 0 );

		// Creates lump matrix
		for ( size_t i = 0; i < size; ++i ) {
			for ( size_t j = 0; j < size; ++j ) {
				M[ i ] += Capacitance(i, j );
			}
			M[ i ] /= DTime;
		}

		return M;
	}

	FenestrationCommon::SparceSquareMatrix< double > ElementsLinear2D::getMassMatrix( const double DTime ) {
		FenestrationCommon::SparceSquareMatrix< double > Capacitance{ NodePool::Instance().maxIndex() };

		// now integrate element matrices into global matrix
		for ( auto & aElement : m_Elements ) {
			auto indexes = aElement->nodeIndexes();
			auto capacitance = aElement->capacitanceMatrix();
			for ( size_t i = 0; i < numOfQuadrilateralNodes; ++i ) {
				for ( size_t j = 0; j < numOfQuadrilateralNodes; ++j ) {
					Capacitance( indexes[ i ] - 1, indexes[ j ] - 1 ) += capacitance(i, j ) / DTime;
				}
			}
		}

		return Capacitance;
	}

	ElementsLinear2D::ElementsLinear2D() : m_Linear( true ) {

	}

	bool ElementsLinear2D::isLinear() const {
		return m_Linear;
	}

	void ElementsLinear2D::updateNodeValues( const std::vector< double > & values,
																					 const Property property ) {
		for ( auto & aBc : m_Elements ) {
			for ( auto i = 0u; i < numOfQuadrilateralNodes; ++i ) {
				auto & node = aBc->getNode( i );
				auto index = node.getNodeNumber();
				node.setProperty( property, values[ index - 1 ] );
			}

		}
	}

	IElementLinear2D *
	ElementsLinear2D::findElement( const Node2D & t_Node1, const Node2D & t_Node2 ) {
		IElementLinear2D * el = nullptr;
		for( auto & element : m_Elements ) {
			if( element->haveBothNodes( t_Node1, t_Node2 ) ) {
				el = element.get();
			}
		}
		return el;
	}


}