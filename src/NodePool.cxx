#include <assert.h>
#include <algorithm>

#include "NodePool.hxx"
#include "Node2D.hxx"

namespace MoisThermFEM {

	NodePool & NodePool::Instance() {
		static NodePool m_Instance;
		return m_Instance;
	}

	Node2D & NodePool::createNode( const std::size_t t_NodeNumber, const double t_x, const double t_y,
	                               const double t_temperature ) {
		auto aNode = Node2D( t_NodeNumber, t_x, t_y, t_temperature );
		m_Nodes.push_back( aNode );
		return m_Nodes.back();
	}

	Node2D & NodePool::getNode( const size_t Index ) {
		assert( Index < m_Nodes.size() );
		return m_Nodes[ Index ];
	}

	size_t NodePool::maxIndex() const {
		Node2D aNode = *max_element( m_Nodes.begin(), m_Nodes.end(),
				[ ]( const Node2D & a, const Node2D & b ) { return a.getNodeNumber() < b.getNodeNumber();
				} );
		return aNode.getNodeNumber();
	}

	std::vector< double > NodePool::nodeProperties( Prop t_Prop ) const {
		std::vector< double > aVector;
		for( const Node2D & aNode : m_Nodes ) {
			aVector.push_back( aNode.getProperty( t_Prop ) );
		}
		return aVector;
	}

	void NodePool::clear() {
		m_Nodes.clear();
	}

}