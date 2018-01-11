#pragma once

#include <vector>

namespace MoisThermFEM {

	struct Node2D;
	enum class Prop;

	class NodePool {
	public:
		static NodePool & Instance();

		Node2D & createNode( std::size_t const t_NodeNumber, double const t_x, double const t_y,
		                     double const t_temperature = 0 );

		Node2D & getNode( std::size_t const Index );

		std::size_t maxIndex() const;

		std::vector< double > nodeProperties( Prop t_Prop ) const;

		void clear();

	private:
		NodePool() = default;

		~NodePool() = default;

		std::vector< Node2D > m_Nodes;
	};

}