#include <assert.h>

#include "Node2D.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////////////////////////
	//   LocalPoint1D
	////////////////////////////////////////////////////////////////////////////

	LocalPoint1D::LocalPoint1D( double const t_ksi )
			: ksi( t_ksi ) {

	}

	LocalPoint1D::LocalPoint1D( LocalPoint1D const & t_LocalPoint ) {
		ksi = t_LocalPoint.ksi;
	}

	////////////////////////////////////////////////////////////////////////////
	//   LocalPoint2D
	////////////////////////////////////////////////////////////////////////////

	LocalPoint2D::LocalPoint2D( double const t_ksi, double const t_eta )
			: ksi( t_ksi ), eta( t_eta ) {
	}

	LocalPoint2D::LocalPoint2D( LocalPoint2D const & t_LocalPoint ) {
		ksi = t_LocalPoint.ksi;
		eta = t_LocalPoint.eta;
	}

	////////////////////////////////////////////////////////////////////////////
	//  Properties
	////////////////////////////////////////////////////////////////////////////

	Property::Property( const double t_Temperature, const double t_Humidity,
	                    const double t_Pressure ) {
		m_Property[ Prop::temperature ] = t_Temperature;
		m_Property[ Prop::humidity ] = t_Humidity;
		m_Property[ Prop::pressure ] = t_Pressure;
	}

	double Property::getValue( Prop t_Property ) const {
		return m_Property.at( t_Property );
	}

	Property::Property( const Property & other ) : m_Property( other.m_Property ) {

	}

	Property & Property::operator=( const Property & other ) {
		m_Property = other.m_Property;
		return *this;
	}

	void Property::setValue( Prop t_Property, double t_Value ) {
		m_Property[ t_Property ] = t_Value;
	}

	////////////////////////////////////////////////////////////////////////////
	//   Node2D
	////////////////////////////////////////////////////////////////////////////

	Node2D::Node2D( const std::size_t t_NodeNumber, const double t_x, const double t_y,
	                const Property & t_Property )
			: m_NodeNumber( t_NodeNumber ), m_x( t_x ), m_y( t_y ), m_Property( t_Property ) {

	}

	Node2D::Node2D( const Node2D & t_Node )
			: m_NodeNumber( t_Node.m_NodeNumber ), m_x( t_Node.m_x ), m_y( t_Node.m_y ),
			  m_Property( t_Node.m_Property ) {

	}

	Node2D & Node2D::operator=( const Node2D & other ) {
		m_NodeNumber = other.m_NodeNumber;
		m_x = other.m_x;
		m_y = other.m_y;
		m_Property = other.m_Property;

		return *this;
	}

	size_t Node2D::getNodeNumber() const {
		return m_NodeNumber;
	}

	double Node2D::X() const {
		return m_x;
	}

	double Node2D::Y() const {
		return m_y;
	}

	double Node2D::getProperty( Prop t_Property ) const {
		return m_Property.getValue( t_Property );
	}

	void Node2D::setProperty( Prop t_Property, double t_value ) {
		m_Property.setValue( t_Property, t_value );
	}

	////////////////////////////////////////////////////////////////////////////
	//   INodesStorage
	////////////////////////////////////////////////////////////////////////////

	INodesStorage::INodesStorage( std::initializer_list< Node2D > t_Nodes )
			: m_Nodes( t_Nodes ) {
	}

	Node2D INodesStorage::getNode( const std::size_t Index ) const {
		assert( Index < m_Nodes.size() );
		return m_Nodes[ Index ];
	}

	std::vector< std::size_t > INodesStorage::getNodeIndexes() const {
		std::vector< std::size_t > indexes;
		for( const auto & aNode : m_Nodes ) {
			indexes.push_back( aNode.getNodeNumber() );
		}
		return indexes;
	}

	////////////////////////////////////////////////////////////////////////////
	//   LineNodes2D
	////////////////////////////////////////////////////////////////////////////

	LineNodes2D::LineNodes2D( const Node2D & t_Node1, const Node2D & t_Node2 )
			: INodesStorage( { t_Node1, t_Node2 } ) {

	}

	////////////////////////////////////////////////////////////////////////////
	//   QuadrilateralNodes2D
	////////////////////////////////////////////////////////////////////////////

	QuadrilateralNodes2D::QuadrilateralNodes2D( const Node2D & t_Node1, const Node2D & t_Node2,
	                                            const Node2D & t_Node3, const Node2D & t_Node4 )
			: INodesStorage( { t_Node1, t_Node2, t_Node3, t_Node4 } ) {

	}
}