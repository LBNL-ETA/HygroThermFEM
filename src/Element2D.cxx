#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "QuadrilateralLocal2D.hxx"

namespace Conrad {

	ElementLinear2D::ElementLinear2D(
		Node2D const & t_Node1,
		Node2D const & t_Node2,
		Node2D const & t_Node3,
		Node2D const & t_Node4,
		double const t_Cond,
		double const t_Rho,
		double const t_Cp ) :
		m_Element( t_Node1, t_Node2, t_Node3, t_Node4 ),
		m_ElementNodes( t_Node1, t_Node2, t_Node3, t_Node4 ),
		m_Cond( t_Cond ), m_Rho( t_Rho ), m_Cp( t_Cp )
	{
		m_Conductivity.resize( 4 );
		for ( auto & item : m_Conductivity ) {
			item.resize( 4 );
		}

		m_RhoCp.resize( 4 );
		for ( auto & item : m_RhoCp ) {
			item.resize( 4 );
		}

		const auto count = IntegrationPoints2D::Instance().count2D();

		for ( unsigned i = 0; i < count; ++i ) {
			auto ConductionInIntegraionPoint = calculateConductionMatrix( i );
			auto RhoCPInIntegrationPoint = calculateRhoCpMatrix( i );
			for ( unsigned j = 0; j < 4; ++j ) {
				for ( unsigned k = 0; k < 4; ++k ) {
					m_Conductivity[ j ][ k ] += ConductionInIntegraionPoint[ j ][ k ];
					m_RhoCp[ j ][ k ] += RhoCPInIntegrationPoint[ j ][ k ];
				}
			}
		}
	}

	ElementLinear2D::ElementLinear2D( ElementLinear2D const & t_Element ) :
		m_Element( t_Element.m_Element ), m_ElementNodes( t_Element.m_ElementNodes ),
        m_Cond( t_Element.m_Cond )
	{
		m_Conductivity.resize( t_Element.m_Conductivity.size() );
		for ( size_t i = 0; i < m_Conductivity.size(); ++i ) {
			for ( auto const aValue : t_Element.m_Conductivity[ i ] ) {
				m_Conductivity[ i ].push_back( aValue );
			}
		}

		m_RhoCp.resize( t_Element.m_RhoCp.size() );
		for ( size_t i = 0; i < m_RhoCp.size(); ++i ) {
			for ( auto const aValue : t_Element.m_RhoCp[ i ] ) {
				m_RhoCp[ i ].push_back( aValue );
			}
		}
	}

	std::vector< size_t >
	ElementLinear2D::nodeIndexes() const
	{
		return m_ElementNodes.getNodeIndexes();
	}

	std::vector< std::vector< double > >
	ElementLinear2D::conductivity() const
	{
		return m_Conductivity;
	}

	std::vector< std::vector< double > >
	ElementLinear2D::RhoCp() const
	{
		return m_RhoCp;
	}

	std::vector< std::vector< double > > ElementLinear2D::calculateConductionMatrix(
		const size_t t_IntegrationPointIndex ) const
	{
		std::vector< std::vector< double > > aMatrix( 4 );
		for ( auto & item : aMatrix ) {
			item.resize( 4 );
		}

		auto DPsiDx = m_Element.DPsiDx( t_IntegrationPointIndex );
		auto DPsiDy = m_Element.DPsiDy( t_IntegrationPointIndex );
		auto det = m_Element.det( t_IntegrationPointIndex );

		for ( size_t i = 0; i < 4; ++i ) {
			for ( size_t j = 0; j < 4; ++j ) {
				aMatrix[ i ][ j ] = ( DPsiDx[ i ] * DPsiDx[ j ] + DPsiDy[ i ] * DPsiDy[ j ] ) * det * m_Cond;
			}
		}

		return aMatrix;
	}

	std::vector< std::vector< double > > ElementLinear2D::calculateRhoCpMatrix(
		size_t const t_IntegrationPointIndex ) const
	{
		std::vector< std::vector< double > > aMatrix( 4 );
		for ( auto & item : aMatrix ) {
			item.resize( 4 );
		}

		auto & aElement = QuadrilateralLinearLocal2D::Instance();

		auto psi = aElement.VPsi( t_IntegrationPointIndex );
		const auto det = m_Element.det( t_IntegrationPointIndex );

		for ( size_t i = 0; i < 4; ++i ) {
			for ( size_t j = 0; j < 4; ++j ) {
				aMatrix[ i ][ j ] = psi[ i ] * psi[ j ] * det * m_Rho * m_Cp;
			}
		}

		return aMatrix;
	}

}
