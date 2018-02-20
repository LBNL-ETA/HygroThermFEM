#include <cmath>
#include <memory>

#include "FEMunique.hxx"
#include "IntegrationPoints.hxx"
#include "Node2D.hxx"

namespace MoisThermFEM {

	////////////////////////////////////////////////////////////////////////////
	///   IIntegrationPoints2D
	////////////////////////////////////////////////////////////////////////////

	std::vector< LocalPoint2D > IIntegrationPoints2D::getPoints() const {
		return m_Points;
	}

	std::vector< double > IIntegrationPoints2D::getWeights() const {
		return m_Weights;
	}

	size_t IIntegrationPoints2D::count() const {
		return m_Points.size();
	}

	////////////////////////////////////////////////////////////////////////////
	///   IIntegrationPoints1D
	////////////////////////////////////////////////////////////////////////////
	std::vector< LocalPoint1D > IIntegrationPoints1D::getPoints() const {
		return m_Points;
	}

	size_t IIntegrationPoints1D::count() const {
		return m_Points.size();
	}

	////////////////////////////////////////////////////////////////////////////
	//   SingleIntegrationPoint1D
	////////////////////////////////////////////////////////////////////////////
	SingleIntegrationPoint1D::SingleIntegrationPoint1D() {
		createPoints();
	}

	void SingleIntegrationPoint1D::createPoints() {
		m_Points.emplace_back( 0 );
	}

	////////////////////////////////////////////////////////////////////////////
	//   SingleIntegrationPoint2D
	////////////////////////////////////////////////////////////////////////////
	SingleIntegrationPoint2D::SingleIntegrationPoint2D() : IIntegrationPoints2D() {
		createPoints();
	}

	void SingleIntegrationPoint2D::createPoints() {
		m_Points.emplace_back( 0, 0 );
		m_Weights.push_back( 2 );
	}

	////////////////////////////////////////////////////////////////////////////
	//   TwoIntegrationPoint1D
	////////////////////////////////////////////////////////////////////////////
	TwoIntegrationPoint1D::TwoIntegrationPoint1D() {
		createPoints();
	}

	void TwoIntegrationPoint1D::createPoints() {
		m_Points.emplace_back( -1 / std::sqrt( 3 ) );
		m_Points.emplace_back( 1 / std::sqrt( 3 ) );
	}

	////////////////////////////////////////////////////////////////////////////
	//   TwoIntegrationPoint2D
	////////////////////////////////////////////////////////////////////////////
	TwoIntegrationPoint2D::TwoIntegrationPoint2D() : IIntegrationPoints2D() {
		createPoints();
	}

	void TwoIntegrationPoint2D::createPoints() {
		const std::vector< double > points = { -1 / std::sqrt( 3 ), 1 / std::sqrt( 3 ) };

		size_t coord1 = 0;
		size_t coord2 = 0;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( 1 );

		coord1 = 1;
		coord2 = 0;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( 1 );

		coord1 = 1;
		coord2 = 1;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( 1 );

		coord1 = 0;
		coord2 = 1;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( 1 );
	}

	////////////////////////////////////////////////////////////////////////////
	//   ThreeIntegrationPoint1D
	////////////////////////////////////////////////////////////////////////////
	ThreeIntegrationPoint1D::ThreeIntegrationPoint1D() {
		createPoints();
	}

	void ThreeIntegrationPoint1D::createPoints() {
		m_Points.emplace_back( -1 / std::sqrt( 3 ) );
		m_Points.emplace_back( 0 );
		m_Points.emplace_back( 1 / std::sqrt( 3 ) );
	}

	////////////////////////////////////////////////////////////////////////////
	//   ThreeIntegrationPoint2D
	////////////////////////////////////////////////////////////////////////////
	ThreeIntegrationPoint2D::ThreeIntegrationPoint2D() : IIntegrationPoints2D() {
		createPoints();
	}

	void ThreeIntegrationPoint2D::createPoints() {
		const std::vector< double > points = { -std::sqrt( 3.0 / 5.0 ), std::sqrt( 3.0 / 5.0 ), 0 };
		const std::vector< double > weights = { 8.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };

		size_t coord1 = 0;
		size_t coord2 = 0;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 1;
		coord2 = 0;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 1;
		coord2 = 1;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 0;
		coord2 = 1;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 2;
		coord2 = 0;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 1;
		coord2 = 2;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 2;
		coord2 = 1;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 0;
		coord2 = 2;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );

		coord1 = 2;
		coord2 = 2;
		m_Points.emplace_back( points[ coord1 ], points[ coord2 ] );
		m_Weights.push_back( weights[ coord1 ] * weights[ coord2 ] );
	}

	////////////////////////////////////////////////////////////////////////////
	//   IntegrationPoints2D
	////////////////////////////////////////////////////////////////////////////

	IntegrationPoints2D & IntegrationPoints2D::Instance() {
		// This should be thread safe in C++11
		static IntegrationPoints2D m_Instance;
		return m_Instance;
	}

	void IntegrationPoints2D::setIntegrationFormula( const IntegrationPointsFormula t_Formula ) {
		switch( t_Formula ) {
			case IntegrationPointsFormula::OnePoint:
				m_IntPoints2D = fem::make_unique< SingleIntegrationPoint2D >();
				m_IntPoints1D = fem::make_unique< SingleIntegrationPoint1D >();
				break;
			case IntegrationPointsFormula::TwoPoints:
				m_IntPoints2D = fem::make_unique< TwoIntegrationPoint2D >();
				m_IntPoints1D = fem::make_unique< TwoIntegrationPoint1D >();
				break;
			case IntegrationPointsFormula::ThreePoints:
				m_IntPoints2D = fem::make_unique< ThreeIntegrationPoint2D >();
				m_IntPoints1D = fem::make_unique< ThreeIntegrationPoint1D >();
				break;
		}
	}

	std::vector< LocalPoint2D > IntegrationPoints2D::getPoints2D() const {
		return m_IntPoints2D->getPoints();
	}

	std::vector< LocalPoint1D > IntegrationPoints2D::getPoints1D() const {
		return m_IntPoints1D->getPoints();
	}

	size_t IntegrationPoints2D::count2D() const {
		return m_IntPoints2D->count();
	}

	size_t IntegrationPoints2D::count1D() const {
		return m_IntPoints1D->count();
	}

	IntegrationPoints2D::IntegrationPoints2D() {
		// Two-point formula is default one
		m_IntPoints2D = fem::make_unique< TwoIntegrationPoint2D >();
		m_IntPoints1D = fem::make_unique< TwoIntegrationPoint1D >();
	}

	IntegrationPoints2D::~IntegrationPoints2D() = default;

}
