#ifndef INTEGRATIONPOINTS_H_6dd73986_7133_4177_bd33_fb194b97a09d
#define INTEGRATIONPOINTS_H_6dd73986_7133_4177_bd33_fb194b97a09d

#include <vector>
#include <memory>

namespace MoisThermFEM {

  struct LocalPoint2D;
  struct LocalPoint1D;

	// Virtual class used in strategy pattern to define number of points for numerical integration
	// in one dimensional space
  ////////////////////////////////////////////////////////////////////////////
  //   IIntegrationPoints1D
  ////////////////////////////////////////////////////////////////////////////
  class IIntegrationPoints1D {
  public:
    virtual ~IIntegrationPoints1D() = default;
    IIntegrationPoints1D() = default;

    virtual std::vector< LocalPoint1D > getPoints() const final;
    virtual size_t count() const final;

  protected:
    virtual void createPoints() = 0;

    std::vector< LocalPoint1D > m_Points;
  };

  ////////////////////////////////////////////////////////////////////////////
  //   IIntegrationPoints2D
  ////////////////////////////////////////////////////////////////////////////

  // Virtual class used in strategy pattern to define number of points for numerical integration
	// in two dimensional space
  class IIntegrationPoints2D {
  public:
    virtual ~IIntegrationPoints2D() = default;
    IIntegrationPoints2D() = default;

    virtual std::vector< LocalPoint2D > getPoints() const final;
    virtual std::vector< double > getWeights() const final;
    virtual size_t count() const final;

  protected:
    virtual void createPoints() = 0;

    std::vector< LocalPoint2D > m_Points;
    std::vector< double > m_Weights;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   SingleIntegrationPoint1D
  ////////////////////////////////////////////////////////////////////////////

  class SingleIntegrationPoint1D : public IIntegrationPoints1D {
  public:
    SingleIntegrationPoint1D();

  private:
    void createPoints() override;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   SingleIntegrationPoint2D
  ////////////////////////////////////////////////////////////////////////////
  class SingleIntegrationPoint2D : public IIntegrationPoints2D {
  public:
    SingleIntegrationPoint2D();

  private:
    void createPoints() override;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   TwoIntegrationPoint1D
  ////////////////////////////////////////////////////////////////////////////
  class TwoIntegrationPoint1D : public IIntegrationPoints1D {
  public:
    TwoIntegrationPoint1D();

  private:
    void createPoints() override;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   TwoIntegrationPoint2D
  ////////////////////////////////////////////////////////////////////////////
  class TwoIntegrationPoint2D : public IIntegrationPoints2D {
  public:
    TwoIntegrationPoint2D();
    
  private:
    void createPoints() override;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   ThreeIntegrationPoint1D
  ////////////////////////////////////////////////////////////////////////////
  class ThreeIntegrationPoint1D : public IIntegrationPoints1D {
  public:
    ThreeIntegrationPoint1D();

  private:
    void createPoints() override;

  };

  ////////////////////////////////////////////////////////////////////////////
  //   ThreeIntegrationPoint2D
  ////////////////////////////////////////////////////////////////////////////
  class ThreeIntegrationPoint2D : public IIntegrationPoints2D {
  public:
    ThreeIntegrationPoint2D();
  
  private:
    void createPoints() override;
  
  };


  enum class IntegrationPointsFormula {
    OnePoint,
    TwoPoints,
    ThreePoints
  };

  ////////////////////////////////////////////////////////////////////////////
  //   IntegrationPoints2D
  ////////////////////////////////////////////////////////////////////////////

  // Interface class for integration points in 2D space
  class IntegrationPoints2D {
  public:
    static IntegrationPoints2D& Instance();
    void setIntegrationFormula( const IntegrationPointsFormula t_Formula );
    
    // Get integration points
    std::vector< LocalPoint2D > getPoints2D() const;
    std::vector< LocalPoint1D > getPoints1D() const;

    // Number of integration points
    size_t count2D() const;
    size_t count1D() const;

  private:
    IntegrationPoints2D();
    ~IntegrationPoints2D();

    // Integration points in local coordinate system
    std::unique_ptr< IIntegrationPoints2D > m_IntPoints2D;
    std::unique_ptr< IIntegrationPoints1D > m_IntPoints1D;

  };



}

#endif