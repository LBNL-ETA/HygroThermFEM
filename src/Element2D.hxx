#pragma once

#include <memory>

#include "Functions.hxx"
#include "Material.hxx"
#include "Node2D.hxx"
#include "Quadrilateral2D.hxx"
#include "SquareMatrix.hxx"

namespace MoisThermFEM {

// Constant that holds number of nodes in certain elementsCreator
const std::size_t numOfQuadrilateralNodes = 4;
// const std::size_t numOfIntegrationPoints = 4;

//////////////////////////////////////////////////////////////////////////////
///  IQLEMatrix2D
//////////////////////////////////////////////////////////////////////////////

// Abstract class that forces users to perform matrix calculation in inherited
// version Depending on equations, matrices will have different calculation
// methods (for example capacitance and conductance matrices have different
// form)
class IQLEMatrix2D {
public:
  IQLEMatrix2D(const QuadrilateralLinearGlobal2D &t_Element);

  // Integrate matrix over all points of integration
  virtual FenestrationCommon::SquareMatrix
  integrate(const std::vector<double> &t_Values) const final;

protected:
  virtual void calculateMatrixInIntegrationPoint(
	const std::vector<double> &t_Values,
	const std::size_t t_IntegrationPointIndex,
	std::vector<std::vector<double>> & t_Matrix ) const final;

  const QuadrilateralLinearGlobal2D &m_Global2D;

  std::vector<FenestrationCommon::SquareMatrix> m_IntegrationMatrix;
};

//////////////////////////////////////////////////////////////////////////////
///  QLEConductance2D
//////////////////////////////////////////////////////////////////////////////

// Class to handle conductance matrix in global coordinate system
class QLEConductance2D : public IQLEMatrix2D {
public:
  QLEConductance2D(const QuadrilateralLinearGlobal2D &t_Element);
};

//////////////////////////////////////////////////////////////////////////////
///  QLEConductanceDerivative2D
//////////////////////////////////////////////////////////////////////////////

// Handles conductance part with derivative term
class QLEConductanceDerivative2D : public IQLEMatrix2D {
public:
  QLEConductanceDerivative2D(const QuadrilateralLinearGlobal2D &t_Element);

  // This updates integration matrix with new derivative values
  void updateIntegrationMatrix(const std::vector<double> &t_Values);

  void clearIntegrationMatrix();
};

//////////////////////////////////////////////////////////////////////////////
///  QLECapacitance2D
//////////////////////////////////////////////////////////////////////////////

// Class to handle capacitance matrix in global coordinate system
class QLECapacitance2D : public IQLEMatrix2D {
public:
  QLECapacitance2D(const QuadrilateralLinearGlobal2D &t_Element);
};

/// Keeping function pointers for QLEConductanceDerivative2D in Elements array
struct DerivativeFunction {
  DerivativeFunction( const iValue & fixedTerm, const iValue & derivativeTerm );

  iValue fixedTerm;
  iValue derivativeTerm;
};

//////////////////////////////////////////////////////////////////////////////
///  IElementLinear2D
//////////////////////////////////////////////////////////////////////////////

/// Class that handles creation of conductance and capacitance matrices in
/// linear 2D world. This class will be used by multiple governing equations
/// since basis of matrix creation are identical with only difference in what
/// coefficients are passed
class IElementLinear2D {
public:
  IElementLinear2D(const Node2D &t_Node1, const Node2D &t_Node2,
                   const Node2D &t_Node3, const Node2D &t_Node4,
                   const Material &t_Material);

  FenestrationCommon::SquareMatrix conductanceMatrix() const;

  FenestrationCommon::SquareMatrix conductanceDerivativeMatrix();

  FenestrationCommon::SquareMatrix capacitanceMatrix() const;

  Node2D &getNode(std::size_t index);

  bool haveBothNodes(const Node2D &t_Node1, const Node2D &t_Node2) const;

  std::vector<std::size_t> nodeIndexes() const;

  const Material &getMaterial() const;

protected:
  /// TODO: This did not work with reference_wrapper and it should. Check later.
  /// Reminder: Introduce pair of curve pointer and Property so that curve knows
  /// what to use
  std::vector<iValue> m_Conductance;
  std::vector<iValue> m_Capacitance;
  std::vector<DerivativeFunction> m_DerivativeConductance;

  const Material &m_Material;

private:

	template <typename T> class NodesVector {
	public:
		explicit NodesVector(const std::initializer_list<T> &__l)
			: vec_(__l), currentIndex(0), passedLast(false) {}

		typename std::vector<T>::const_iterator begin() const { return vec_.begin(); }

		typename std::vector<T>::const_iterator end() const { return vec_.end(); }

		std::size_t size() const { return vec_.size(); }

		T & operator[](const std::size_t index) { return vec_[index]; }
		const T & operator[](const std::size_t index) const { return vec_[index]; }

		/// Keeps iterating over unique elements of the vector
		T &current() { return vec_[currentIndex]; }

		bool last() { return passedLast; }

		T &previous() {

			auto validIndex = checkPrevIndex(currentIndex);

			while (vec_[validIndex] == vec_[currentIndex]) {
				validIndex = checkPrevIndex(validIndex);
			}

			return vec_[validIndex];
		}

		T &next() {
			auto validIndex = checkNextIndex(currentIndex);
			while (vec_[validIndex] == vec_[currentIndex]) {
				validIndex = checkNextIndex(validIndex);
			}

			return vec_[validIndex];
		}

		void moveToNext() {
			auto nextIndex = checkNextIndex(currentIndex);
			while (vec_[nextIndex] == vec_[currentIndex]) {
				nextIndex = checkNextIndex(nextIndex);
			}
			passedLast = nextIndex < currentIndex;
			currentIndex = nextIndex;
		}

	private:
		std::size_t checkNextIndex(const std::size_t index) const {
			auto validIndex = index;

			if (validIndex != vec_.size() - 1) {
				validIndex = index + 1;
			} else {
				validIndex = 0;
			}

			return validIndex;
		}

		std::size_t checkPrevIndex(const std::size_t index) const {
			auto validIndex = index;

			if (validIndex != 0) {
				validIndex = index - 1;
			} else {
				validIndex = vec_.size() - 1;
			}

			return validIndex;
		}
		std::vector<T> vec_;
		std::size_t currentIndex{0};
		bool passedLast{false};
	};

  NodesVector<Node2D> m_Node;

  QuadrilateralLinearGlobal2D m_Global2D;
  QLECapacitance2D m_QLECapacitance2D;
  QLEConductance2D m_QLEConductance2D;
  /// This one depends on functions and must be stored for every
  /// DerivativeConductance submatrix
  std::vector<QLEConductanceDerivative2D> m_QLEDerivativeConductance;
};

//////////////////////////////////////////////////////////////////////////////
///  ElementThermalLinear2D
//////////////////////////////////////////////////////////////////////////////

class ElementThermalLinear2D : public IElementLinear2D {
public:
  ElementThermalLinear2D(const Node2D &t_Node1, const Node2D &t_Node2,
                         const Node2D &t_Node3, const Node2D &t_Node4,
                         const Material &mat);
};

//////////////////////////////////////////////////////////////////////////////
///  ElementMoistureLinear2D
//////////////////////////////////////////////////////////////////////////////

class ElementMoistureLinear2D : public IElementLinear2D {
public:
  ElementMoistureLinear2D(const Node2D &t_Node1, const Node2D &t_Node2,
                          const Node2D &t_Node3, const Node2D &t_Node4,
                          const Material &mat);
};

} // namespace MoisThermFEM
