#pragma once

#include <Eigen/Sparse>
#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>
#include <vector>

namespace FenestrationCommon {
//////////////////////////////////////////////////////////////////////////////
/// Square matrices
//////////////////////////////////////////////////////////////////////////////
template <class T> class SparceSquareMatrix {
public:
  explicit SparceSquareMatrix(const std::size_t size)
      : m_size(size), m_Matrix(size, std::vector<double>(size, 0)) {}

  explicit SparceSquareMatrix(
      const std::initializer_list<std::vector<T>> &t_input)
      : m_size(t_input.size()),
        m_Matrix(m_size, std::vector<double>(m_size, 0)) {
    auto i = 0u;
    for (auto vec : t_input) {
      for (auto j = 0u; j < vec.size(); ++j) {
        m_Matrix[i][j] = vec[j];
      }
      ++i;
    }
  }

  explicit SparceSquareMatrix(const std::vector<std::vector<T>> &t_input)
      : m_size(t_input.size()),
        m_Matrix(m_size, std::vector<double>(m_size, 0)) {
    auto i = 0u;
    for (auto vec : t_input) {
      for (auto j = 0u; j < vec.size(); ++j) {
        m_Matrix[i][j] = vec[j];
      }
      ++i;
    }
  }

  explicit SparceSquareMatrix(const std::vector<std::vector<T>> &&t_input)
      : m_size(t_input.size()) {
    auto i = 0u;
    for (auto vec : t_input) {
      for (auto j = 0u; j < vec.size(); ++j) {
        m_Matrix[i][j] = vec[j];
      }
      ++i;
    }
  }

  std::size_t size() const { return m_size; }

  std::vector<Eigen::Triplet<double>> triplets() const {
    std::vector<Eigen::Triplet<double>> tri;
    for (auto i = 0; i < m_size; ++i) {
      for (auto j = 0; j < m_size; ++j) {
        if (m_Matrix[i][j] != 0) {
          tri.emplace_back(i, j, m_Matrix[i][j]);
        }
      }
    }
    return tri;
  }

  void setZeros() { m_Matrix.clear(); }

  void setIdentity() {
    setZeros();
    for (unsigned i = 0; i < m_size; ++i) {
      m_Matrix[i][i] = 1;
    }
  }

  void setDiagonal(const std::vector<T> &t_Values) {
    if (t_Values.size() != m_size) {
      throw std::runtime_error("Supplied vector size mismatch matrix size");
    }
    setZeros();
    for (unsigned i = 0; i < m_size; ++i) {
      m_Matrix[i][i] = t_Values[i];
    }
  }

  /// Makes equivalent upper triangular matrix. It also changes vector so that
  /// it can be used in solver
  void makeUpperTriangular(std::vector<T> &rightSide) {
    if (m_size != rightSide.size()) {
      throw std::runtime_error("Provided vector is not same size as matrix.");
    }

    std::vector<std::set<size_t>> RowIndex(m_size, std::set<size_t>());

    /// build indexer for fast sparse solution
    for (auto i = 0u; i < m_size; ++i) {
      for (auto j = 0u; j < m_size; ++j) {
        if (m_Matrix[i][j] != 0) {
          RowIndex[i].insert(j);
        }
      }
    }

    for (auto i = 0u; i < m_size - 1; ++i) {
      for (auto j = i + 1; j < m_size; ++j)
      // for(auto & j : m_ColumnIndex[i])
      {
        if (m_Matrix[j][i] != 0)
        // if(equation >= i + 1 )
        {
          const auto coeff = -m_Matrix[i][i] / m_Matrix[j][i];
          rightSide[j] = coeff * rightSide[j] + rightSide[i];

          // Need to iterate over union of two rows because we should sum all
          // non-zeros.
          std::vector<size_t> unionRange;
          std::set_union(RowIndex[i].begin(), RowIndex[i].end(),
                         RowIndex[j].begin(), RowIndex[j].end(),
                         std::inserter(unionRange, unionRange.begin()));

          // for(auto k = 0u; k < m_size; ++k)
          for (auto &row : unionRange) {
            // m_Matrix[j][k] = coeff * m_Matrix[j][k] + m_Matrix[i][k];
            m_Matrix[j][row] = coeff * m_Matrix[j][row] + m_Matrix[i][row];
          }
        }
      }
    }
  };

  SparceSquareMatrix<T> LU() const {
    auto D{SparceSquareMatrix<T>(m_size)};
    D = *this;

    for (auto k = 0u; k <= m_size - 2; ++k) {
      for (auto j = k + 1; j <= m_size - 1; ++j) {
        auto x = D(j, k) / D(k, k);
        for (auto i = k; i <= m_size - 1; ++i) {
          D(j, i) = D(j, i) - x * D(k, i);
        }
        D(j, k) = x;
      }
    }

    return D;
  }

  SparceSquareMatrix<T> inverse() const {
    // return LU decomposed matrix of current matrix
    auto aLU = LU();

    // find the inverse
    auto inverse = SparceSquareMatrix<T>(m_size);
    std::vector<double> d(m_size);
    std::vector<double> y(m_size);

    const auto size = int(m_size - 1);

    for (auto m = 0; m <= size; ++m) {
      fill(d.begin(), d.end(), 0);
      fill(y.begin(), y.end(), 0);
      d[m] = 1;
      for (auto i = 0; i <= size; ++i) {
        double x = 0;
        for (auto j = 0; j <= i - 1; ++j) {
          x = x + aLU(i, j) * y[j];
        }
        y[i] = (d[i] - x);
      }

      for (auto i = size; i >= 0; --i) {
        auto x = 0.0;
        for (auto j = i + 1; j <= size; ++j) {
          x = x + aLU(i, j) * inverse(j, m);
        }
        inverse(i, m) = (y[i] - x) / aLU(i, i);
      }
    }

    return inverse;
  }

  SparceSquareMatrix<T> mmultRows(const std::vector<T> &t_vector) const {
    SparceSquareMatrix<T> result{t_vector.size()};

    for (auto i = 0u; i < t_vector.size(); ++i) {
      for (auto j = 0u; j < t_vector.size(); ++j) {
        result(i, j) = m_Matrix[i][j] * t_vector[i];
      }
    }
    return result;
  }

  SparceSquareMatrix<T> addDiagonal(const std::vector<T> &t_Vector) {
    if (m_size != t_Vector.size()) {
      std::runtime_error("Matrix and vector have different sizes.");
    }

    SparceSquareMatrix<T> aMatrix{m_size};
    for (auto i = 0u; i < m_size; ++i) {
      for (auto j = 0u; j < m_size; ++j) {
        aMatrix(i, j) = m_Matrix[i][j];
      }
      aMatrix(i, i) += t_Vector[i];
    }

    return aMatrix;
  }

  T operator()(const size_t i, const size_t j) const { return m_Matrix[i][j]; }

  T &operator()(const size_t i, const size_t j) { return m_Matrix[i][j]; }

  SparceSquareMatrix<T> &operator+(const SparceSquareMatrix<T> &rhs) {
    if (rhs.m_size != m_size) {
      throw std::runtime_error("Matrices must be identical in size.");
    }

    for (unsigned i = 0; i < m_size; ++i) {
      for (unsigned j = 0; j < m_size; ++j) {
        m_Matrix[i][j] = m_Matrix[i][j] + rhs.m_Matrix[i][j];
      }
    }

    return *this;
  }

  SparceSquareMatrix<T> &operator+=(const SparceSquareMatrix<T> &rhs) {
    return operator+(rhs);
  }

  SparceSquareMatrix<T> &operator-(const SparceSquareMatrix<T> &rhs) {
    if (rhs.m_size != m_size) {
      throw std::runtime_error("Matrices must be identical in size.");
    }

    for (unsigned i = 0; i < m_size; ++i) {
      for (unsigned j = 0; j < m_size; ++j) {
        m_Matrix[i][j] = m_Matrix[i][j] - rhs.m_Matrix[i][j];
      }
    }

    return *this;
  }

  SparceSquareMatrix<T> operator*(const SparceSquareMatrix<T> &rhs) {
    if (m_size != rhs.size()) {
      throw std::runtime_error("Matrices must be identical in size.");
    }

    SparceSquareMatrix<T> temp(m_size);

    for (auto i = 0u; i < m_size; ++i) {
      for (auto k = 0u; k < m_size; ++k) {
        for (auto j = 0u; j < m_size; ++j) {
          temp(i, j) = temp(i, j) + m_Matrix[i][k] * rhs.m_Matrix[k][j];
        }
      }
    }

    return temp;
  }

  friend SparceSquareMatrix<T> operator*(const SparceSquareMatrix<T> &lhs,
                                         const double t_Value) {
    SparceSquareMatrix aMatrix{lhs.size()};

    for (auto i = 0u; i < lhs.size(); ++i) {
      for (auto k = 0u; k < lhs.size(); ++k) {
        aMatrix(i, k) += t_Value * lhs.m_Matrix[i][k];
      }
    }

    return aMatrix;
  }

  friend SparceSquareMatrix<T> operator*(const double t_Value,
                                         const SparceSquareMatrix<T> &lhs) {
    return lhs * t_Value;
  }

  friend std::vector<T> operator*(const std::vector<T> &t_vector,
                                  const SparceSquareMatrix &t_matrix) {
    if (t_vector.size() != t_matrix.size()) {
      throw std::runtime_error("Vector and matrix have incompatible sizes.");
    }

    std::vector<T> result(t_vector.size());
    for (auto i = 0u; i < t_vector.size(); ++i) {
      for (auto j = 0u; j < t_vector.size(); ++j) {
        result[i] += t_matrix.m_Matrix[j][i] * t_vector[j];
      }
    }

    return result;
  }

  friend std::vector<T> operator*(const SparceSquareMatrix &t_matrix,
                                  const std::vector<T> &t_vector) {
    if (t_vector.size() != t_matrix.size()) {
      throw std::runtime_error("Vector and matrix have incompatible sizes.");
    }

    std::vector<T> result(t_vector.size());
    for (auto i = 0u; i < t_vector.size(); ++i) {
      for (auto j = 0u; j < t_vector.size(); ++j) {
        result[i] += t_matrix.m_Matrix[i][j] * t_vector[j];
      }
    }

    return result;
  }

protected:
  size_t m_size;
  std::vector<std::vector<double>> m_Matrix;
};

} // namespace FenestrationCommon
