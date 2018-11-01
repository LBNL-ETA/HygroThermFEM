#pragma once

#include <cstddef>
#include <vector>

#pragma warning(push, 0)
#include <Eigen/Core>
#include <Eigen/SparseCore>
#pragma warning(pop)

namespace FenestrationCommon {
// Works only with double
class SquareMatrix {
public:
  explicit SquareMatrix(const std::size_t m_size);
  // Constructor with triplets should be used for huge matrices. For example
  // 100,000 x 100,000 will not fit into memory of vector<vector<double>>
  explicit SquareMatrix(const std::size_t m_size,
                        const std::vector<Eigen::Triplet<double>> &tripletList);
  explicit SquareMatrix(
      const std::initializer_list<std::vector<double>> &tInput);
  explicit SquareMatrix(const std::vector<std::vector<double>> &tInput);
  explicit SquareMatrix(const std::vector<std::vector<double>> &&tInput);

  std::size_t size() const;
  void setZeros();
  void setIdentity();
  void setDiagonal(const std::vector<double> &tInput);
  SquareMatrix addDiagonal(const std::vector<double> &tInput) const;

  SquareMatrix inverse() const;

  double operator()(const std::size_t i, const std::size_t j) const;
  double &operator()(const std::size_t i, const std::size_t j);

  SquareMatrix mmultRows(const std::vector<double> &tInput) const;

  Eigen::SparseMatrix<double> getSparseMatrix() const;

  SquareMatrix operator*(const double value) const;
  friend SquareMatrix operator*(const double value, const SquareMatrix &other);
  friend SquareMatrix operator*(const SquareMatrix &first,
                                const SquareMatrix &second);
  SquareMatrix &operator*=(const SquareMatrix &other);
  SquareMatrix &operator+(const SquareMatrix &other);
  SquareMatrix &operator+=(const SquareMatrix &other);
  SquareMatrix &operator-(const SquareMatrix &other);
  SquareMatrix &operator-=(const SquareMatrix &other);

  std::vector<double> operator*(const std::vector<double> &tVec) const;
  friend std::vector<double> operator*(const std::vector<double> &first,
                                       const SquareMatrix &second);

  std::vector<std::vector<double>> toVector() const;

private:
  // explicit SquareMatrix(Eigen::MatrixXd && tMatrix);
  explicit SquareMatrix(Eigen::SparseMatrix<double> &&tMatrix);
  std::size_t m_size;
  Eigen::SparseMatrix<double> m_Matrix;
  // Eigen::MatrixXd m_Matrix;
};

} // namespace FenestrationCommon
