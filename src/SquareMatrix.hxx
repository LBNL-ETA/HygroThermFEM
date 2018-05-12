#pragma once

#include <cstddef>
#include <vector>
#pragma warning(push,0)
#include <Eigen/Core>
#pragma warning(pop)

namespace FenestrationCommon
{
	class CLinearSolver;

    // Works only with double
    class SquareMatrix
    {
    	friend class CLinearSolver;
    public:
        explicit SquareMatrix(const std::size_t m_size);
        explicit SquareMatrix(const std::initializer_list<std::vector<double>> & tInput);
        explicit SquareMatrix(const std::vector<std::vector<double>> & tInput);
        explicit SquareMatrix(const std::vector<std::vector<double>> && tInput);

        std::size_t size() const;
        void setZeros();
        void setIdentity();
        void setDiagonal(const std::vector<double> & tInput);

        SquareMatrix inverse() const;

        double operator()(const std::size_t i, const std::size_t j) const;
        double & operator()(const std::size_t i, const std::size_t j);

        SquareMatrix mmultRows(const std::vector<double> & tInput) const;
        SquareMatrix addDiagonal(const std::vector<double> & tInput) const;

        friend SquareMatrix operator*(const SquareMatrix & first, const SquareMatrix & second);
        SquareMatrix & operator*=(const SquareMatrix & other);
        SquareMatrix & operator+(const SquareMatrix & other);
        SquareMatrix & operator+=(const SquareMatrix & other);
        SquareMatrix & operator-(const SquareMatrix & other);
        SquareMatrix & operator-=(const SquareMatrix & other);

        std::vector<double> operator*(const std::vector<double> & tVec) const;
		friend std::vector<double> operator*(const std::vector<double> & first, const SquareMatrix & second);
		friend SquareMatrix operator*(const double value, const SquareMatrix & matrix);
		friend SquareMatrix operator*(const SquareMatrix & matrix, const double value);


    private:
        explicit SquareMatrix(Eigen::MatrixXd && tMatrix);
		explicit SquareMatrix(Eigen::MatrixXd & tMatrix);
        std::size_t m_size;
        Eigen::MatrixXd m_Matrix;
    };

}   // namespace FenestrationCommon