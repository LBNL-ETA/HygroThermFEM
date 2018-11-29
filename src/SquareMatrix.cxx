#include "SquareMatrix.hxx"

#pragma warning(push, 0)
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#pragma warning(pop)

namespace FenestrationCommon
{
    SquareMatrix::SquareMatrix(const std::size_t size) : m_size(size), m_Matrix(size, size)
    {}

    SquareMatrix::SquareMatrix(const std::size_t tSize,
                               const std::vector<Eigen::Triplet<double>> & tripletList) :
        m_size(tSize),
        m_Matrix(m_size, m_size)
    {
        m_Matrix.setFromTriplets(tripletList.begin(), tripletList.end());
        m_size = m_Matrix.innerSize();
    }

    SquareMatrix::SquareMatrix(const std::initializer_list<std::vector<double>> & tInput) :
        m_size(tInput.size()),
        m_Matrix(m_size, m_size)
    {
        // std::vector<Eigen::Triplet<double>> tripletList;
        auto i = 0u;
        for(auto vec : tInput)
        {
            for(auto j = 0u; j < vec.size(); ++j)
            {
                // tripletList.emplace_back(i, j, vec[j]);
                m_Matrix.coeffRef(i, j) = vec[j];
            }
            ++i;
        }
        // m_Matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    SquareMatrix::SquareMatrix(const std::vector<std::vector<double>> & tInput) :
        m_size(tInput.size()),
        m_Matrix(m_size, m_size)
    {
        // std::vector<Eigen::Triplet<double>> tripletList;
        for(auto i = 0u; i < tInput.size(); ++i)
        {
            for(auto j = 0u; j < tInput.size(); ++j)
            {
                if(tInput[i][j] != 0)
                {
                    // tripletList.emplace_back(i, j, tInput[i][j]);
                    m_Matrix.coeffRef(i, j) = tInput[i][j];
                }
            }
        }
        // m_Matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    SquareMatrix::SquareMatrix(const std::vector<std::vector<double>> && tInput) :
        m_size(tInput.size()),
        m_Matrix(m_size, m_size)
    {
        // std::vector<Eigen::Triplet<double>> tripletList;
        for(auto i = 0u; i < tInput.size(); ++i)
        {
            for(auto j = 0u; j < tInput.size(); ++j)
            {
                if(tInput[i][j] != 0)
                {
                    // tripletList.emplace_back(i, j, tInput[i][j]);
                    m_Matrix.coeffRef(i, j) = tInput[i][j];
                }
            }
        }
        // m_Matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    std::size_t SquareMatrix::size() const
    {
        return m_size;
    }

    void SquareMatrix::setZeros()
    {
        m_Matrix.setZero();
    }

    void SquareMatrix::setIdentity()
    {
        m_Matrix.setIdentity();
    }

    void SquareMatrix::setDiagonal(const std::vector<double> & tInput)
    {
        if(tInput.size() != m_size)
        {
            throw std::runtime_error("Matrix and vector must be same size.");
        }
        m_Matrix.setZero();

        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(m_size);
        for(auto i = 0u; i < m_size; ++i)
        {
            tripletList.emplace_back(i, i, tInput[i]);
        }
        m_Matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    SquareMatrix SquareMatrix::addDiagonal(const std::vector<double> & tInput) const
    {
        SquareMatrix mat{*this};
        SquareMatrix diag{m_size};
        diag.setDiagonal(tInput);
        return mat + diag;
    }

    SquareMatrix::SquareMatrix(Eigen::SparseMatrix<double> && tMatrix) :
        m_size(tMatrix.innerSize()),
        m_Matrix(tMatrix)
    {}

    SquareMatrix SquareMatrix::inverse() const
    {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(m_Matrix);
        Eigen::SparseMatrix<double> I(m_size, m_size);
        I.setIdentity();
        const auto A_inv = solver.solve(I);
        return SquareMatrix{A_inv};
    }

    double SquareMatrix::operator()(const std::size_t i, const std::size_t j) const
    {
        return m_Matrix.coeff(i, j);
    }

    double & SquareMatrix::operator()(const std::size_t i, const std::size_t j)
    {
        return m_Matrix.coeffRef(i, j);
    }

    SquareMatrix & SquareMatrix::operator*=(const SquareMatrix & other)
    {
        *this = *this * other;
        return *this;
    }

    SquareMatrix & SquareMatrix::operator+(const SquareMatrix & other)
    {
        m_Matrix = m_Matrix + other.m_Matrix;
        return *this;
    }

    SquareMatrix & SquareMatrix::operator+=(const SquareMatrix & other)
    {
        return operator+(other);
    }

    SquareMatrix & SquareMatrix::operator-(const SquareMatrix & other)
    {
        m_Matrix = m_Matrix - other.m_Matrix;
        return *this;
    }

    SquareMatrix & SquareMatrix::operator-=(const SquareMatrix & other)
    {
        return operator-(other);
    }

    SquareMatrix SquareMatrix::mmultRows(const std::vector<double> & tInput) const
    {
        Eigen::VectorXd vec = Eigen::VectorXd::Map(tInput.data(), tInput.size());
        SquareMatrix res{m_Matrix * vec.asDiagonal()};
        return res;
    }

    SquareMatrix SquareMatrix::operator*(const double value) const
    {
        SquareMatrix res{m_Matrix * value};
        return res;
    }

    std::vector<double> SquareMatrix::operator*(const std::vector<double> & tVec) const
    {
        const Eigen::VectorXd vec = Eigen::VectorXd::Map(tVec.data(), tVec.size());
        Eigen::VectorXd res = m_Matrix * vec;
        return std::vector<double>(res.data(), res.data() + res.rows() * res.cols());
    }

    SquareMatrix operator*(const double value, const SquareMatrix & other)
    {
        auto res{other * value};
        return res;
    }

    SquareMatrix operator*(const SquareMatrix & first, const SquareMatrix & second)
    {
        return SquareMatrix{first.m_Matrix * second.m_Matrix};
    }

    std::vector<double> operator*(const std::vector<double> & first, const SquareMatrix & second)
    {
        Eigen::VectorXd vec = Eigen::VectorXd::Map(first.data(), first.size());
        Eigen::VectorXd res = vec.transpose() * second.m_Matrix;
        return std::vector<double>(res.data(), res.data() + res.rows() * res.cols());
    }

    Eigen::SparseMatrix<double> SquareMatrix::getSparseMatrix() const
    {
        return m_Matrix;
    }

    std::vector<std::vector<double>> SquareMatrix::toVector() const
    {
		std::vector<std::vector<double>> result(m_size, std::vector<double>(m_size));
		for ( size_t i = 0u; i < m_size; ++i ) {
			for ( size_t j = 0u; j < m_size; ++j ) {
				result[i][j] = m_Matrix.coeff(i, j);
			}
		}
        return result;
    }

}   // namespace FenestrationCommon
