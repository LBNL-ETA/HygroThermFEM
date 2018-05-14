#include "SquareMatrix.hxx"

#include <Eigen/Dense>
#include <Eigen/LU>

namespace FenestrationCommon
{
    SquareMatrix::SquareMatrix(const std::size_t m_size) : m_size(m_size), m_Matrix(m_size, m_size)
    {
		m_Matrix.setZero();
    }

    SquareMatrix::SquareMatrix(const std::initializer_list<std::vector<double>> & tInput) :
        m_size(tInput.size()),
        m_Matrix(m_size, m_size)
    {
        auto i = 0u;
        for(auto vec : tInput)
        {
            for(auto j = 0u; j < vec.size(); ++j)
            {
                m_Matrix(i, j) = vec[j];
            }
            ++i;
        }
    }

    SquareMatrix::SquareMatrix(const std::vector<std::vector<double>> & tInput) :
        m_size(tInput.size()),
        m_Matrix(m_size, m_size)
    {
        for(auto i = 0u; i < tInput.size(); ++i)
        {
            for(auto j = 0u; j < tInput.size(); ++j)
            {
                m_Matrix(i, j) = tInput[i][j];
            }
        }
    }

    SquareMatrix::SquareMatrix(const std::vector<std::vector<double>> && tInput) :
        m_size(tInput.size()),
        m_Matrix(m_size, m_size)
    {
        for(auto i = 0u; i < tInput.size(); ++i)
        {
            for(auto j = 0u; j < tInput.size(); ++j)
            {
                m_Matrix(i, j) = tInput[i][j];
            }
        }
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
        for(auto i = 0u; i < m_size; ++i)
        {
            m_Matrix(i, i) = tInput[i];
        }
    }

    SquareMatrix::SquareMatrix(Eigen::MatrixXd && tMatrix) :
        m_size(std::size_t(tMatrix.innerSize())),
        m_Matrix(std::move(tMatrix))
    {}

	SquareMatrix::SquareMatrix( Eigen::MatrixXd & tMatrix ) :
		m_size(std::size_t(tMatrix.innerSize())),
		m_Matrix(std::move(tMatrix))
	{
	}

	SquareMatrix SquareMatrix::inverse() const
	{
		return SquareMatrix{ m_Matrix.inverse() };
	}

	double SquareMatrix::operator()( const std::size_t i, const std::size_t j ) const {
		return m_Matrix.coeffRef(i,j);
	}

	double & SquareMatrix::operator()( const std::size_t i, const std::size_t j ) {
		return m_Matrix(i,j);
	}

	SquareMatrix operator*( const SquareMatrix & first, const SquareMatrix & second ) {
		return SquareMatrix{ first.m_Matrix * second.m_Matrix };
	}

	SquareMatrix operator*( const SquareMatrix & other, const double value ) {
		return SquareMatrix{other.m_Matrix * value};
	}

	SquareMatrix operator*( const double value, const SquareMatrix & matrix ) {
		return SquareMatrix{matrix.m_Matrix * value};
	}

	SquareMatrix & SquareMatrix::operator*=( const SquareMatrix & other ) {
		m_Matrix = m_Matrix * other.m_Matrix;
		return *this;
	}

	SquareMatrix & SquareMatrix::operator+( const SquareMatrix & other ) {
		m_Matrix = m_Matrix + other.m_Matrix;
		return *this;
	}

	SquareMatrix & SquareMatrix::operator+=( const SquareMatrix & other ) {
		return operator+(other);
	}

	SquareMatrix & SquareMatrix::operator-( const SquareMatrix & other ) {
		m_Matrix = m_Matrix - other.m_Matrix;
		return *this;
	}

	SquareMatrix & SquareMatrix::operator-=( const SquareMatrix & other ) {
		return operator-(other);
	}

	SquareMatrix SquareMatrix::mmultRows( const std::vector< double > & tInput ) const {
    	Eigen::VectorXd vec = Eigen::VectorXd::Map(tInput.data(), tInput.size());
		const auto res = m_Matrix * vec.asDiagonal();
		return SquareMatrix{res};
	}

	std::vector< double > SquareMatrix::operator*( const std::vector< double > & tVec ) const
	{
		Eigen::VectorXd vec = Eigen::VectorXd::Map(tVec.data(), tVec.size());
		Eigen::VectorXd res = m_Matrix * vec;
		return std::vector< double >(res.data(), res.data() + res.rows() * res.cols());
	}

	std::vector< double > operator*( const std::vector< double > & first, const SquareMatrix & second ) {
		Eigen::VectorXd vec = Eigen::VectorXd::Map(first.data(), first.size());
		Eigen::VectorXd res = vec.transpose() * second.m_Matrix;
		return std::vector< double >(res.data(), res.data() + res.rows() * res.cols());
	}

	SquareMatrix SquareMatrix::addDiagonal( const std::vector< double > & tInput ) const {
		Eigen::VectorXd vec = Eigen::VectorXd::Map(tInput.data(), tInput.size());
		Eigen::MatrixXd res = m_Matrix + Eigen::MatrixXd{vec.asDiagonal()};
		return SquareMatrix{res};
	}

}   // namespace FenestrationCommon