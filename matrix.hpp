#include <algorithm>
#include <valarray>
#include <vector>

template <class element_type>
class MatrixVA {
public:
    MatrixVA(const size_t n):m_size(n) {
	m_storage.resize(n*n);
    }
    MatrixVA(const size_t n, const element_type val):m_size(n) {
	m_storage.resize(n*n);
	m_storage = val;
    }
    // Accurate copy constructor
    MatrixVA(const size_t n, const MatrixVA& from):m_size(n) {
	m_storage.resize(n*n);
	const size_t from_size = from.m_size;
	if (n == from_size) {
	    m_storage = from.m_storage;
	}
	else {
	    const int to = std::min(n, from_size);
	    for (int i = 0; i < to; i++)
		for (int j = 0; j < to; j++)
		    m_storage[i*n + j] = from.m_storage[i*from_size + j];
	}
    }

    MatrixVA<element_type>& operator+=(const MatrixVA<element_type> &rhs)
    {
	this->m_storage += rhs.m_storage;
	return *this;
    }

    MatrixVA<element_type> operator+(const MatrixVA<element_type> &other)
    {
	return MatrixVA(*this) += other;
    }

    MatrixVA<element_type>& operator-=(const MatrixVA<element_type> &rhs)
    {
	this->m_storage -= rhs.m_storage;
	return *this;
    }

    MatrixVA<element_type> operator-(const MatrixVA<element_type> &other)
    {
	return MatrixVA(*this) -= other;
    }

    MatrixVA<element_type> operator*(const MatrixVA<element_type> &other)
    {
	int i, k, j;
	const int n = this->m_size;
	MatrixVA res(n);
	for (i = 0; i < n; i++) {
	    for (k = 0; k < n; k++) {
		const element_type aink = m_storage[i*n + k];
		for (j = 0; j < n; j++)
		    res(i,j) += aink * other(k,j);
	    }
	}
	return res;
    }

    element_type& operator()(const int row, const int col)
    {
	return m_storage[row * m_size + col];
    }

    const element_type& operator()(const int row, const int col) const
    {
	return m_storage[row * m_size + col];
    }

private:
    std::valarray<element_type> m_storage;
    size_t m_size;
};

template <class element_type>
class MatrixVect {
public:
    MatrixVect(const size_t n):m_size(n) {
	m_storage.resize(n*n);
    }
    MatrixVect(const size_t n, const element_type val):m_size(n) {
	m_storage.resize(n*n);
	std::fill(m_storage.begin(), m_storage.end(), val);
    }

    element_type& operator()(const int row, const int col)
    {
	return m_storage[row * m_size + col];
    }

private:
    std::vector<element_type> m_storage;
    size_t m_size;
};
