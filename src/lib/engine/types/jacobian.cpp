#include "gridfire/engine/types/jacobian.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>

namespace gridfire {
    NetworkJacobian::NetworkJacobian(
        const Eigen::SparseMatrix<double>& jacobianMatrix,
        const std::function<fourdst::atomic::Species(size_t)> &indexToSpeciesFunc
    ): m_jacobianMatrix(jacobianMatrix) {
        for (size_t i = 0; i < jacobianMatrix.rows(); ++i) {
            fourdst::atomic::Species species = indexToSpeciesFunc(i);
            m_speciesToIndexMap[species] = i;
        }

        if (m_jacobianMatrix.rows() == 0 || m_jacobianMatrix.cols() == 0) {
            m_rank = 0;
        } else {
            Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
            solver.compute(m_jacobianMatrix);
            m_rank = solver.rank();
        }
    }

    double NetworkJacobian::operator()(const fourdst::atomic::Species &row, const fourdst::atomic::Species &col) const {
        if (!m_speciesToIndexMap.contains(row) || !m_speciesToIndexMap.contains(col)) {
            throw std::out_of_range("Species not found in NetworkJacobian operator().");
        }
        const size_t i = m_speciesToIndexMap.at(row);
        const size_t j = m_speciesToIndexMap.at(col);
        return this->operator()(i, j);
    }

    double NetworkJacobian::operator()(const size_t i, const size_t j) const {
        if (i >= m_jacobianMatrix.rows() || j >= m_jacobianMatrix.cols()) {
            throw std::out_of_range("Index out of bounds in NetworkJacobian operator().");
        }
        return m_jacobianMatrix.coeff(i, j);
    }

    std::tuple<size_t, size_t> NetworkJacobian::shape() const {
        return {m_jacobianMatrix.rows(), m_jacobianMatrix.cols()};
    }

    size_t NetworkJacobian::nnz() const {
        return m_jacobianMatrix.nonZeros();
    }

    size_t NetworkJacobian::rank() const {
        return m_rank;
    }

    bool NetworkJacobian::singular() const {
        const size_t rows = m_jacobianMatrix.rows();
        const size_t cols = m_jacobianMatrix.cols();
        const size_t minDim = (rows < cols) ? rows : cols;
        return m_rank < minDim;
    }


}