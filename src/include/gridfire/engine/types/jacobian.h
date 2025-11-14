#pragma once

#include "fourdst/atomic/atomicSpecies.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>

#include <tuple>
#include <functional>
#include <unordered_map>

namespace gridfire {
    class NetworkJacobian {
    public:
        explicit NetworkJacobian(
            const Eigen::SparseMatrix<double>& jacobianMatrix,
            const std::function<fourdst::atomic::Species(size_t)> &indexToSpeciesFunc
        );

        double operator()(
            const fourdst::atomic::Species& row,
            const fourdst::atomic::Species& col
        ) const;

        double operator()(
            size_t i,
            size_t j
        ) const;

        std::tuple<size_t, size_t> shape() const;
        size_t rank() const;
        size_t nnz() const;

        bool singular() const;

    private:
        Eigen::SparseMatrix<double> m_jacobianMatrix;
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap;
        size_t m_rank;
    };
}