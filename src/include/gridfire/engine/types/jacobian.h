#pragma once

#include "fourdst/atomic/atomicSpecies.h"
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>

#include <tuple>
#include <functional>
#include <unordered_map>

namespace gridfire {
    using JacobianEntry = std::pair<std::pair<fourdst::atomic::Species, fourdst::atomic::Species>, double>;

    class NetworkJacobian {
    public:
        explicit NetworkJacobian(
            const Eigen::SparseMatrix<double>& jacobianMatrix,
            const std::function<fourdst::atomic::Species(size_t)> &indexToSpeciesFunc
        );

        NetworkJacobian(const NetworkJacobian& jacobian);
        NetworkJacobian(NetworkJacobian&& jacobian) noexcept;
        NetworkJacobian& operator=(NetworkJacobian&& jacobian) noexcept;

        double operator()(
            const fourdst::atomic::Species& row,
            const fourdst::atomic::Species& col
        ) const;

        double operator()(
            size_t i,
            size_t j
        ) const;

        void set(
            const fourdst::atomic::Species& row,
            const fourdst::atomic::Species& col,
            double value
        );

        void set(
            size_t i,
            size_t j,
            double value
        );

        void set(
            const JacobianEntry &entry
        );

        std::tuple<size_t, size_t> shape() const;
        size_t rank() const;
        size_t nnz() const;

        bool singular() const;

        [[nodiscard]] std::vector<JacobianEntry> infs() const;
        [[nodiscard]] std::vector<JacobianEntry> nans() const;

        [[nodiscard]] Eigen::SparseMatrix<double> data() const;
        [[nodiscard]] const std::unordered_map<fourdst::atomic::Species, size_t>& mapping() const;

    private:
        Eigen::SparseMatrix<double> m_jacobianMatrix;
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap;
        size_t m_rank;
    };
}