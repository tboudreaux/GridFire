#include "gridfire/engine/types/jacobian.h"

#include <format>
#include <fstream>
#include <print>
#include <chrono>
#include <ranges>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>

#include "quill/LogMacros.h"

namespace gridfire::engine {
    NetworkJacobian::NetworkJacobian(
        const Eigen::SparseMatrix<double>& jacobianMatrix,
        const std::function<fourdst::atomic::Species(size_t)> &indexToSpeciesFunc
    ): m_jacobianMatrix(jacobianMatrix) {
        for (long int i = 0; i < jacobianMatrix.rows(); ++i) {
            fourdst::atomic::Species species = indexToSpeciesFunc(i);
            m_speciesToIndexMap[species] = i;
        }

    }

    NetworkJacobian::NetworkJacobian(
        const NetworkJacobian &jacobian
    ) : m_jacobianMatrix(jacobian.m_jacobianMatrix),
        m_speciesToIndexMap(jacobian.m_speciesToIndexMap),
        m_rank(jacobian.m_rank)
    {}

    NetworkJacobian::NetworkJacobian(
        NetworkJacobian &&jacobian
    ) noexcept : m_jacobianMatrix(std::move(jacobian.m_jacobianMatrix)),
                 m_speciesToIndexMap(std::move(jacobian.m_speciesToIndexMap)),
                 m_rank(jacobian.m_rank)
    {}

    NetworkJacobian & NetworkJacobian::operator=(
        NetworkJacobian &&jacobian
    ) noexcept {
        if (this != &jacobian) {
            m_jacobianMatrix = std::move(jacobian.m_jacobianMatrix);
            m_speciesToIndexMap = std::move(jacobian.m_speciesToIndexMap);
            m_rank = jacobian.m_rank;
        }
        return *this;
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
        if (i >= static_cast<size_t>(m_jacobianMatrix.rows()) || j >= static_cast<size_t>(m_jacobianMatrix.cols())) {
            throw std::out_of_range(std::format("Index ({}, {}) out of bounds in NetworkJacobian operator() for jacobian of shape ({}, {}).", i, j, m_jacobianMatrix.rows(), m_jacobianMatrix.cols()));
        }
        return m_jacobianMatrix.coeff(i, j);
    }

    void NetworkJacobian::set(const fourdst::atomic::Species &row, const fourdst::atomic::Species &col, const double value) {
        if (!m_speciesToIndexMap.contains(row) || !m_speciesToIndexMap.contains(col)) {
            throw std::out_of_range("Species not found in NetworkJacobian set().");
        }
        const size_t i = m_speciesToIndexMap.at(row);
        const size_t j = m_speciesToIndexMap.at(col);
        set(i, j, value);
    }

    void NetworkJacobian::set(const size_t i, const size_t j, const double value) {
        if (i >= static_cast<size_t>(m_jacobianMatrix.rows()) || j >= static_cast<size_t>(m_jacobianMatrix.cols())) {
            throw std::out_of_range(std::format("Index ({}, {}) out of bounds in NetworkJacobian set() for jacobian of shape ({}, {}).", i, j, m_jacobianMatrix.rows(), m_jacobianMatrix.cols()));
        }
        m_jacobianMatrix.coeffRef(i, j) = value;
    }

    void NetworkJacobian::set(const JacobianEntry &entry) {
        set(entry.first.first, entry.first.second, entry.second);
    }

    std::tuple<size_t, size_t> NetworkJacobian::shape() const {
        return {m_jacobianMatrix.rows(), m_jacobianMatrix.cols()};
    }

    size_t NetworkJacobian::nnz() const {
        return m_jacobianMatrix.nonZeros();
    }

    size_t NetworkJacobian::rank() const {
        if (m_rank) {
            return m_rank.value();
        }
        if (m_jacobianMatrix.rows() == 0 || m_jacobianMatrix.cols() == 0) {
            m_rank = 0;
        } else {
            Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
            solver.compute(m_jacobianMatrix);
            m_rank = solver.rank();
        }
        return m_rank.value();
    }

    bool NetworkJacobian::singular() const {
        const size_t rows = m_jacobianMatrix.rows();
        const size_t cols = m_jacobianMatrix.cols();
        const size_t minDim = (rows < cols) ? rows : cols;
        return m_rank < minDim;
    }

    std::vector<JacobianEntry> NetworkJacobian::infs() const {
        std::vector<JacobianEntry> infs;
        for (int k=0; k<m_jacobianMatrix.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_jacobianMatrix,k); it; ++it) {
                if (std::isinf(it.value())) {
                    fourdst::atomic::Species rowSpecies = std::ranges::find_if(
                        m_speciesToIndexMap,
                        [it](const auto& pair) {
                            return pair.second == static_cast<size_t>(it.row());
                        })->first;

                    fourdst::atomic::Species colSpecies = std::ranges::find_if(
                        m_speciesToIndexMap,
                        [it](const auto& pair) {
                            return pair.second == static_cast<size_t>(it.col());
                        })->first;

                    infs.emplace_back(std::make_pair(rowSpecies, colSpecies), it.value());
                }
            }
        }
        return infs;
    }

    std::vector<JacobianEntry> NetworkJacobian::nans() const {
        std::vector<JacobianEntry> nans;
        for (int k=0; k<m_jacobianMatrix.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(m_jacobianMatrix,k); it; ++it) {
                if (std::isnan(it.value())) {
                    fourdst::atomic::Species rowSpecies = std::ranges::find_if(
                        m_speciesToIndexMap,
                        [it](const auto& pair) {
                            return pair.second == static_cast<size_t>(it.row());
                        })->first;

                    fourdst::atomic::Species colSpecies = std::ranges::find_if(
                        m_speciesToIndexMap,
                        [it](const auto& pair) {
                            return pair.second == static_cast<size_t>(it.col());
                        })->first;

                    nans.emplace_back(std::make_pair(rowSpecies, colSpecies), it.value());
                }
            }
        }
        return nans;

    }

    Eigen::SparseMatrix<double> NetworkJacobian::data() const {
        return m_jacobianMatrix;
    }

    const std::unordered_map<fourdst::atomic::Species, size_t> & NetworkJacobian::mapping() const {
        return m_speciesToIndexMap;
    }

    void NetworkJacobian::to_csv(const std::string &filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            const size_t numSpecies = m_jacobianMatrix.rows();
            size_t i = 0;
            file << "# Species: ";
            for (size_t k = 0; k < numSpecies; ++k) {
                fourdst::atomic::Species species = std::ranges::find_if(
                    m_speciesToIndexMap,
                    [k](const auto& pair) {
                        return pair.second == k;
                    })->first;
                file << species.name();
                if (i < numSpecies - 1) {
                    file << ",";
                }
                i++;
            }
            const std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
            const std::time_t now_c = std::chrono::system_clock::to_time_t(now);

            file << "\n";
            file << "# Rows: " << numSpecies << ", Columns: " << numSpecies << "\n";
            file << "# Generated on " << std::ctime(&now_c) << "\n";
            file << "# Generated by GridFire NetworkJacobian::to_csv\n";

            for (size_t row = 0; row < numSpecies; ++row) {
                for (size_t col = 0; col < numSpecies; ++col) {
                    file << this->operator()(row, col);
                    if (col < numSpecies - 1) {
                        file << ",";
                    }
                }
                file << "\n";
            }
            file.close();
        }
    }

    NetworkJacobian regularize_jacobian(
        const NetworkJacobian& jacobian,
        const fourdst::composition::CompositionAbstract& comp,
        const std::optional<quill::Logger*> logger
    ) {
        const std::vector<JacobianEntry> infs = jacobian.infs();
        const std::vector<JacobianEntry> nans = jacobian.nans();
        if (infs.size() == 0 && nans.size() == 0) {
            return jacobian;
        }

        NetworkJacobian newJacobian = jacobian;
        for (const auto& [iSp, dSp] : infs | std::views::keys) {
            if (comp.getMolarAbundance(iSp) < MIN_ABUNDANCE_TO_CONTRIBUTE_TO_JACOBIAN || comp.getMolarAbundance(dSp) < MIN_ABUNDANCE_TO_CONTRIBUTE_TO_JACOBIAN) {
                newJacobian.set(iSp, dSp, 0.0);
                if (logger) {
                    LOG_TRACE_L1(logger.value(), "Regularized Jacobian entry ({}, {}) from inf to 0.0 due to low abundance.", iSp.name(), dSp.name());
                }
            }
        }
        for (const auto& [iSp, dSp] : nans | std::views::keys) {
            if (comp.getMolarAbundance(iSp) < MIN_ABUNDANCE_TO_CONTRIBUTE_TO_JACOBIAN || comp.getMolarAbundance(dSp) < MIN_ABUNDANCE_TO_CONTRIBUTE_TO_JACOBIAN) {
                newJacobian.set(iSp, dSp, 0.0);
                if (logger) {
                    LOG_TRACE_L1(logger.value(), "Regularized Jacobian entry ({}, {}) from inf to 0.0 due to low abundance.", iSp.name(), dSp.name());
                }
            }
        }

        return newJacobian;
    }
}
