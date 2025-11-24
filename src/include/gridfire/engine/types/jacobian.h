/**
 * @file jacobian.h
 * @brief Wrapper for handling the network Jacobian matrix in GridFire. Currently uses Eigen's SparseMatrix.
 */
#pragma once

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/composition/composition_abstract.h"
#include "quill/Logger.h"
#include <Eigen/SparseCore>

#include <tuple>
#include <functional>
#include <unordered_map>
#include <optional>

namespace gridfire::engine {
    constexpr double MIN_ABUNDANCE_TO_CONTRIBUTE_TO_JACOBIAN = 1e-100; ///< Minimum abundance for a species to contribute to the Jacobian
    using JacobianEntry = std::pair<std::pair<fourdst::atomic::Species, fourdst::atomic::Species>, double>; ///< Represents an entry in the Jacobian matrix

    /**
     * @class NetworkJacobian
     * @brief A wrapper class for handling the network Jacobian matrix.
     *
     * This class encapsulates an Eigen::SparseMatrix to represent the Jacobian matrix of a reaction network.
     * It provides methods for accessing and modifying matrix entries using species identifiers, as well as
     * utility functions for analyzing the matrix (e.g., checking for singularity, counting non-zero entries).
     */
    class NetworkJacobian {
    public:
        /**
         * @brief Constructs a NetworkJacobian with the given sparse matrix and species index mapping function.
         *
         * @param jacobianMatrix The sparse matrix representing the Jacobian.
         * @param indexToSpeciesFunc A function that maps matrix indices to species identifiers.
         */
        explicit NetworkJacobian(
            const Eigen::SparseMatrix<double>& jacobianMatrix,
            const std::function<fourdst::atomic::Species(size_t)> &indexToSpeciesFunc
        );

        /**
         * @brief Copy constructor for NetworkJacobian.
         *
         * @param jacobian The NetworkJacobian instance to copy from.
         */
        NetworkJacobian(const NetworkJacobian& jacobian);

        /**
         * @brief Copy assignment operator for NetworkJacobian.
         *
         * @param jacobian The NetworkJacobian instance to copy from.
         * @return Reference to the assigned NetworkJacobian instance.
         */
        NetworkJacobian(NetworkJacobian&& jacobian) noexcept;

        /**
         * @brief Move constructor for NetworkJacobian.
         *
         * @param jacobian The NetworkJacobian instance to move from.
         * @return Reference to the assigned NetworkJacobian instance.
         */
        NetworkJacobian& operator=(NetworkJacobian&& jacobian) noexcept;

        /**
         * @brief Accesses the value at the specified row and column corresponding to the given species.
         * @param row Species for the row
         * @param col Species for the column
         * @return value at the specified position
         */
        double operator()(
            const fourdst::atomic::Species& row,
            const fourdst::atomic::Species& col
        ) const;

        /**
         * @brief Accesses the value at the specified row and column indices.
         * @param i Row index
         * @param j Column index
         * @return value at the specified position
         */
        double operator()(
            size_t i,
            size_t j
        ) const;

        /**
         * @brief Sets the value at the specified row and column corresponding to the given species.
         * @param row Row species
         * @param col Column species
         * @param value Value to set
         */
        void set(
            const fourdst::atomic::Species& row,
            const fourdst::atomic::Species& col,
            double value
        );

        /**
         * @brief Sets the value at the specified row and column indices.
         * @param i Row index
         * @param j Column index
         * @param value Value to set
         */
        void set(
            size_t i,
            size_t j,
            double value
        );

        /**
         * @brief Sets the value in the Jacobian matrix based on a JacobianEntry.
         * @param entry The JacobianEntry containing row species, column species, and value.
         */
        void set(
            const JacobianEntry &entry
        );

        /**
         * @brief Retrieves the shape of the Jacobian matrix as a tuple (rows, columns).
         * @return A tuple containing the number of rows and columns.
         */
        std::tuple<size_t, size_t> shape() const;

        /**
         * @brief Retrieves the rank of the Jacobian matrix.
         * @return The rank of the matrix.
         * @note Rank is computed using QR decomposition and cached for efficiency. The rank is computed the first time
         * this method is called and stored for subsequent calls. If any set operation is performed on the matrix, the cached rank is invalidated.
         */
        size_t rank() const;

        /**
         * @brief Get the number of non-zero entries in the Jacobian matrix.
         * @return The number of non-zero entries.
         */
        size_t nnz() const;

        /**
         * @brief Checks if the Jacobian matrix is singular.
         * @return True if the matrix is singular, false otherwise.
         */
        bool singular() const;

        /**
         * @brief Retrieves all entries in the Jacobian matrix that are infinite.
         * @return A vector of JacobianEntry representing infinite entries.
         */
        [[nodiscard]] std::vector<JacobianEntry> infs() const;

        /**
         * @brief Retrieves all entries in the Jacobian matrix that are NaN (Not a Number).
         * @return A vector of JacobianEntry representing NaN entries.
         */
        [[nodiscard]] std::vector<JacobianEntry> nans() const;

        /**
         * @brief Provides access to the underlying sparse matrix data.
         * @return The Eigen::SparseMatrix representing the Jacobian.
         */
        [[nodiscard]] Eigen::SparseMatrix<double> data() const;

        /**
         * @brief Provides access to the species-to-index mapping.
         * @return A constant reference to the species-to-index mapping.
         */
        [[nodiscard]] const std::unordered_map<fourdst::atomic::Species, size_t>& mapping() const;

        /**
         * @brief Exports the Jacobian matrix to a CSV file.
         * @param filename The name of the CSV file to export to.
         */
        void to_csv(const std::string& filename) const;

    private:
        Eigen::SparseMatrix<double> m_jacobianMatrix;
        std::unordered_map<fourdst::atomic::Species, size_t> m_speciesToIndexMap;

        mutable std::optional<size_t> m_rank = std::nullopt;
    };

    /**
     * @brief Regularizes the given Jacobian matrix based on the provided composition.
     *
     * This function applies regularization techniques to the Jacobian matrix to improve its numerical stability.
     * Specifically any (row, column) entries corresponding to species who's abundance is below a threshold of
     * MIN_ABUNDANCE_TO_CONTRIBUTE_TO_JACOBIAN in the provided composition will be set to zero if those entries
     * were either infinite or NaN.
     *
     * @param jacobian The NetworkJacobian to be regularized.
     * @param comp The composition used for regularization.
     * @param logger Optional logger for logging regularization steps and information.
     * @return A new NetworkJacobian instance representing the regularized Jacobian matrix.
     */
    NetworkJacobian regularize_jacobian(
        const NetworkJacobian& jacobian,
        const fourdst::composition::CompositionAbstract& comp,
        std::optional<quill::Logger*> logger = std::nullopt
    );
}