#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_graph.h"

#include "unsupported/Eigen/NonLinearOptimization"

namespace gridfire {
    struct QSECacheConfig {
        double T9_tol; ///< Absolute tolerance to produce the same hash for T9.
        double rho_tol; ///< Absolute tolerance to produce the same hash for rho.
        double Yi_tol; ///< Absolute tolerance to produce the same hash for species abundances.
    };

    struct QSECacheKey {
        double m_T9;
        double m_rho;
        std::vector<double> m_Y; ///< Note that the ordering of Y must match the dynamic species indices in the view.

        std::size_t m_hash = 0; ///< Precomputed hash value for this key.

        // TODO: We should probably sort out how to adjust these from absolute to relative tolerances.
        QSECacheConfig m_cacheConfig = {
            1e-3, // Default tolerance for T9
            1e-1, // Default tolerance for rho
            1e-3  // Default tolerance for species abundances
        };

        QSECacheKey(
            const double T9,
            const double rho,
            const std::vector<double>& Y
        );

        size_t hash() const;

        static long bin(double value, double tol);

        bool operator==(const QSECacheKey& other) const;

    };
}

// Needs to be in this order (splitting gridfire namespace up) to avoid some issues with forward declarations and the () operator.
namespace std {
    template <>
    struct hash<gridfire::QSECacheKey> {
        /**
         * @brief Computes the hash of a QSECacheKey.
         * @param key The QSECacheKey to hash.
         * @return The pre-computed hash value of the key.
         */
        size_t operator()(const gridfire::QSECacheKey& key) const noexcept {
            // The hash is pre-computed, so we just return it.
            return key.m_hash;
        }
    };
} // namespace std

namespace gridfire {
    class MultiscalePartitioningEngineView final: public DynamicEngine, public EngineView<DynamicEngine> {
        typedef std::tuple<std::vector<fourdst::atomic::Species>, std::vector<size_t>, std::vector<fourdst::atomic::Species>, std::vector<size_t>> QSEPartition;
    public:
        explicit MultiscalePartitioningEngineView(GraphEngine& baseEngine);

        [[nodiscard]] const std::vector<fourdst::atomic::Species> & getNetworkSpecies() const override;

        [[nodiscard]] std::expected<StepDerivatives<double>, expectations::StaleEngineError> calculateRHSAndEnergy(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const override;

        void generateJacobianMatrix(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] double getJacobianMatrixEntry(
            int i_full,
            int j_full
        ) const override;

        void generateStoichiometryMatrix() override;

        [[nodiscard]] int getStoichiometryMatrixEntry(
            int speciesIndex,
            int reactionIndex
        ) const override;

        [[nodiscard]] double calculateMolarReactionFlow(
            const reaction::Reaction &reaction,
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] const reaction::LogicalReactionSet & getNetworkReactions() const override;

        void setNetworkReactions(
            const reaction::LogicalReactionSet &reactions
        ) override;

        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesTimescales(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] std::expected<std::unordered_map<fourdst::atomic::Species, double>, expectations::StaleEngineError> getSpeciesDestructionTimescales(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        fourdst::composition::Composition update(
            const NetIn &netIn
        ) override;

        bool isStale(const NetIn& netIn) override;

        void setScreeningModel(
            screening::ScreeningType model
        ) override;

        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;

        const DynamicEngine & getBaseEngine() const override;

        std::vector<std::vector<size_t>> analyzeTimescalePoolConnectivity(
            const std::vector<std::vector<size_t>> &timescale_pools,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const;

        void partitionNetwork(
            const std::vector<double>& Y,
            double T9,
            double rho
        );

        void partitionNetwork(
            const NetIn& netIn
        );

        void exportToDot(
            const std::string& filename,
            const std::vector<double>& Y,
            const double T9,
            const double rho
        ) const;

        [[nodiscard]] int getSpeciesIndex(const fourdst::atomic::Species &species) const override;

        [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const override;

        [[nodiscard]] PrimingReport primeEngine(const NetIn &netIn) override;

        [[nodiscard]] std::vector<fourdst::atomic::Species> getFastSpecies() const;
        [[nodiscard]] const std::vector<fourdst::atomic::Species>& getDynamicSpecies() const;

        fourdst::composition::Composition equilibrateNetwork(
            const std::vector<double> &Y,
            double T9,
            double rho
        );

        fourdst::composition::Composition equilibrateNetwork(
            const NetIn &netIn
        );


    private:
        struct QSEGroup {
            std::set<size_t> species_indices; ///< Indices of all species in this group.
            bool is_in_equilibrium = false;      ///< Flag set by flux analysis.
            std::set<size_t> algebraic_indices; ///< Indices of algebraic species in this group.
            std::set<size_t> seed_indices; ///< Indices of dynamic species in this group.
            double mean_timescale; ///< Mean timescale of the group.

            bool operator<(const QSEGroup& other) const;
            bool operator>(const QSEGroup& other) const;
            bool operator==(const QSEGroup& other) const;
            bool operator!=(const QSEGroup& other) const;

            friend std::ostream& operator<<(std::ostream& os, const QSEGroup& group) {
                os << "QSEGroup(species_indices: [";
                int count = 0;
                for (const auto& idx : group.species_indices) {
                    os << idx;
                    if (count < group.species_indices.size() - 1) {
                        os << ", ";
                    }
                    count++;
                }
                count = 0;
                os << "], is_in_equilibrium: " << group.is_in_equilibrium
                   << ", mean_timescale: " << group.mean_timescale
                   << " s, algebraic_indices: [";
                for (const auto& idx : group.algebraic_indices) {
                    os << idx;
                    if (count < group.algebraic_indices.size() - 1) {
                        os << ", ";
                    }
                    count++;
                }
                count = 0;
                os << "], seed_indices: [";
                for (const auto& idx : group.seed_indices) {
                    os << idx;
                    if (count < group.seed_indices.size() - 1) {
                        os << ", ";
                    }
                    count++;
                }
                os << "])";
                return os;
            }
        };

        struct EigenFunctor {
            using InputType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
            using OutputType = Eigen::Matrix<double, Eigen::Dynamic, 1>;
            using JacobianType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
            enum {
                InputsAtCompileTime = Eigen::Dynamic,
                ValuesAtCompileTime = Eigen::Dynamic
            };

            MultiscalePartitioningEngineView* m_view;
            const std::vector<size_t>& m_qse_solve_indices;
            const std::vector<double>& m_Y_full_initial;
            const double m_T9;
            const double m_rho;
            const Eigen::VectorXd& m_Y_scale;

            EigenFunctor(
                MultiscalePartitioningEngineView& view,
                const std::vector<size_t>& qse_solve_indices,
                const std::vector<double>& Y_full_initial,
                const double T9,
                const double rho,
                const Eigen::VectorXd& Y_scale
            ) :
            m_view(&view),
            m_qse_solve_indices(qse_solve_indices),
            m_Y_full_initial(Y_full_initial),
            m_T9(T9),
            m_rho(rho),
            m_Y_scale(Y_scale) {}

            [[nodiscard]] int values() const { return m_qse_solve_indices.size(); }
            [[nodiscard]] int inputs() const { return m_qse_solve_indices.size(); }

            int operator()(const InputType& v_qse, OutputType& f_qse) const;
            int df(const InputType& v_qse, JacobianType& J_qse) const;
        };


        struct CacheStats {
            enum class operators {
                CalculateRHSAndEnergy,
                GenerateJacobianMatrix,
                CalculateMolarReactionFlow,
                GetSpeciesTimescales,
                GetSpeciesDestructionTimescales,
                Other,
                All
            };

            std::map<operators, std::string> operatorsNameMap = {
                {operators::CalculateRHSAndEnergy, "calculateRHSAndEnergy"},
                {operators::GenerateJacobianMatrix, "generateJacobianMatrix"},
                {operators::CalculateMolarReactionFlow, "calculateMolarReactionFlow"},
                {operators::GetSpeciesTimescales, "getSpeciesTimescales"},
                {operators::GetSpeciesDestructionTimescales, "getSpeciesDestructionTimescales"},
                {operators::Other, "other"}
            };

            size_t m_hit = 0;
            size_t m_miss = 0;

            std::map<operators, size_t> m_operatorHits = {
                {operators::CalculateRHSAndEnergy, 0},
                {operators::GenerateJacobianMatrix, 0},
                {operators::CalculateMolarReactionFlow, 0},
                {operators::GetSpeciesTimescales, 0},
                {operators::GetSpeciesDestructionTimescales, 0},
                {operators::Other, 0}
            };


            friend std::ostream& operator<<(std::ostream& os, const CacheStats& stats) {
                os << "CacheStats(hit: " << stats.m_hit << ", miss: " << stats.m_miss << ")";
                return os;
            }

            void hit(const operators op=operators::Other) {
                if (op==operators::All) {
                    throw std::invalid_argument("Cannot use 'All' as an operator for hit/miss.");
                }
                m_hit++;
                m_operatorHits[op]++;
            }
            void miss(const operators op=operators::Other) {
                if (op==operators::All) {
                    throw std::invalid_argument("Cannot use 'All' as an operator for hit/miss.");
                }
                m_miss++;
                m_operatorHits[op]++;
            }

            [[nodiscard]] size_t hits(const operators op=operators::All) const {
                if (op==operators::All) {
                    return m_hit;
                }
                return m_operatorHits.at(op);
            }
            [[nodiscard]] size_t misses(const operators op=operators::All) const {
                if (op==operators::All) {
                    return m_miss;
                }
                return m_operatorHits.at(op);
            }
        };


    private:
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
        GraphEngine& m_baseEngine; ///< The base engine to which this view delegates calculations.
        std::vector<QSEGroup> m_qse_groups; ///< The list of identified equilibrium groups.
        std::vector<fourdst::atomic::Species> m_dynamic_species; ///< The simplified set of species presented to the solver.
        std::vector<size_t> m_dynamic_species_indices; ///< Indices mapping the dynamic species back to the master engine's list.
        std::vector<fourdst::atomic::Species> m_algebraic_species; ///< Species that are algebraic in the QSE groups.
        std::vector<size_t> m_algebraic_species_indices; ///< Indices of algebraic species in the full network.

        std::vector<size_t> m_activeSpeciesIndices; ///< Indices of active species in the full network.
        std::vector<size_t> m_activeReactionIndices; ///< Indices of active reactions in the full network.

        // TODO: Enhance the hashing for the cache to consider not just T and rho but also the current abundance in some careful way that automatically ignores small changes (i.e. network should only be repartitioned sometimes)
        std::unordered_map<QSECacheKey, std::vector<double>> m_qse_abundance_cache; ///< Cache for QSE abundances based on T9 and rho.
        mutable CacheStats m_cacheStats; ///< Statistics for the QSE abundance cache.




    private:
        std::vector<std::vector<size_t>> partitionByTimescale(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        ) const;

        std::unordered_map<size_t, std::vector<size_t>> buildConnectivityGraph(
            const std::unordered_set<size_t> &fast_reaction_indices
        ) const;

        std::vector<QSEGroup> validateGroupsWithFluxAnalysis(
            const std::vector<QSEGroup> &candidate_groups,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const;

        std::vector<double> solveQSEAbundances(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        );

        size_t identifyMeanSlowestPool(
            const std::vector<std::vector<size_t>>& pools,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const;

        std::unordered_map<size_t, std::vector<size_t>> buildConnectivityGraph(
            const std::vector<size_t>& species_pool
        ) const;

        std::vector<QSEGroup> constructCandidateGroups(
            const std::vector<std::vector<size_t>>& candidate_pools,
            const std::vector<double>& Y,
            double T9, double rho
        ) const;
    };
}

