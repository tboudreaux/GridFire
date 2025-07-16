#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_graph.h"

#include "unsupported/Eigen/NonLinearOptimization"

namespace gridfire {
    class MultiscalePartitioningEngineView final: public DynamicEngine, public EngineView<DynamicEngine> {
        typedef std::tuple<std::vector<fourdst::atomic::Species>, std::vector<size_t>, std::vector<fourdst::atomic::Species>, std::vector<size_t>> QSEPartition;
    public:
        explicit MultiscalePartitioningEngineView(GraphEngine& baseEngine);

        [[nodiscard]] const std::vector<fourdst::atomic::Species> & getNetworkSpecies() const override;

        [[nodiscard]] StepDerivatives<double> calculateRHSAndEnergy(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        void generateJacobianMatrix(
            const std::vector<double> &Y_dynamic,
            double T9,
            double rho
        ) override;

        [[nodiscard]] double getJacobianMatrixEntry(
            int i,
            int j
        ) const override;

        void generateStoichiometryMatrix() override;

        [[nodiscard]] int getStoichiometryMatrixEntry(
            int speciesIndex,
            int reactionIndex
        ) const override;

        [[nodiscard]] double calculateMolarReactionFlow(
            const reaction::Reaction &reaction,
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        [[nodiscard]] const reaction::LogicalReactionSet & getNetworkReactions() const override;

        [[nodiscard]] std::unordered_map<fourdst::atomic::Species, double> getSpeciesTimescales(
            const std::vector<double> &Y,
            double T9,
            double rho
        ) const override;

        void update(
            const NetIn &netIn
        ) override;

        void setScreeningModel(
            screening::ScreeningType model
        ) override;

        [[nodiscard]] screening::ScreeningType getScreeningModel() const override;

        const DynamicEngine & getBaseEngine() const override;

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
                   << ", algebraic_indices: [";
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

    private:
        quill::Logger* m_logger = LogManager::getInstance().getLogger("log");
        GraphEngine& m_baseEngine; ///< The base engine to which this view delegates calculations.
        std::vector<QSEGroup> m_qse_groups; ///< The list of identified equilibrium groups.
        std::vector<fourdst::atomic::Species> m_dynamic_species; ///< The simplified set of species presented to the solver.
        std::vector<size_t> m_dynamic_species_indices; ///< Indices mapping the dynamic species back to the master engine's list.
        std::vector<fourdst::atomic::Species> m_algebraic_species; ///< Species that are algebraic in the QSE groups.
        std::vector<size_t> m_algebraic_species_indices; ///< Indices of algebraic species in the full network.
        std::unordered_map<size_t, std::vector<size_t>> m_connectivity_graph;

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

        std::vector<QSEGroup> constructCandidateGroups(
            const std::vector<std::vector<size_t>>& timescale_pools,
            const std::vector<double>& Y,
            double T9,
            double rho
        ) const;
    };
}
