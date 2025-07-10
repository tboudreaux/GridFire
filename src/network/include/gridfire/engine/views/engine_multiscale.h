#pragma once

#include "gridfire/engine/engine_abstract.h"
#include "gridfire/engine/views/engine_view_abstract.h"
#include "gridfire/engine/engine_graph.h"

#include "Eigen/Dense"
#include "unsupported/Eigen/NonLinearOptimization"

namespace gridfire {
    class MultiscalePartitioningEngineView final: public DynamicEngine, public EngineView<DynamicEngine> {
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
            double rho,
            double dt_control
        );

        void partitionNetwork(
            const NetIn& netIn,
            double dt_control
        );

        void exportToDot(const std::string& filename) const;

        [[nodiscard]] int getSpeciesIndex(const fourdst::atomic::Species &species) const override;

        [[nodiscard]] std::vector<double> mapNetInToMolarAbundanceVector(const NetIn &netIn) const override;

        std::vector<fourdst::atomic::Species> getFastSpecies() const;
        const std::vector<fourdst::atomic::Species>& getDynamicSpecies() const;

        void equilibrateNetwork(
            const std::vector<double>& Y,
            double T9,
            double rho,
            double dt_control
        );

        void equilibrateNetwork(
            const NetIn& netIn,
            const double dt_control
        );


    private:
        struct QSEGroup {
            std::vector<size_t> species_indices; ///< Indices of all species in this group.
            size_t seed_nucleus_index;           ///< Index of the one species that will be evolved dynamically.
            bool is_in_equilibrium = false;      ///< Flag set by flux analysis.
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

            int values() const { return m_qse_solve_indices.size(); }
            int inputs() const { return m_qse_solve_indices.size(); }

            int operator()(const InputType& v_qse, OutputType& f_qse) const;
            int df(const InputType& v_qse, JacobianType& J_qse) const;
        };

    private:
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        GraphEngine& m_baseEngine; ///< The base engine to which this view delegates calculations.
        std::vector<QSEGroup> m_qse_groups; ///< The list of identified equilibrium groups.
        std::vector<fourdst::atomic::Species> m_dynamic_species; ///< The simplified set of species presented to the solver.
        std::vector<size_t> m_dynamic_species_indices; ///< Indices mapping the dynamic species back to the master engine's list.
        std::unordered_map<size_t, std::vector<size_t>> m_connectivity_graph;

    private:
        std::unordered_set<size_t> identifyFastReactions(
            const std::vector<double>& Y_full,
            double T9,
            double rho,
            double dt_control
        ) const;

        void buildConnectivityGraph(
            const std::unordered_set<size_t>& fast_reaction_indices
        );

        void findConnectedComponents();

        void validateGroupsWithFluxAnalysis(
            const std::vector<double>& Y,
            double T9,
            double rho
        );

        std::pair<std::vector<fourdst::atomic::Species>, std::vector<size_t>> identifyDynamicSpecies(
            const std::vector<double>& Y,
            const std::vector<QSEGroup>& qse_groups,
            double T9,
            double rho
        ) const;

        std::vector<double> solveQSEAbundances(
            const std::vector<double> &Y_full,
            double T9,
            double rho
        );
    };
}
