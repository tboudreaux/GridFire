// Validation: GridFire against MESA (via bbq) on identical networks and burn histories.
//
// For each case in cases.json with a reference in reference/<case>/ (see generate_reference.py), GridFire integrates
// the same chain of constant (T, rho) steps that bbq ran, starting each step from the previous step's composition,
// and compares mass fractions after every step (in dex) and the nuclear energy released (cumulatively, and in total).
// Tolerances are set per case in cases.json from the measured agreement, with notes on the known causes of the
// differences (mainly MESA's curated rate choices versus GridFire's REACLIB rates).
//
// Cases without reference data are skipped. A reference produced for a different network (reactions.txt differs from
// the network GridFire builds today) fails: regenerate the references.
//
// A per-case comparison table is written to mesa_comparison_<case>.csv in the working directory for inspection.
#include "mesa_case.h"
#include "gf_test_support.h"

#include "gridfire/screening/screening_types.h"

#include "fourdst/atomic/species.h"
#include "nlohmann/json.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using gridfire::test::mesa::BurnCase;

namespace {
    struct ReferenceTable {
        std::vector<std::string> columns;
        std::vector<std::vector<double>> rows;

        [[nodiscard]] size_t column(const std::string& name) const {
            const auto it = std::ranges::find(columns, name);
            if (it == columns.end()) throw std::runtime_error("reference has no column " + name);
            return static_cast<size_t>(it - columns.begin());
        }
    };

    ReferenceTable read_reference(const std::filesystem::path& path) {
        std::ifstream in(path);
        ReferenceTable table;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::stringstream ss(line);
            std::string cell;
            if (table.columns.empty()) {
                while (std::getline(ss, cell, ',')) table.columns.push_back(cell);
                continue;
            }
            std::vector<double> row;
            while (std::getline(ss, cell, ',')) row.push_back(std::stod(cell));
            table.rows.push_back(std::move(row));
        }
        return table;
    }

    std::vector<std::string> read_lines(const std::filesystem::path& path) {
        std::ifstream in(path);
        std::vector<std::string> lines;
        for (std::string line; std::getline(in, line);) {
            if (!line.empty()) lines.push_back(line);
        }
        return lines;
    }

    /// Agreement thresholds, set per case in cases.json ("tolerances"); see the notes there for why each is what it is.
    struct Tolerances {
        double species_dex = 0.1;     ///< Max |log10 X_GridFire - log10 X_MESA| for mass fractions above species_floor
        double species_floor = 1e-3;  ///< Mass fractions below this (in both codes) are reported but not asserted
        std::map<std::string, double> species_dex_overrides; ///< Per-species species_dex (documented exceptions)
        double energy_rel = 0.05;     ///< Max relative error of the cumulative energy, once it exceeds energy_floor
        double energy_floor = 0.1;    ///< Fraction of the final cumulative energy before which energy is not asserted
        double final_energy_rel = 0.02; ///< Max relative error of the total energy released over the whole burn
    };

    Tolerances tolerances_for(const std::string& case_name) {
        Tolerances tol;
        std::ifstream in(gridfire::test::mesa::validation_dir() / "cases.json");
        const auto doc = nlohmann::json::parse(in);
        for (const auto& c : doc.at("cases")) {
            if (c.at("name") != case_name || !c.contains("tolerances")) continue;
            const auto& t = c.at("tolerances");
            tol.species_dex = t.value("species_dex", tol.species_dex);
            tol.species_floor = t.value("species_floor", tol.species_floor);
            tol.energy_rel = t.value("energy_rel", tol.energy_rel);
            tol.energy_floor = t.value("energy_floor", tol.energy_floor);
            tol.final_energy_rel = t.value("final_energy_rel", tol.final_energy_rel);
            if (t.contains("species_dex_overrides")) {
                for (const auto& [sp, v] : t.at("species_dex_overrides").items()) tol.species_dex_overrides[sp] = v.get<double>();
            }
        }
        return tol;
    }

    std::vector<BurnCase> all_cases() {
        return gridfire::test::mesa::load_cases(gridfire::test::mesa::validation_dir() / "cases.json");
    }

    class MesaComparison : public ::testing::TestWithParam<BurnCase> {
    protected:
        static void SetUpTestSuite() { gridfire::test::quiet_logging(); }
    };
}

TEST_P(MesaComparison, AgreesWithMESA) {
    const BurnCase& burn_case = GetParam();
    const auto ref_dir = gridfire::test::mesa::validation_dir() / "reference" / burn_case.name;
    if (!std::filesystem::exists(ref_dir / "reference.csv")) {
        GTEST_SKIP() << "No MESA reference for " << burn_case.name << "; run tests/validation/mesa/generate_reference.py";
    }

    // Reactions MESA cannot evaluate were removed from the case when the reference was made; remove them here too
    const auto excluded = gridfire::test::mesa::read_excluded_reactions(ref_dir / "excluded_reactions.txt");
    auto network = gridfire::test::mesa::build_network(burn_case, excluded);

    // --- The reference must have been produced for this exact network ---
    const auto ids = gridfire::test::mesa::sorted_reaction_ids(network);
    const auto ref_ids = read_lines(ref_dir / "reactions.txt");
    if (ids != ref_ids) {
        std::set<std::string> now(ids.begin(), ids.end()), then(ref_ids.begin(), ref_ids.end());
        std::string added, removed;
        for (const auto& id : now) if (!then.contains(id)) added += " " + id;
        for (const auto& id : then) if (!now.contains(id)) removed += " " + id;
        FAIL() << "GridFire's network for " << burn_case.name << " differs from the one the MESA reference was made for; "
               << "regenerate the references.\n  added:" << added << "\n  removed:" << removed;
    }

    // --- ...and with the same screening: GridFire's weak screening is MESA's 'salpeter', bare is MESA's '' ---
    {
        std::ifstream provenance_in(ref_dir / "provenance.json");
        ASSERT_TRUE(provenance_in) << "missing provenance.json for " << burn_case.name;
        const auto provenance = nlohmann::json::parse(provenance_in);
        const std::string mesa_screening = provenance.value("screening_mode", "<not recorded>");
        const auto gf_screening = network.engine->getScreeningModel(*network.context);
        const std::string expected = gf_screening == gridfire::screening::ScreeningType::WEAK ? "salpeter" : "";
        ASSERT_EQ(mesa_screening, expected)
            << "The MESA reference for " << burn_case.name << " was made with screening_mode '" << mesa_screening
            << "', which does not match GridFire's screening; regenerate the references.";
    }

    const ReferenceTable ref = read_reference(ref_dir / "reference.csv");
    ASSERT_GE(ref.rows.size(), 2u);
    const size_t c_dt = ref.column("dt"), c_energy = ref.column("energy_step");
    std::vector<std::pair<fourdst::atomic::Species, size_t>> tracked;
    for (size_t c = 4; c < ref.columns.size(); ++c) {
        tracked.emplace_back(fourdst::atomic::species.at(ref.columns[c]), c);
    }

    // --- Initial compositions must agree before anything is integrated ---
    for (const auto& [sp, c] : tracked) {
        const double gf = network.initial_composition.contains(sp) ? network.initial_composition.getMassFraction(sp) : 0.0;
        ASSERT_NEAR(gf, ref.rows[0][c], 1e-12) << "initial mass fraction of " << sp.name();
    }

    // --- Start from MESA's initial state, nucleus for nucleus ---
    // MESA's mass fractions are nucleon fractions, X_i = A_i Y_i, whereas GridFire's are true mass fractions,
    // X_i = m_i Y_i / sum_j m_j Y_j. The same list of X therefore means different molar abundances (0.8% fewer protons
    // in GridFire for a solar mixture), which shows up directly as a difference in energy released. Give GridFire
    // Y_i = X_i / A_i, MESA's abundances, and compare species as A_i Y_i, MESA's X, below.
    const auto& network_species = network.engine->getNetworkSpecies(*network.context);
    std::vector<double> Y0(network_species.size(), 0.0);
    for (const auto& [sp, c] : tracked) {
        const auto it = std::ranges::find(network_species, sp);
        if (it != network_species.end()) Y0[static_cast<size_t>(it - network_species.begin())] = ref.rows[0][c] / sp.a();
    }
    const auto nucleon_fraction = [](const fourdst::composition::CompositionAbstract& comp, const fourdst::atomic::Species& sp) {
        return comp.contains(sp) ? sp.a() * comp.getMolarAbundance(sp) : 0.0;
    };

    // --- Integrate the same step chain ---
    const Tolerances tol = tolerances_for(burn_case.name);
    double final_energy_mesa = 0.0;
    for (size_t r = 1; r < ref.rows.size(); ++r) final_energy_mesa += ref.rows[r][c_energy];

    std::ofstream table(std::format("mesa_comparison_{}.csv", burn_case.name));
    table << "# X columns are nucleon fractions A_i Y_i (MESA's mass fraction convention) for both codes\n";
    table << "age,energy_mesa,energy_gridfire,energy_rel_err";
    for (const auto& [sp, c] : tracked) table << ",X_mesa_" << sp.name() << ",X_gridfire_" << sp.name();
    table << "\n";

    std::map<std::string, double> worst_dex;
    double worst_energy_err = 0.0, cumulative_mesa = 0.0, cumulative_gridfire = 0.0;
    int reported = 0;
    const auto report = [&](const std::string& message) { if (reported++ < 20) ADD_FAILURE() << message; };

    fourdst::composition::Composition comp(network_species, Y0);
    for (size_t r = 1; r < ref.rows.size(); ++r) {
        const auto& row = ref.rows[r];
        gridfire::solver::PointSolverContext sctx(*network.context);
        gridfire::solver::PointSolver solver(*network.engine);
        gridfire::NetIn in;
        in.composition = comp;
        in.temperature = std::pow(10.0, burn_case.logT);
        in.density = std::pow(10.0, burn_case.logRho);
        in.energy = 0.0;
        in.tMax = row[c_dt];
        in.dt0 = std::min(1e-12, in.tMax);
        const gridfire::NetOut out = solver.evaluate(sctx, in);
        comp = out.composition;

        const double e_ref = row[c_energy];
        cumulative_mesa += e_ref;
        cumulative_gridfire += out.energy;
        table << std::format("{:.10e},{:.10e},{:.10e},{:.4e}", row[0], e_ref, out.energy, gridfire::test::rel_diff(out.energy, e_ref));
        if (cumulative_mesa >= tol.energy_floor * final_energy_mesa) {
            const double err = gridfire::test::rel_diff(cumulative_gridfire, cumulative_mesa);
            worst_energy_err = std::max(worst_energy_err, err);
            if (err > tol.energy_rel) {
                report(std::format("{} step {} (age {:.4e} s): cumulative energy GridFire {:.6e} vs MESA {:.6e} erg/g (rel {:.3f} > {})",
                                   burn_case.name, r, row[0], cumulative_gridfire, cumulative_mesa, err, tol.energy_rel));
            }
        }

        for (const auto& [sp, c] : tracked) {
            const double x_ref = row[c];
            const double x_gf = nucleon_fraction(comp, sp);
            table << std::format(",{:.10e},{:.10e}", x_ref, x_gf);
            if (std::max(x_ref, x_gf) < tol.species_floor) continue;
            const std::string name(sp.name());
            const double dex = std::abs(std::log10(std::max(x_gf, 1e-300)) - std::log10(std::max(x_ref, 1e-300)));
            worst_dex[name] = std::max(worst_dex[name], dex);
            const double limit = tol.species_dex_overrides.contains(name) ? tol.species_dex_overrides.at(name) : tol.species_dex;
            if (dex > limit) {
                report(std::format("{} step {} (age {:.4e} s): X({}) GridFire {:.6e} vs MESA {:.6e} ({:.3f} dex > {})",
                                   burn_case.name, r, row[0], name, x_gf, x_ref, dex, limit));
            }
        }
        table << "\n";
    }

    const double final_err = gridfire::test::rel_diff(cumulative_gridfire, cumulative_mesa);
    EXPECT_LE(final_err, tol.final_energy_rel) << std::format("{}: total energy GridFire {:.6e} vs MESA {:.6e} erg/g",
                                                               burn_case.name, cumulative_gridfire, cumulative_mesa);

    // --- Human-readable summary ---
    std::string summary = std::format("{}: total energy rel err {:.2e} (tol {}), worst cumulative energy rel err {:.2e} (tol {})\n",
                                      burn_case.name, final_err, tol.final_energy_rel, worst_energy_err, tol.energy_rel);
    std::vector<std::pair<std::string, double>> sorted(worst_dex.begin(), worst_dex.end());
    std::ranges::sort(sorted, [](const auto& a, const auto& b) { return a.second > b.second; });
    for (const auto& [name, dex] : sorted) summary += std::format("  {:<7} worst {:.3f} dex\n", name, dex);
    std::cout << summary;
    RecordProperty("total_energy_rel_err", std::format("{:.3e}", final_err));
}

INSTANTIATE_TEST_SUITE_P(Cases, MesaComparison, ::testing::ValuesIn(all_cases()),
                         [](const ::testing::TestParamInfo<BurnCase>& info) { return info.param.name; });
