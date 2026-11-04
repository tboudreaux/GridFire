# GridFire vs MESA validation

Once GridFire's output starts changing, the question that matters most is whether it still agrees with MESA. This
suite answers it: GridFire and MESA (driven through [bbq](https://github.com/rjfarmer/bbq)) burn **the same network**
through **the same history**, and the results are compared step by step.

## How it works

* `cases.json` defines constant (T, ρ) burns. Each case is a chain of steps of duration `10^t` with `t` on a uniform
  grid, each starting from the previous step's composition (bbq's hydrostatic mode).
* `gf_mesa_setup` (built with the tests) builds each case's GridFire network and exports it as a MESA `.net` file,
  together with its sorted reaction list and initial composition. Both GridFire and MESA therefore integrate exactly
  GridFire's network.
* `generate_reference.py` runs bbq on those files and writes `reference/<case>/` (MESA output in GridFire species
  names, the `.net` file, the reaction list, the inlist, and `provenance.json` recording MESA version, bbq version and
  the hashes of the patches applied). MESA cannot evaluate a few exotic reactions GridFire includes (reverse rates
  missing from its reaclib index, weak decays missing from its tables) and reports them as non-fatal errors while
  continuing with bad rates. The generator scans every bbq log, removes such reactions from the case in **both** codes,
  and records them with the reason in `excluded_reactions.txt`; the test rebuilds the network with the same exclusions.
* `mesa_comparison_test` (`meson test --suite validation`) rebuilds each network, **fails if it no longer matches the
  reference's reaction list** (regenerate the references), integrates the same step chain in GridFire, and compares mass
  fractions after every step (in dex, above a mass-fraction floor), the cumulative energy released, and the total
  energy. Each case's `"tolerances"` in `cases.json` are regression guards set at about 1.5x the measured agreement,
  with notes on the known causes of the remaining differences (mostly MESA's curated rates versus GridFire's REACLIB
  rates). Cases without reference data are skipped.
* The test writes `mesa_comparison_<case>.csv` in its working directory; `plot_comparison.py` turns those into figures.

Both codes use weak (Salpeter 1954) screening: GridFire's default `engine.screening = "weak"` and MESA's
`screening_mode = 'salpeter'`, set in bbq's `&nuclear` group. bbq reads `screening_mode` only from `&nuclear`; anywhere
else it is silently ignored and bbq uses its default, Chugunov screening (this affected the GridFire paper's inlists,
which set it in a `&controls` group). The comparison test fails if `provenance.json` records a different screening from
GridFire's. bbq runs with `BBQ_DISABLE_EOS=1`, as in the paper.

MESA's mass fractions are nucleon fractions (X = A·Y), GridFire's are true mass fractions (X = m·Y / Σ m·Y), so the same
list of X gives different molar abundances (0.8% fewer protons in GridFire for a solar mixture, which shows up as ~0.6%
less energy). The comparison test therefore starts GridFire from MESA's molar abundances, Y = X/A, and compares species
as A·Y, MESA's X.

To compare unscreened instead, MESA r23.05.1 needs `mesa-r23.05.1-no-screening.patch` (`rates/private/screen.f90`):
without it MESA's `screening_mode = ''` sets every screened rate to zero, because the screening context's temperature is
left uninitialised. `mesa_rates/` (a MESA `net` driver) and `gf_mesa_rates` print per-reaction rates, Q values and
neutrino energies at one state for both codes, for tracking down differences reaction by reaction.

## Regenerating references (needs MESA + bbq)

References only need regenerating when GridFire's network for a case changes or the MESA setup changes; they can be
committed or uploaded as CI artifacts so that machines without MESA still run the comparison.

1. Get MESA r23.05.1 (a version bbq supports; on GitHub the tag is `23.05.1`):
   `GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 --branch 23.05.1 https://github.com/MESAHub/mesa.git && git lfs pull`.
2. Apply `mesa-r23.05.1-net_burn.patch` (this folder) to `net/private/net_burn.f90` (and, for unscreened runs,
   `mesa-r23.05.1-no-screening.patch`), then build MESA with `MESA_DIR` set. The patch skips the EOS call unless `BBQ_DISABLE_EOS` is exactly `False` (as in the GridFire paper, so timings
   are comparable) and, unlike the paper's patch, sets `eta` to zero instead of leaving it uninitialised when skipped.
   (The paper's `validation/ManuscriptFigures/BBQCompare/diffs/net_burn.f90.patch` was made against a different MESA
   release and does not apply to r23.05.1.) Creating an empty `skip_test` file in `$MESA_DIR` skips MESA's module
   self-tests, whose stored outputs can differ in the last digits with newer compilers.
3. Build bbq against it with `lib_bbq.f90.patch` and `lib_hydrostatic.f90.patch` from
   `validation/ManuscriptFigures/BBQCompare/diffs` applied (timing instrumentation).
4. From a GridFire build directory:

   ```bash
   python ../tests/validation/mesa/generate_reference.py --bbq /path/to/bbq --setup tests/validation/mesa/gf_mesa_setup
   meson test --suite validation
   python ../tests/validation/mesa/plot_comparison.py
   ```
