# LBM-CDE Source-Feedback TRT Derivation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (- [ ]) syntax for tracking.

**Goal:** Derive, verify, and document a D2Q9/D2Q9 LBM-CDE source-feedback TRT algorithm, including the BGK limit, effective relaxation rates, boundary magic analysis, and D2Q9 fourth-order temperature conditions.

**Architecture:** Use one shared symbolic kernel for lattice moments, parity projectors, collision matrices, and series coefficients. Prove second-order recovery analytically, verify constant-coefficient fourth-order terms through two independent routes, and treat boundary residuals as separate symbolic problems. Integrate only reviewed results into a Chinese Markdown/LaTeX report.

**Tech Stack:** Python 3, SymPy 1.14, NumPy 2.4, mpmath 1.3, standard-library unittest, XeLaTeX, Poppler.

## Global Constraints

- The main model is D2Q9 flow plus D2Q9 temperature; D2Q5 is verification-only.
- The four user-supplied PDFs in the original workspace are authoritative source material and remain untracked.
- Xs/2DRBOpenmpLBMCDE.F90 and Xs/2DRBOpenaccLBMCDE.F90 are prior reproductions and never theoretical evidence.
- Flow and temperature have separate nominal TRT pairs.
- Source terms use parity-specific trapezoidal factors; one shared sourcePref is forbidden.
- Nominal TRT rates and source-feedback effective moment rates must never be conflated.
- A boundary constant is called magic only if all independent residual coefficients vanish.
- Fourth-order isotropy and complete fourth-order cancellation are reported separately.
- Dubois-Lallemand D2Q9 printed coefficients are audited, not silently corrected.
- Every major derivation receives an implementer pass, an adversarial review, and an independent limiting-case check.
- Production Fortran solvers are not modified in this plan.
- Every RED step must end in an expected assertion failure; delayed imports or a
  temporary test helper must convert a missing implementation into `FAIL`, not
  leave an import, syntax, or collection `ERROR`.

---

### Task 1: Lattice Algebra and Evidence Ledger

**Files:**
- Create: tools/derivation/__init__.py
- Create: tools/derivation/lattice.py
- Create: tests/derivation/test_lattice.py
- Create: docs/derivation/evidence-ledger.md

**Interfaces:**
- Produces: Lattice, d2q5(), d2q9(), parity_projectors(), raw_moment(), hermite_moment().
- Consumes: no project code.

- [ ] **Step 1: Write failing lattice tests**

~~~python
class D2Q9MomentTests(unittest.TestCase):
    def test_opposite_is_involution(self):
        lat = d2q9()
        self.assertEqual(
            tuple(lat.opposite[lat.opposite[i]] for i in range(9)),
            tuple(range(9)),
        )

    def test_weights_and_lambda_constraints(self):
        lat = d2q9()
        self.assertEqual(sum(lat.weights), 1)
        self.assertEqual(sum(lat.lambda_t), 0)
        self.assertEqual(raw_moment(lat.lambda_t, lat.velocities, (2, 0)), Rational(1, 3))
        self.assertEqual(raw_moment(lat.lambda_t, lat.velocities, (0, 2)), Rational(1, 3))
        self.assertEqual(raw_moment(lat.lambda_t, lat.velocities, (1, 1)), 0)

    def test_parity_projectors(self):
        p_plus, p_minus = parity_projectors(d2q9())
        self.assertEqual(p_plus * p_plus, p_plus)
        self.assertEqual(p_minus * p_minus, p_minus)
        self.assertEqual(p_plus * p_minus, zeros(9))
~~~

- [ ] **Step 2: Run the tests and confirm failure**

Run: python -m unittest tests.derivation.test_lattice -v  
Expected: the delayed-import assertion reports the missing lattice behavior as
`FAIL`, with zero syntax or collection errors.

- [ ] **Step 3: Implement exact lattice definitions**

~~~python
@dataclass(frozen=True)
class Lattice:
    velocities: tuple[tuple[int, int], ...]
    weights: tuple[Expr, ...]
    lambda_t: tuple[Expr, ...]
    opposite: tuple[int, ...]
    cs2: Expr

def d2q9() -> Lattice: ...
def d2q5() -> Lattice: ...
def parity_projectors(lattice: Lattice) -> tuple[Matrix, Matrix]: ...
def raw_moment(coeffs, velocities, powers: tuple[int, int]) -> Expr: ...
~~~

Use exact SymPy Rational values. Record every paper equation used in evidence-ledger.md with PDF filename, printed page, equation number, and purpose.

- [ ] **Step 4: Run tests**

Run: python -m unittest tests.derivation.test_lattice -v  
Expected: all lattice tests pass.

- [ ] **Step 5: Commit**

~~~powershell
git add -- tools/derivation tests/derivation docs/derivation/evidence-ledger.md
git commit -m "test: establish exact D2Q lattice algebra"
~~~

### Task 2: Parity-Resolved Sources and TRT Collision

**Files:**
- Create: tools/derivation/sources.py
- Create: tools/derivation/collision.py
- Create: tests/derivation/test_sources_collision.py
- Create: docs/derivation/chapters/02-trt-source-design.md

**Interfaces:**
- Consumes: Lattice and parity_projectors from Task 1.
- Produces: flow_source(), scalar_source(), trt_collision(), bgk_collision(), source_moment_table().

- [ ] **Step 1: Write failing source-moment and BGK-limit tests**

~~~python
class SourceCollisionTests(unittest.TestCase):
    def test_flow_source_parity_moments(self):
        s_plus, s_minus = flow_source(symbolic_state())
        self.assertEqual(moment(s_minus, 0), 0)
        self.assertEqual(moment(s_minus, 1, axis=0), Fx)
        self.assertEqual(moment(s_plus, 1, axis=0), 0)

    def test_scalar_source_even_second_moment(self):
        r_plus, _ = scalar_source(symbolic_state())
        self.assertEqual(moment(r_plus, 0), Q)
        self.assertEqual(moment(r_plus, 2, axes=(0, 0)), cs2 * Q)

    def test_equal_rates_recover_bgk_componentwise(self):
        trt = trt_collision(sample_populations(), sfp=s, sfm=s)
        bgk = bgk_collision(sample_populations(), s=s)
        self.assertEqual(simplify(trt - bgk), zeros(9, 1))
~~~

- [ ] **Step 2: Run tests and confirm failure**

Run: python -m unittest tests.derivation.test_sources_collision -v  
Expected: assertion failure for the unimplemented source/collision behavior.

- [ ] **Step 3: Implement parity-resolved sources**

Implement the flow and scalar source formulas from the approved design. The returned objects must contain raw source, plus/minus projections, and exact moment tables through third order.

- [ ] **Step 4: Implement transformed TRT collision**

~~~python
def trt_collision(h_tilde, h_eq, source_plus, source_minus, s_plus, s_minus, dt=1):
    return (
        h_tilde
        - s_plus * P_PLUS * (h_tilde - h_eq)
        - s_minus * P_MINUS * (h_tilde - h_eq)
        + dt * ((1 - s_plus/2) * source_plus + (1 - s_minus/2) * source_minus)
    )
~~~

Add explicit tests that net momentum and net heat source after one collision are dt*F and dt*Q when half-source nonequilibrium moments are used.

- [ ] **Step 5: Run tests and write derivation chapter**

Run: python -m unittest tests.derivation.test_sources_collision -v  
Expected: all tests pass.  
Chapter must show the projection algebra, source moment table, trapezoidal operator derivation, and componentwise BGK limit.

- [ ] **Step 6: Commit**

~~~powershell
git add -- tools/derivation tests/derivation docs/derivation/chapters/02-trt-source-design.md
git commit -m "feat: derive parity-resolved TRT sources"
~~~

### Task 3: Effective Rates and Second-Order Recovery

**Files:**
- Create: tools/derivation/effective_rates.py
- Create: tests/derivation/test_effective_rates.py
- Create: docs/derivation/chapters/03-effective-rates.md
- Create: docs/derivation/chapters/04-second-order-recovery.md

**Interfaces:**
- Consumes: source and collision operators from Task 2.
- Produces: scalar_flux_effective_rate(), shear_effective_rate(), second_order_residual_table().

- [ ] **Step 1: Write failing effective-rate tests**

~~~python
class EffectiveRateTests(unittest.TestCase):
    def test_scalar_flux_feedback(self):
        sigma_eff = scalar_flux_effective_rate(sigma_o, chi, pressure=0)
        self.assertEqual(simplify(sigma_eff - (1-chi)*sigma_o), 0)

    def test_diffusivity_is_pressure_independent_at_second_order(self):
        sigma_eff = scalar_flux_effective_rate(sigma_o, chi, pressure=P)
        diffusivity = simplify((cs2 + P) * sigma_eff)
        self.assertEqual(diffusivity, cs2 * (1-chi) * sigma_o)

    def test_bgk_feedback_limit(self):
        self.assertEqual(
            shear_effective_rate(sigma, chi_s=0),
            sigma,
        )
~~~

- [ ] **Step 2: Run tests and confirm failure**

Run: python -m unittest tests.derivation.test_effective_rates -v  
Expected: assertion failure for the unimplemented effective-rate behavior.

- [ ] **Step 3: Implement local feedback elimination**

Eliminate the LBM-CDE nonequilibrium gradient and strain closures algebraically before computing effective relaxation parameters. Return both the effective physical-moment rate and the unchanged ghost rates.

- [ ] **Step 4: Derive second-order moment equations**

Write explicit order-by-order equations for continuity, momentum, and temperature. The residual table must track p*grad(T), T*F, Q*u, u*F, time-varying Q/F, and the first omitted order.

- [ ] **Step 5: Run tests and perform adversarial equation review**

Run: python -m unittest tests.derivation.test_effective_rates -v  
Expected: all tests pass.  
Independent reviewer must check conservation, half-source terms, and pressure cancellation.

- [ ] **Step 6: Commit**

~~~powershell
git add -- tools/derivation tests/derivation docs/derivation/chapters/03-effective-rates.md docs/derivation/chapters/04-second-order-recovery.md
git commit -m "feat: derive effective rates and second-order recovery"
~~~

### Task 4: D2Q5 Reference Generator

**Files:**
- Create: tools/derivation/series.py
- Create: tools/derivation/d2q5_reference.py
- Create: tests/derivation/test_d2q5_reference.py
- Create: docs/derivation/chapters/A-d2q5-verifier.md

**Interfaces:**
- Produces: matrix_series(), hydro_eigen_series(), modified_equation_coefficients(), d2q5_coefficients().
- Serves as the independent reference gate for Task 5.

- [ ] **Step 1: Write failing D2Q5 reference tests**

~~~python
class D2Q5ReferenceTests(unittest.TestCase):
    def test_isotropy_relation(self):
        k40, k22 = d2q5_coefficients(alpha, h1, h3, h4)
        relation = simplify((k22 - 2*k40).subs(h4, 1/(6*h1)))
        self.assertEqual(relation, 0)

    def test_quartic_trt_point(self):
        coeffs = d2q5_coefficients(alpha, 1/sqrt(12), 1/sqrt(3), 1/sqrt(3))
        self.assertEqual(tuple(map(simplify, coeffs)), (0, 0))
~~~

- [ ] **Step 2: Run tests and confirm failure**

Run: python -m unittest tests.derivation.test_d2q5_reference -v  
Expected: assertion failure for the unimplemented series behavior.

- [ ] **Step 3: Implement two coefficient routes**

Route A expands the hydrodynamic eigenvalue of the exact amplification matrix. Route B implements the recursive moment/Taylor formulas. Both return k2, C40, and C22 in exact SymPy form.

- [ ] **Step 4: Reproduce the paper equations**

The chapter must reproduce Dubois-Lallemand D2Q5 equations (40)-(42) and (55), show that isotropy differs from cancellation, and map the Wang/Luo symbols without importing any D2Q5 constant into D2Q9.

- [ ] **Step 5: Run tests**

Run: python -m unittest tests.derivation.test_d2q5_reference -v  
Expected: all tests pass and both routes agree.

- [ ] **Step 6: Commit**

~~~powershell
git add -- tools/derivation tests/derivation docs/derivation/chapters/A-d2q5-verifier.md
git commit -m "test: reproduce D2Q5 fourth-order reference"
~~~

### Task 5: D2Q9 Fourth-Order Generator and Dubois Audit

**Files:**
- Create: tools/derivation/d2q9_temperature.py
- Create: tools/derivation/fourier_verify.py
- Create: tests/derivation/test_d2q9_fourth_order.py
- Create: docs/derivation/chapters/06-d2q9-fourth-order.md

**Interfaces:**
- Consumes: series kernel from Task 4.
- Produces: d2q9_amplification(), d2q9_coefficients(),
  solve_d2q9_quartic_conditions(), dubois_printed_residual(), directional_fit().

- [ ] **Step 1: Write failing D2Q9 tests**

~~~python
class D2Q9FourthOrderTests(unittest.TestCase):
    def test_printed_k22_conflict_is_explicit(self):
        residual = dubois_printed_residual(
            sigma_o=1/sqrt(12), sigma_e=1/sqrt(3), xi=Rational(1,3)
        )
        self.assertNotEqual(simplify(residual), 0)

    def test_matrix_and_recursive_routes_agree(self):
        matrix = d2q9_coefficients(generic_parameters(), route="matrix")
        recursive = d2q9_coefficients(generic_parameters(), route="recursive")
        self.assertEqual(
            tuple(simplify(a-b) for a, b in zip(matrix, recursive)),
            (0, 0),
        )

    def test_solved_quartic_family_has_no_q4_term(self):
        point = solve_d2q9_quartic_conditions(physical_branch=True)
        c40, c22 = d2q9_coefficients(point, route="matrix")
        self.assertEqual((simplify(c40), simplify(c22)), (0, 0))

    def test_axis_and_diagonal_agree_at_cancellation(self):
        point = solve_d2q9_quartic_conditions(physical_branch=True)
        axis, diagonal = directional_fit(point)
        self.assertLess(abs(axis.q4), 1e-40)
        self.assertLess(abs(diagonal.q4), 1e-40)
~~~

- [ ] **Step 2: Run tests and confirm failure**

Run: python -m unittest tests.derivation.test_d2q9_fourth_order -v  
Expected: assertion failure for the unimplemented D2Q9 behavior.

- [ ] **Step 3: Implement source-feedback D2Q9 matrices**

Support three cases separately: chi=0, exact external gradient source, and local nonequilibrium gradient feedback. Keep sigma_flux_eff, nominal odd ghost sigma, and even sigma as distinct arguments.

- [ ] **Step 4: Generate and cross-check C40/C22**

Use exact symbolic coefficients where tractable and at least 80-digit numerical directional fits as an independent check. Report C22-2*C40 and the pair (C40,C22) separately.

- [ ] **Step 5: Audit the printed equation**

Document the exact nonzero residual obtained from the printed k22 formula and the zero residual obtained from the regenerated matrix/Taylor routes. Do not silently alter the source paper.

- [ ] **Step 6: Run tests and commit**

Run: python -m unittest tests.derivation.test_d2q9_fourth_order -v  
Expected: all tests pass.

~~~powershell
git add -- tools/derivation tests/derivation docs/derivation/chapters/06-d2q9-fourth-order.md
git commit -m "feat: regenerate D2Q9 fourth-order conditions"
~~~

### Task 6: Boundary Residuals and Magic Conditions

**Files:**
- Create: tools/derivation/boundary.py
- Create: tests/derivation/test_boundary.py
- Create: docs/derivation/chapters/05-boundary-magic.md

**Interfaces:**
- Consumes: effective-rate functions and lattice definitions.
- Produces: velocity_bb_residual(), temperature_abb_residual(), adiabatic_bb_residual(), classify_magic().

- [ ] **Step 1: Write failing boundary tests**

~~~python
class BoundaryResidualTests(unittest.TestCase):
    def test_classical_velocity_limit_after_documented_mapping(self):
        residual = velocity_bb_residual(
            chi_s=0, force_mode="uniform", normalization="glt2002"
        )
        condition = residual.solve_zero_conditions()
        self.assertEqual(condition[sigma_even*sigma_odd], Rational(3,16))

    def test_temperature_manufactured_limit_is_solved_not_assumed(self):
        residual = temperature_abb_residual(
            pressure=0, velocity=0, force=0, q_mode="constant"
        )
        calibration = residual.solve_zero_conditions()
        self.assertEqual(
            tuple(map(simplify, residual.substitute(calibration).coefficients)),
            tuple(0 for _ in residual.coefficients),
        )

    def test_pressure_wall_value_leaves_independent_residual(self):
        residual = temperature_abb_residual(pressure="variable")
        self.assertTrue(residual.has_independent("normal_pressure_gradient"))
~~~

- [ ] **Step 2: Run tests and confirm failure**

Run: python -m unittest tests.derivation.test_boundary -v  
Expected: assertion failure for the unimplemented boundary behavior.

- [ ] **Step 3: Implement chain Taylor residuals**

Derive residuals for velocity halfway BB, temperature Dirichlet ABB, and adiabatic BB through the first nonzero wall-error order. Keep pressure, source, time, and tangential derivatives as independent symbolic jets.

- [ ] **Step 4: Classify each result**

classify_magic() returns one of: universal_magic, restricted_calibration, boundary_correction_required, or no_single_magic. Corner rules are assessed separately and cannot be fixed by a scalar relaxation product.

- [ ] **Step 5: Run tests and write chapter**

Run: python -m unittest tests.derivation.test_boundary -v  
Expected: all tests pass.  
Chapter must state every assumption beside each candidate constant.

- [ ] **Step 6: Commit**

~~~powershell
git add -- tools/derivation tests/derivation docs/derivation/chapters/05-boundary-magic.md
git commit -m "feat: derive source-aware boundary residuals"
~~~

### Task 7: Parameter Compatibility and Algorithm Specification

**Files:**
- Create: tools/derivation/parameters.py
- Create: tests/derivation/test_parameters.py
- Create: docs/derivation/chapters/07-parameter-compatibility.md
- Create: docs/derivation/chapters/08-algorithm.md

**Interfaces:**
- Consumes: fourth-order and boundary conditions.
- Produces: solve_parameter_family(), admissibility_report(), algorithm pseudocode.

- [ ] **Step 1: Write failing parameter tests**

~~~python
class ParameterTests(unittest.TestCase):
    def test_rates_stay_in_open_stability_interval(self):
        family = solve_parameter_family(target_kappa=Rational(1,1000))
        self.assertTrue(all(0 < rate < 2 for rate in family.rates))

    def test_incompatible_constraints_are_reported_not_hidden(self):
        report = admissibility_report(overconstrained_example())
        self.assertEqual(report.status, "no_feasible_solution")
        self.assertTrue(report.minimal_extension)
~~~

- [ ] **Step 2: Run tests and confirm failure**

Run: python -m unittest tests.derivation.test_parameters -v  
Expected: assertion failure for the unimplemented parameter behavior.

- [ ] **Step 3: Implement exact and numerical solvers**

Solve physical transport, boundary, and fourth-order constraints in priority order. Return symbolic families when available and high-precision interval-checked numerical solutions otherwise.

- [ ] **Step 4: Write implementation-ready algorithm**

The chapter fixes this order: half-source macroscopic reconstruction; equilibrium; local nonequilibrium closures; raw parity sources; parity-specific collision; streaming; BB/ABB/adiabatic/corner boundaries; diagnostics. Include a BGK switch and invariants for source moments and parity.

- [ ] **Step 5: Run tests and commit**

Run: python -m unittest tests.derivation.test_parameters -v  
Expected: all tests pass.

~~~powershell
git add -- tools/derivation tests/derivation docs/derivation/chapters/07-parameter-compatibility.md docs/derivation/chapters/08-algorithm.md
git commit -m "feat: solve TRT compatibility and specify algorithm"
~~~

### Task 8: Integrate and Cross-Review the Complete Report

**Files:**
- Create: docs/lbm-cde-trt-derivation.md
- Create: docs/lbm-cde-trt-derivation.tex
- Create: docs/derivation/review-matrix.md
- Modify: all chapter drafts only as required by review findings.

**Interfaces:**
- Consumes: all reviewed chapters, scripts, and test outputs.
- Produces: one coherent Chinese report and a source-to-claim review matrix.

- [ ] **Step 1: Assemble the report**

Required chapter order: executive conclusion levels; notation and target PDEs; inverse TRT design; trapezoidal discrete LBE; effective rates; second-order recovery; BGK limit; boundary residuals; D2Q9 fourth-order analysis; compatibility; algorithm pseudocode; verification; limitations; references; D2Q5 appendix.

- [ ] **Step 2: Build a requirement-to-evidence matrix**

Every major claim is labeled strictly_proved, restricted_model, numerical_evidence, or unresolved. Each row links the report equation, source paper equation, verification function, and test name.

- [ ] **Step 3: Run the full derivation suite**

Run: python -m unittest discover -s tests/derivation -v  
Expected: all tests pass with zero failures and zero errors.

- [ ] **Step 4: Conduct three-agent cross-review**

Reviewer A checks TRT/CE/source algebra. Reviewer B checks all boundary assumptions and constants. Reviewer C checks D2Q5/D2Q9 high-order algebra and reproduced program outputs. Critical and important findings are fixed and re-reviewed.

- [ ] **Step 5: Commit**

~~~powershell
git add -- docs/lbm-cde-trt-derivation.md docs/lbm-cde-trt-derivation.tex docs/derivation
git commit -m "docs: integrate LBM-CDE TRT derivation report"
~~~

### Task 9: Build, Render, and Audit the PDF

**Files:**
- Create: tools/derivation/build_report.ps1
- Create: output/pdf/lbm-cde-trt-derivation.pdf
- Create: output/pdf/rendered/ page PNGs as temporary QA artifacts; do not commit PNGs.

**Interfaces:**
- Consumes: docs/lbm-cde-trt-derivation.tex.
- Produces: final PDF with verified layout.

- [ ] **Step 1: Write the deterministic build script**

~~~powershell
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force output/pdf | Out-Null
xelatex -interaction=nonstopmode -halt-on-error -output-directory output/pdf docs/lbm-cde-trt-derivation.tex
xelatex -interaction=nonstopmode -halt-on-error -output-directory output/pdf docs/lbm-cde-trt-derivation.tex
~~~

- [ ] **Step 2: Build the PDF**

Run: powershell -ExecutionPolicy Bypass -File tools/derivation/build_report.ps1  
Expected: exit 0 and output/pdf/lbm-cde-trt-derivation.pdf exists.

- [ ] **Step 3: Render every PDF page**

Run: pdftoppm -png output/pdf/lbm-cde-trt-derivation.pdf output/pdf/rendered/page  
Expected: one PNG per PDF page.

- [ ] **Step 4: Visually inspect and revise**

Inspect every page for clipped equations, broken Chinese glyphs, overfull boxes, unreadable tables, duplicated numbering, or citation problems. Rebuild after every layout change.

- [ ] **Step 5: Run final verification**

Run: python -m unittest discover -s tests/derivation -v  
Expected: all tests pass.  
Run: git diff --check  
Expected: no whitespace errors.  
Run: pdfinfo output/pdf/lbm-cde-trt-derivation.pdf  
Expected: valid PDF metadata and nonzero page count.

- [ ] **Step 6: Commit**

~~~powershell
git add -- tools/derivation/build_report.ps1 output/pdf/lbm-cde-trt-derivation.pdf
git commit -m "docs: publish verified LBM-CDE TRT report"
~~~

## Final Review Gate

- Generate one whole-branch review package from commit 4d88c15 to HEAD.
- Verify the full objective requirement-by-requirement against current files and fresh command output.
- Fix every critical or important review finding and re-run its covering tests.
- Keep the goal active until the report, programs, tests, and rendered PDF all satisfy the design completion standard.
