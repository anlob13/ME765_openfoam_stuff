# periodicHill fixes — exact changes



Every change below was applied **identically to both `mlSim/` and `spaldingSim/`**
(except the MATLAB scripts, which are shared at the top level) — both cases must
share the same physical setup for the comparison between them to mean anything.

---

## 1. `mlSim/Allrun` and `spaldingSim/Allrun`

**Old:**
```sh
runParallel simpleFoam
```

**New:**
```sh
runParallel foamRun
```

**Why:** `system/controlDict` already declared `application foamRun; solver
incompressibleFluid;` — the modern, transient, PIMPLE-based solver framework
(matching `deltaT`, `adjustTimeStep`, `maxCo` in the same file) — and
`constant/momentumTransport` was set to `LES { model WALE; ... }`. But `Allrun`
was still invoking `simpleFoam`, OpenFOAM's steady-state, RAS-only SIMPLE-algorithm
solver, which has never supported LES in any OpenFOAM version and does not read
the transient controls in `controlDict` at all. Running an LES/WALE closure
through a steady solver produces a non-physical mean field for an inherently
unsteady, separated flow like periodic hill — this is the direct cause of the
LES results being "completely off," independent of the wall model. `foamRun`
correctly dispatches on `momentumTransport`'s `simulationType`, so this same
line is also correct for RAS runs (e.g. if `momentumTransport` is switched back
to RAS for a RANS comparison).

---

## 2. `mlSim/constant/transportProperties` and `spaldingSim/constant/transportProperties`

**Old:**
```
nu              2.643e-6;
```

**New:**
```
nu              9.438414e-05;
```

**Why:** The paper specifies Re_H = U_b·H/ν = 10,595. This case's own convention
sets H=1 and U_b=1 (see `blockMeshDict`'s `scale 1` + vertices normalized by hill
height, and `fvOptions`' `meanVelocityForce` with `Ubar (1 0 0)`), so
ν = 1/Re_H = 9.438414×10⁻⁵. The old value (2.643×10⁻⁶) implies an effective
Re_H ≈ 378,400 — about 36× higher than the paper's target, and far beyond what
this mesh's resolution (120×40×60, matching the paper's own "fine mesh") could
resolve even before any of the other fixes here.

---

## 3. `mlSim/system/controlDict` and `spaldingSim/system/controlDict`

**Old:**
```
endTime         0.5;            // set based on several flow-through times + averaging window

deltaT          1e-5;         // initial guess — adjustTimeStep will correct it
writeControl    adjustableRunTime;
writeInterval   0.01;
```

**New:**
```
endTime         2000;

deltaT          1e-5;         // initial guess — adjustTimeStep will correct it
writeControl    adjustableRunTime;
writeInterval   40;            // same endTime/writeInterval ratio as before (~50 snapshots)
```

**Why:** The paper averages the periodic hill case over T_f = 200 flow-through
times (vs. T_f = 50 for their plain channel case — periodic hill needs more
because of slow, large-scale unsteadiness that plain near-wall turbulence
doesn't have: the reattachment point wanders over time, and the separated shear
layer flaps at low frequency, both with much longer correlation times than
ordinary near-wall eddies). One flow-through time here = L_x/U_b = 9/1 = 9 (time
units, this case's own convention). So the averaging window alone needs
200 × 9 = 1800. On top of that, ~22 flow-through times (200 time units) are
reserved as spin-up, to let the flow forget its initial condition and reach
statistical stationarity before averaging starts (matches `timeStart 200` in
`system/functions`, change #4 below) — so exactly 1800 of simulated time
actually gets averaged. Total: 200 + 1800 = 2000. The old value (0.5) is
roughly 0.06 flow-through times — nowhere near enough to develop the flow at
all, let alone average statistics over it.

`writeInterval` was scaled up in the same proportion as `endTime` (0.5/0.01 = 50
snapshots before → 2000/40 = 50 snapshots now), purely to keep disk usage
comparable; it doesn't affect the statistics themselves (see change #4 - the
averaging itself now happens every timestep, independent of how often results
get written to disk).

---

## 4. `mlSim/system/functions` and `spaldingSim/system/functions`

**Old:**
```
wallShearStress
{
    type            wallShearStress;
    libs            ("libfieldFunctionObjects.so");
    solver          incompressibleFluid;
    patches
    (
        hills
    );

    executeControl  writeTime;
    writeControl    writeTime;
}
```
(no `fieldAverage` object existed at all)

**New:**
```
wallShearStress
{
    type            wallShearStress;
    libs            ("libfieldFunctionObjects.so");
    solver          incompressibleFluid;
    patches
    (
        hills
    );

    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;
}

fieldAverage1
{
    type            fieldAverage;
    libs            ("libfieldFunctionObjects.so");
    timeStart       200;
    executeControl  timeStep;
    executeInterval 1;
    writeControl    writeTime;

    fields
    (
        U
        {
            mean        on;
            prime2Mean  on;
            base        time;
        }
        p
        {
            mean        on;
            prime2Mean  off;
            base        time;
        }
        wallShearStress
        {
            mean        on;
            prime2Mean  off;
            base        time;
        }
    );
}
```

**Why, two separate things:**

- **No time-averaging existed at all before.** The MATLAB post-processing
  scripts were reading raw, instantaneous `wallShearStress`/`U` at a single
  time (`0.1` or similar). For LES, an instantaneous field is one chaotic
  turbulent snapshot, not the converged mean quantity Spalding's law actually
  describes — comparing the two isn't a meaningful comparison, independent of
  whether the CFD itself is even right. `fieldAverage1` accumulates a proper
  running time-average of U, p, and wallShearStress starting at `timeStart 200`
  (the end of the spin-up period from change #3), so only the statistically
  stationary portion of the run contributes to the mean.
- **`wallShearStress`'s own `executeControl` had to change too** (from
  `writeTime` to every timestep). It was only recomputing a fresh value at each
  `writeTime` (previously every 0.01, now every 40 time units) — if
  `fieldAverage` executed every timestep while `wallShearStress` only updated
  every ~40 time units, `fieldAverage` would just be re-averaging the same
  stale value over and over between updates (a coarse, sparse-sample estimate
  wearing the clothes of a proper continuous average, not one). Computing it
  every timestep is cheap relative to the flow solve itself, and it still only
  *writes* the instantaneous field at `writeTime` (via `writeControl`, left
  unchanged), so disk usage for that field is the same as before.

---

## 5. `wallShearCompFine.m`

**Old:**
```matlab
time = '0.09';

cases(1).name  = 'Spalding';
cases(1).dir   = 'spaldingSim';
cases(1).time  ='0.1';
cases(1).color = [0.0000 0.4470 0.7410];

cases(2).name  = 'ML';
cases(2).dir   = 'mlSim';
cases(2).time  = '0.1';
cases(2).color = [0.8500 0.3250 0.0980];
```
```matlab
file = fullfile( ...
    cases(c).dir, ...
    cases(c).time, ...
    'wallShearStress');
```

**New:**
```matlab
cases(1).name  = 'Spalding';
cases(1).dir   = 'spaldingSim';
cases(1).time  = '2000';
cases(1).color = [0.0000 0.4470 0.7410];

cases(2).name  = 'ML';
cases(2).dir   = 'mlSim';
cases(2).time  = '2000';
cases(2).color = [0.8500 0.3250 0.0980];
```
```matlab
file = fullfile( ...
    cases(c).dir, ...
    cases(c).time, ...
    'wallShearStressMean');
```

**Why:** Matches change #3 (new `endTime`) and change #4 (the new mean field
that actually exists to read now). Was reading `wallShearStress` (instantaneous)
at `t≈0.1` (an early, barely-past-transient instant under the old, wrong
`endTime=0.5`); now reads `wallShearStressMean` (the properly time-averaged
field) at `t=2000` (the new, correct final/averaged time).

---

## 6. `velocityCompFine.m`

**Old:**
```matlab
time = '0.1';
```
```matlab
Ufile = fullfile( ...
    cases(c).dir, ...
    cases(c).time, ...
    'U');
```

**New:**
```matlab
time = '2000';
```
```matlab
Ufile = fullfile( ...
    cases(c).dir, ...
    cases(c).time, ...
    'UMean');
```

**Why:** Same reasoning as change #5, for the velocity-profile comparison
script instead of the skin-friction one. `C` (cell centres, a fixed geometric
field) is intentionally left reading the plain, non-averaged field — cell
positions don't fluctuate in time, so there's nothing to average.

---

## 7. `mlSim/periodicHillVelocityProfiles.m` and `spaldingSim/periodicHillVelocityProfiles.m`

**Old:**
```matlab
caseDir = '.';
timeName = '0.08';
```
```matlab
fprintf('Reading U...\n');

Ufile = fullfile(caseDir,timeName,'U');
```

**New:**
```matlab
caseDir = '.';
timeName = '2000';
```
```matlab
fprintf('Reading UMean...\n');

Ufile = fullfile(caseDir,timeName,'UMean');
```

**Why:** Same reasoning as changes #5/#6 — these are each case's own
single-case velocity-profile script (as opposed to `velocityCompFine.m`, which
compares both cases together), with the identical instantaneous-vs-time-and
stale-time issue, fixed the same way.

---

## Summary of files touched

```
mlSim/Allrun
mlSim/constant/transportProperties
mlSim/system/controlDict
mlSim/system/functions
mlSim/periodicHillVelocityProfiles.m
spaldingSim/Allrun
spaldingSim/constant/transportProperties
spaldingSim/system/controlDict
spaldingSim/system/functions
spaldingSim/periodicHillVelocityProfiles.m
wallShearCompFine.m
velocityCompFine.m
```

Not touched (already matched the paper): `system/blockMeshDict` in both cases
(domain 9×3.035×4.5 H, mesh 120×40×60 — both already exactly correct), and the
wall model's own C++ code (`turbModel.nu(patchi)` already reads the case's own
fluid properties correctly, not a hardcoded constant).
