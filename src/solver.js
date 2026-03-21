const G0 = 9.80665;
const BISECT_TOL = 1e-10;
const BISECT_MAX_ITER = 200;
const BRACKET_MAX_ITER = 100;
const DV_EPSILON = 1e-6;

// Bisection search for λ such that f(λ) = target.
// f must be monotonically increasing in λ.
function bisect(f, target, lo = BISECT_TOL, hi = 1) {
  let bracket = 0;
  while (f(hi) < target) {
    hi *= 2;
    if (++bracket > BRACKET_MAX_ITER) return { lambda: hi, converged: false };
  }

  for (let i = 0; i < BISECT_MAX_ITER; i++) {
    const mid = (lo + hi) / 2;
    const val = f(mid);
    if (Math.abs(val - target) < BISECT_TOL) return { lambda: mid, converged: true };
    if (val < target) lo = mid;
    else hi = mid;
  }

  return { lambda: (lo + hi) / 2, converged: false };
}

// Optimal delta-v split for the given stages via Lagrange multipliers.
//
// The optimality condition ∂/∂(Δv_i) log(1 + k_i) = λ gives
// MR_i = (1 - 1/(λ·Isp_i·g₀)) / s_i per stage, so the full
// solve reduces to bisecting on λ until Σ Δv_i = Δv_total.
//
// Returns { dvs: number[], lambda: number, converged: boolean }.
function maxDeltaV(isps, sfs) {
  return isps.reduce((sum, isp, i) => sum + isp * G0 * Math.log(1 / sfs[i]), 0);
}

function solveUnconstrainedSplit(isps, sfs, totalDv) {
  if (totalDv > maxDeltaV(isps, sfs)) {
    throw new Error('Requested delta-v exceeds what these stages can physically provide');
  }

  function dvForLambda(lambda) {
    return isps.map((isp, i) => {
      const mr = (1 - 1 / (lambda * isp * G0)) / sfs[i];
      return mr > 1 ? isp * G0 * Math.log(mr) : 0;
    });
  }

  function dvSum(lambda) {
    return dvForLambda(lambda).reduce((a, b) => a + b, 0);
  }

  const { lambda, converged } = bisect(dvSum, totalDv);
  return { dvs: dvForLambda(lambda), lambda, converged };
}

// Marginal cost of stage i at a given delta-v:
// ∂/∂(Δv_i) log(1+k_i) = 1 / (Isp_i · g₀ · (1 - MR_i · s_i))
function marginalCost(isp, sf, dv) {
  const mr = Math.exp(dv / (isp * G0));
  return 1 / (isp * G0 * (1 - mr * sf));
}

// Compute per-stage masses from a delta-v split, working top-down.
function computeStageMasses(dvs, isps, pmfs, sfs, payload) {
  const stages = [];
  let carry = payload;

  for (let i = dvs.length - 1; i >= 0; i--) {
    const mr = Math.exp(dvs[i] / (isps[i] * G0));
    const wetMass = dvs[i] > DV_EPSILON ? (carry * (mr - 1)) / (1 - mr * sfs[i]) : 0;
    const propellantMass = pmfs[i] * wetMass;

    stages[i] = {
      specificImpulse: isps[i],
      propellantMassFraction: pmfs[i],
      deltaV: dvs[i],
      massRatio: mr,
      payloadMass: carry,
      totalMass: wetMass,
      propellantMass,
      structuralMass: wetMass - propellantMass,
    };

    carry += wetMass;
  }

  return { stages, totalMass: carry };
}

// Active-set method: pin stages that fall below the minimum,
// then re-optimize the remainder across unpinned stages.
// After each solve, also check dual variables (μ_i) for pinned
// stages so that we can release any where μ_i < 0 (overconstrained).
function solveConstrainedSplit(isps, sfs, totalDv, minDv) {
  const n = isps.length;
  const pinned = new Array(n).fill(false);
  let dvs;
  let converged = true;

  for (let iter = 0; iter < 2 * n; iter++) {
    const pinnedDv = pinned.reduce((sum, p) => sum + (p ? minDv : 0), 0);
    const freeDv = totalDv - pinnedDv;
    const freeIndices = isps.map((_, i) => i).filter(i => !pinned[i]);

    if (freeIndices.length === 0) {
      dvs = new Array(n).fill(totalDv / n);
      break;
    }

    const freeIsps = freeIndices.map(i => isps[i]);
    const freeSfs = freeIndices.map(i => sfs[i]);

    const solution = solveUnconstrainedSplit(freeIsps, freeSfs, freeDv);
    if (!solution.converged) converged = false;

    dvs = new Array(n);
    freeIndices.forEach((idx, j) => {
      dvs[idx] = solution.dvs[j];
    });
    for (let i = 0; i < n; i++) {
      if (pinned[i]) dvs[i] = minDv;
    }

    let changed = false;

    // Pin free stages that fell below the minimum
    for (const i of freeIndices) {
      if (dvs[i] < minDv) {
        pinned[i] = true;
        changed = true;
      }
    }

    // Release pinned stages with negative dual variable:
    // μ_i = marginalCost(i, minDv) - λ_free
    // μ_i < 0 means this stage is cheaper than the free stages
    // at the boundary, so it should get more than minDv.
    for (let i = 0; i < n; i++) {
      if (pinned[i] && marginalCost(isps[i], sfs[i], minDv) < solution.lambda) {
        pinned[i] = false;
        changed = true;
      }
    }

    if (!changed) break;
  }

  return { dvs, converged };
}

function validate(payload, totalDeltaV, stages, minContribution) {
  if (!Number.isFinite(payload) || payload <= 0) {
    throw new Error('Payload must be a positive number');
  }
  if (!Number.isFinite(totalDeltaV) || totalDeltaV <= 0) {
    throw new Error('Total delta-v must be a positive number');
  }
  if (!Array.isArray(stages) || stages.length === 0) {
    throw new Error('At least one stage is required');
  }
  for (let i = 0; i < stages.length; i++) {
    const { specificImpulse, propellantMassFraction } = stages[i];
    if (!Number.isFinite(specificImpulse) || specificImpulse <= 0) {
      throw new Error(`Stage ${i + 1}: specific impulse must be a positive number`);
    }
    if (
      !Number.isFinite(propellantMassFraction) ||
      propellantMassFraction <= 0 ||
      propellantMassFraction >= 1
    ) {
      throw new Error(
        `Stage ${i + 1}: propellant mass fraction must be between 0 and 1 (exclusive)`
      );
    }
  }
  const maxContribution = 100 / stages.length;
  if (
    !Number.isFinite(minContribution) ||
    minContribution < 0 ||
    minContribution > maxContribution
  ) {
    throw new Error(
      `Minimum contribution must be between 0 and ${Math.floor(maxContribution)} for ${stages.length} stages`
    );
  }
}

export function optimize(payload, totalDeltaV, stages, minContribution = 0) {
  validate(payload, totalDeltaV, stages, minContribution);

  const isps = stages.map(s => s.specificImpulse);
  const pmfs = stages.map(s => s.propellantMassFraction);
  const sfs = pmfs.map(pmf => 1 - pmf);
  const minDv = (minContribution / 100) * totalDeltaV;

  const { dvs, converged } =
    minDv > 0
      ? solveConstrainedSplit(isps, sfs, totalDeltaV, minDv)
      : solveUnconstrainedSplit(isps, sfs, totalDeltaV);

  if (!converged) {
    throw new Error('Solver failed to converge');
  }

  const { stages: result, totalMass } = computeStageMasses(
    minDv > 0 ? dvs : dvs,
    isps,
    pmfs,
    sfs,
    payload
  );

  return { payload, totalDeltaV, stages: result, totalMass };
}
