const G0 = 9.80665;

// Optimal delta-v split via Lagrange multipliers.
//
// The optimality condition d/d(Δv_i) log(1 + k_i) = λ gives
// MR_i = (1 - 1/(λ·Isp_i·g₀)) / s_i per stage, so the full
// solve reduces to bisecting on λ until Σ Δv_i = Δv_total.
export function optimize(payload, totalDeltaV, stages) {
  const isps = stages.map((s) => s.specificImpulse);
  const pmfs = stages.map((s) => s.propellantMassFraction);
  const sfs = pmfs.map((pmf) => 1 - pmf);

  function splitForLambda(lambda) {
    return isps.map((isp, i) => {
      const mr = (1 - 1 / (lambda * isp * G0)) / sfs[i];
      return mr > 1 ? isp * G0 * Math.log(mr) : 0;
    });
  }

  function dvSum(lambda) {
    return splitForLambda(lambda).reduce((a, b) => a + b, 0);
  }

  let lo = 1e-10;
  let hi = 1;
  while (dvSum(hi) < totalDeltaV) hi *= 2;

  for (let i = 0; i < 200; i++) {
    const mid = (lo + hi) / 2;
    const sum = dvSum(mid);
    if (Math.abs(sum - totalDeltaV) < 1e-10) { lo = mid; break; }
    if (sum < totalDeltaV) lo = mid;
    else hi = mid;
  }

  const dvs = splitForLambda(lo);
  const result = [];
  let carry = payload;

  for (let i = stages.length - 1; i >= 0; i--) {
    const mr = Math.exp(dvs[i] / (isps[i] * G0));
    const wetMass = dvs[i] > 1e-6 ? carry * (mr - 1) / (1 - mr * sfs[i]) : 0;
    const propellantMass = pmfs[i] * wetMass;

    result[i] = {
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

  return { payload, totalDeltaV, stages: result, totalMass: carry };
}
