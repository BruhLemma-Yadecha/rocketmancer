import { describe, it, expect } from 'vitest';
import { optimize } from './solver';

describe('optimize', () => {
  it('2-stage: matches hand-verified reference values', () => {
    const result = optimize(1000, 4760.08, [
      { specificImpulse: 307.36, propellantMassFraction: 0.83 },
      { specificImpulse: 348.81, propellantMassFraction: 0.87 },
    ]);

    const totalDv = result.stages.reduce((sum, s) => sum + s.deltaV, 0);
    expect(totalDv).toBeCloseTo(4760.08, 2);
    expect(result.totalMass).toBeCloseTo(6441.86, 0);

    expect(result.stages[0].deltaV).toBeCloseTo(1344.57, 0);
    expect(result.stages[0].totalMass).toBeCloseTo(2793.05, 0);
    expect(result.stages[1].deltaV).toBeCloseTo(3415.51, 0);
    expect(result.stages[1].totalMass).toBeCloseTo(2648.81, 0);

    expect(result.stages[1].payloadMass).toBeCloseTo(1000, 2);
    expect(result.stages[0].payloadMass).toBeCloseTo(1000 + result.stages[1].totalMass, 2);
  });

  it('single stage gets all the delta-v', () => {
    const result = optimize(500, 3000, [{ specificImpulse: 350, propellantMassFraction: 0.9 }]);

    expect(result.stages).toHaveLength(1);
    expect(result.stages[0].deltaV).toBeCloseTo(3000, 2);
    expect(result.stages[0].payloadMass).toBeCloseTo(500, 2);
    expect(result.totalMass).toBeGreaterThan(500);
  });

  it('favours the higher-Isp stage', () => {
    const result = optimize(1000, 5000, [
      { specificImpulse: 250, propellantMassFraction: 0.85 },
      { specificImpulse: 450, propellantMassFraction: 0.85 },
    ]);

    expect(result.stages[1].deltaV).toBeGreaterThan(result.stages[0].deltaV);
  });

  it('total mass = payload + sum of stage wet masses', () => {
    const result = optimize(800, 6000, [
      { specificImpulse: 300, propellantMassFraction: 0.85 },
      { specificImpulse: 350, propellantMassFraction: 0.88 },
      { specificImpulse: 400, propellantMassFraction: 0.9 },
    ]);

    const stagesMass = result.stages.reduce((sum, s) => sum + s.totalMass, 0);
    expect(result.totalMass).toBeCloseTo(result.payload + stagesMass, 4);
  });

  it('propellant + structural = wet mass per stage', () => {
    const result = optimize(1000, 4760.08, [
      { specificImpulse: 307.36, propellantMassFraction: 0.83 },
      { specificImpulse: 348.81, propellantMassFraction: 0.87 },
    ]);

    for (const stage of result.stages) {
      expect(stage.propellantMass + stage.structuralMass).toBeCloseTo(stage.totalMass, 4);
      expect(stage.propellantMass).toBeCloseTo(stage.propellantMassFraction * stage.totalMass, 4);
    }
  });
});

describe('optimize with minContribution', () => {
  const twoStage = [
    { specificImpulse: 307.36, propellantMassFraction: 0.83 },
    { specificImpulse: 348.81, propellantMassFraction: 0.87 },
  ];

  // Saturn V-like: stage 1 has much lower Isp, unconstrained gives it 0
  const mixedIsp = [
    { specificImpulse: 304, propellantMassFraction: 0.94 },
    { specificImpulse: 421, propellantMassFraction: 0.92 },
    { specificImpulse: 421, propellantMassFraction: 0.89 },
  ];

  it('minContribution=0 matches unconstrained result', () => {
    const unconstrained = optimize(1000, 4760.08, twoStage);
    const withZero = optimize(1000, 4760.08, twoStage, 0);

    for (let i = 0; i < twoStage.length; i++) {
      expect(withZero.stages[i].deltaV).toBeCloseTo(unconstrained.stages[i].deltaV, 8);
      expect(withZero.stages[i].totalMass).toBeCloseTo(unconstrained.stages[i].totalMass, 8);
    }
  });

  it('non-binding constraint does not change results', () => {
    // Both stages get >28% unconstrained, so 5% floor should be non-binding
    const unconstrained = optimize(1000, 4760.08, twoStage);
    const withFloor = optimize(1000, 4760.08, twoStage, 5);

    for (let i = 0; i < twoStage.length; i++) {
      expect(withFloor.stages[i].deltaV).toBeCloseTo(unconstrained.stages[i].deltaV, 8);
    }
  });

  it('pins low-Isp stage to minimum when unconstrained gives it zero', () => {
    const result = optimize(140000, 9400, mixedIsp, 10);

    // Stage 1 should get at least 10% of total delta-v
    expect(result.stages[0].deltaV).toBeGreaterThanOrEqual(9400 * 0.1 - 1);
    // All stages should contribute
    for (const stage of result.stages) {
      expect(stage.deltaV).toBeGreaterThan(0);
    }
    // Delta-v should still sum to total
    const totalDv = result.stages.reduce((sum, s) => sum + s.deltaV, 0);
    expect(totalDv).toBeCloseTo(9400, 2);
  });

  it('mass invariants hold with constraints active', () => {
    const result = optimize(140000, 9400, mixedIsp, 15);

    // Total mass = payload + sum of stage masses
    const stagesMass = result.stages.reduce((sum, s) => sum + s.totalMass, 0);
    expect(result.totalMass).toBeCloseTo(result.payload + stagesMass, 4);

    // Per-stage: propellant + structural = wet mass
    for (const stage of result.stages) {
      expect(stage.propellantMass + stage.structuralMass).toBeCloseTo(stage.totalMass, 4);
    }
  });

  it('constrained result has higher total mass than unconstrained', () => {
    const unconstrained = optimize(140000, 9400, mixedIsp, 0);
    const constrained = optimize(140000, 9400, mixedIsp, 10);

    expect(constrained.totalMass).toBeGreaterThanOrEqual(unconstrained.totalMass);
  });
});

describe('optimize validation', () => {
  const validStage = [{ specificImpulse: 300, propellantMassFraction: 0.85 }];

  it('rejects non-positive payload', () => {
    expect(() => optimize(0, 3000, validStage)).toThrow('Payload');
    expect(() => optimize(-1, 3000, validStage)).toThrow('Payload');
  });

  it('rejects non-positive delta-v', () => {
    expect(() => optimize(1000, 0, validStage)).toThrow('delta-v');
    expect(() => optimize(1000, -1, validStage)).toThrow('delta-v');
  });

  it('rejects invalid stage parameters', () => {
    expect(() =>
      optimize(1000, 3000, [{ specificImpulse: 0, propellantMassFraction: 0.85 }])
    ).toThrow('specific impulse');
    expect(() =>
      optimize(1000, 3000, [{ specificImpulse: 300, propellantMassFraction: 1 }])
    ).toThrow('propellant mass fraction');
    expect(() =>
      optimize(1000, 3000, [{ specificImpulse: 300, propellantMassFraction: 0 }])
    ).toThrow('propellant mass fraction');
  });

  it('rejects empty stages array', () => {
    expect(() => optimize(1000, 3000, [])).toThrow('At least one stage');
  });

  it('rejects out-of-range minContribution', () => {
    expect(() => optimize(1000, 3000, validStage, -1)).toThrow('Minimum contribution');
    // 1 stage → max is 100%, but 2 stages → max is 50%
    const twoStages = [
      { specificImpulse: 300, propellantMassFraction: 0.85 },
      { specificImpulse: 350, propellantMassFraction: 0.85 },
    ];
    expect(() => optimize(1000, 3000, twoStages, 51)).toThrow('Minimum contribution');
  });
});
