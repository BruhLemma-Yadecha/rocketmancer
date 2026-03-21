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
    const result = optimize(500, 3000, [
      { specificImpulse: 350, propellantMassFraction: 0.9 },
    ]);

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
      { specificImpulse: 400, propellantMassFraction: 0.90 },
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
