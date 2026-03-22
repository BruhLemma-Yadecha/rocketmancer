import { color } from '../palette';
import FormattedInput from './FormattedInput';

export default function StagesInput({ stages, updateStage }) {
  return (
    <div className="card fade-in">
      <div className="card-title">Stages</div>

      <table>
        <thead>
          <tr>
            <th></th>
            <th className={`color-${color(0)}`}>Specific Impulse (s)</th>
            <th className={`color-${color(1)}`}>Mass Fraction</th>
          </tr>
        </thead>
        <tbody>
          {stages.map((stage, i) => (
            <tr key={i}>
              <td>
                <span className={`stage-badge bg-${color(i)}`}>{i + 1}</span>
              </td>
              <td>
                <FormattedInput
                  min="0"
                  step="0.01"
                  decimals={2}
                  value={stage.specificImpulse}
                  onChange={v => updateStage(i, 'specificImpulse', v)}
                />
              </td>
              <td>
                <FormattedInput
                  min="0"
                  max="1"
                  step="0.01"
                  decimals={3}
                  value={stage.propellantMassFraction}
                  onChange={v => updateStage(i, 'propellantMassFraction', v)}
                />
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
