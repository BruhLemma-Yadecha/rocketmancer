import { color } from '../palette';

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
                <input
                  type="number"
                  min="0"
                  step="0.01"
                  value={stage.specificImpulse}
                  onChange={e => updateStage(i, 'specificImpulse', parseFloat(e.target.value) || 0)}
                />
              </td>
              <td>
                <input
                  type="number"
                  min="0"
                  max="1"
                  step="0.01"
                  value={stage.propellantMassFraction}
                  onChange={e =>
                    updateStage(i, 'propellantMassFraction', parseFloat(e.target.value) || 0)
                  }
                />
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
