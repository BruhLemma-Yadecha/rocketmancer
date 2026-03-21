import { color, fmt } from '../palette';

const PROPERTIES = [
  { name: '\u0394v', key: 'deltaV', unit: 'm/s' },
  { name: 'Mass Ratio', key: 'massRatio', unit: null },
  { name: 'Payload Mass', key: 'payloadMass', unit: 'kg' },
  { name: 'Wet Mass', key: 'totalMass', unit: 'kg' },
  { name: 'Dry Mass', key: 'structuralMass', unit: 'kg' },
  { name: 'Propellant Mass', key: 'propellantMass', unit: 'kg' },
];

const STATS = [
  { label: '\u0394v', value: r => r.totalDeltaV.toFixed(0), unit: 'm/s' },
  { label: 'Payload', value: r => r.payload.toFixed(0), unit: 'kg' },
  { label: 'Total Mass', value: r => r.totalMass.toFixed(2), unit: 'kg' },
];

export default function Results({ result, name }) {
  return (
    <div className="card fade-in results-card">
      <div className="card-title">Results</div>
      {result ? (
        <>
          <div className="result-header">
            <span className="rocket-name">{name}</span>
            <div className="stats-bar">
              {STATS.map((stat, i) => (
                <div className="stat" key={stat.label}>
                  <div className={`stat-label color-${color(i)}`}>{stat.label}</div>
                  <div className="stat-value">
                    {stat.value(result)}
                    {stat.unit && <span className="stat-unit">{stat.unit}</span>}
                  </div>
                </div>
              ))}
            </div>
          </div>

          <table>
            <thead>
              <tr>
                <th className="color-dim">Stage</th>
                {PROPERTIES.map((prop, i) => (
                  <th key={prop.key}>
                    <span className={`color-${color(i)}`}>{prop.name}</span>
                    {prop.unit && <span className="prop-unit"> ({prop.unit})</span>}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {result.stages.map((stage, i) => (
                <tr key={i}>
                  <td>
                    <span className={`stage-badge bg-${color(i)}`}>{i + 1}</span>
                  </td>
                  {PROPERTIES.map(prop => (
                    <td key={prop.key}>
                      <span className="value-pill">{fmt(stage[prop.key])}</span>
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        </>
      ) : (
        <div style={{ textAlign: 'center', padding: '4rem 0', color: 'var(--text-muted)' }}>
          Enter valid parameters to see results
        </div>
      )}
    </div>
  );
}
