import { useRef, useEffect, useState, useCallback } from 'react';
import { color, fmt } from '../palette';

const PROPERTIES = [
  { name: '\u0394v', key: 'deltaV', unit: 'm/s', decimals: 2 },
  { name: 'Mass Ratio', key: 'massRatio', unit: null, decimals: 4 },
  { name: 'Payload Mass', key: 'payloadMass', unit: 'kg', decimals: 2 },
  { name: 'Wet Mass', key: 'totalMass', unit: 'kg', decimals: 2 },
  { name: 'Dry Mass', key: 'structuralMass', unit: 'kg', decimals: 2 },
  { name: 'Propellant Mass', key: 'propellantMass', unit: 'kg', decimals: 2 },
];

const STATS = [
  { label: '\u0394v', value: r => fmt(r.totalDeltaV, 0), unit: 'm/s' },
  { label: 'Payload', value: r => fmt(r.payload, 0), unit: 'kg' },
  { label: 'Total Mass', value: r => fmt(r.totalMass, 2), unit: 'kg' },
];

function CopyPill({ value, formatted }) {
  const [copied, setCopied] = useState(false);

  const handleClick = useCallback(() => {
    navigator.clipboard.writeText(String(value)).then(() => {
      setCopied(true);
      setTimeout(() => setCopied(false), 1000);
    });
  }, [value]);

  return (
    <span
      className={`value-pill copyable${copied ? ' copied' : ''}`}
      onClick={handleClick}
      title="Click to copy"
    >
      {copied ? 'Copied' : formatted}
    </span>
  );
}

export default function Results({ result, name }) {
  const tableRef = useRef(null);

  useEffect(() => {
    if (!tableRef.current) return;
    const pills = tableRef.current.querySelectorAll('.value-pill');
    pills.forEach(p => (p.style.minWidth = ''));

    const cols = PROPERTIES.length;
    for (let c = 0; c < cols; c++) {
      let max = 0;
      for (let r = 0; r < result?.stages.length; r++) {
        const pill = pills[r * cols + c];
        if (pill) max = Math.max(max, pill.scrollWidth);
      }
      for (let r = 0; r < result?.stages.length; r++) {
        const pill = pills[r * cols + c];
        if (pill) pill.style.minWidth = `${max}px`;
      }
    }
  }, [result]);

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

          <table ref={tableRef}>
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
                      <CopyPill
                        value={stage[prop.key]}
                        formatted={fmt(stage[prop.key], prop.decimals)}
                      />
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
