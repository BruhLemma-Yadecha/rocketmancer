import { color } from '../palette';

const PARAMS = [
  { label: 'Name', key: 'name', type: 'text' },
  { label: 'Stages', key: 'stageCount', type: 'number', min: 1, max: 10, step: 1 },
  { label: '\u0394v (m/s)', key: 'totalDeltaV', type: 'number', min: 0, step: 0.01 },
  { label: 'Payload (kg)', key: 'payload', type: 'number', min: 0, step: 0.01 },
];

export default function Parameters({ name, setName, payload, setPayload, totalDeltaV, setTotalDeltaV, stageCount, setStageCount }) {
  const values = { name, stageCount, totalDeltaV, payload };
  const setters = {
    name: setName,
    stageCount: setStageCount,
    totalDeltaV: setTotalDeltaV,
    payload: setPayload,
  };

  return (
    <div className="card fade-in">
      <div className="card-title">Parameters</div>

      <div className="param-grid">
        {PARAMS.map((param, i) => (
          <div className="param-cell" key={param.key}>
            <label className={`color-${color(i)}`}>{param.label}</label>
            <input
              type={param.type}
              min={param.min}
              max={param.max}
              step={param.step}
              value={values[param.key]}
              onChange={(e) =>
                setters[param.key](param.type === 'number' ? parseFloat(e.target.value) || 0 : e.target.value)
              }
            />
          </div>
        ))}
      </div>
    </div>
  );
}
