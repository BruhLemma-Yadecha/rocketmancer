export function SettingsToggle({ open, setOpen }) {
  return (
    <button
      className="settings-toggle"
      onClick={() => setOpen(!open)}
      aria-label={open ? 'Close settings' : 'Open settings'}
    >
      <svg width="18" height="18" viewBox="0 0 20 20" fill="currentColor">
        <path
          fillRule="evenodd"
          d="M11.49 3.17c-.38-1.56-2.6-1.56-2.98 0a1.532 1.532 0 01-2.286.948c-1.372-.836-2.942.734-2.106 2.106.54.886.062 2.042-.947 2.287-1.561.379-1.561 2.6 0 2.978a1.532 1.532 0 01.947 2.287c-.836 1.372.734 2.942 2.106 2.106a1.532 1.532 0 012.287.947c.379 1.561 2.6 1.561 2.978 0a1.533 1.533 0 012.287-.947c1.372.836 2.942-.734 2.106-2.106a1.533 1.533 0 01.947-2.287c1.561-.379 1.561-2.6 0-2.978a1.532 1.532 0 01-.947-2.287c.836-1.372-.734-2.942-2.106-2.106a1.532 1.532 0 01-2.287-.947zM10 13a3 3 0 100-6 3 3 0 000 6z"
          clipRule="evenodd"
        />
      </svg>
    </button>
  );
}

function Setting({ name, description, unit, children }) {
  return (
    <div className="setting-row">
      <div className="setting-info">
        <span className="setting-name">{name}</span>
        <span className="setting-desc">{description}</span>
      </div>
      <div className="setting-input">
        {children}
        {unit && <span className="setting-unit">{unit}</span>}
      </div>
    </div>
  );
}

export function SettingsPanel({ minContribution, setMinContribution, stageCount }) {
  const maxContribution = Math.floor(100 / stageCount);

  return (
    <div className="card settings-card">
      <div className="card-title">Settings</div>

      <Setting
        name="Minimum stage contribution"
        description={`Force each stage to provide at least this percentage of total delta-v (max ${maxContribution}% with ${stageCount} stages). Useful for rockets with mixed propellant types where the optimizer would otherwise assign zero to lower-Isp stages. Setting this above zero constrains the solver and may produce suboptimal results in the global problem space.`}
        unit="%"
      >
        <input
          type="number"
          min="0"
          max={maxContribution}
          step="1"
          value={minContribution}
          onChange={e =>
            setMinContribution(
              Math.max(0, Math.min(maxContribution, parseFloat(e.target.value) || 0))
            )
          }
        />
      </Setting>
    </div>
  );
}
