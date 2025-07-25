import DisplayProperty from './DisplayProperty';

const MASS_UNIT = 't';
const VELOCITY_UNIT = 'm/s';
const SPECIFIC_IMPULSE_UNIT = 's';
const DECIMAL_PLACES = 4;

const Rocket = ({ rocket, rocketName }) => {
  if (!rocket) {
    return (
      <div className="acrylic rounded-xl p-8 animate-fade-in">
        <div className="flex items-center justify-center h-64">
          <div className="text-center">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-400 mx-auto mb-4"></div>
            <p className="text-white/80">Calculating rocket parameters...</p>
          </div>
        </div>
      </div>
    );
  }

  const properties = [
    { name: 'Delta V', type: 'velocity' },
    { name: 'Mass Ratio', type: 'ratio' },
    { name: 'Payload Mass', type: 'mass' },
    { name: 'Wet Mass', type: 'mass' },
    { name: 'Dry Mass', type: 'mass' },
    { name: 'Structural Mass', type: 'mass' },
    { name: 'Propellant Mass', type: 'mass' },
    { name: 'Exhaust Velocity', type: 'velocity' },
    { name: 'Specific Impulse', type: 'time' },
    { name: 'Propellant Mass Fraction', type: 'percentage' },
  ];

  return (
    <div className="card-glass animate-fade-in">
      {/* Compact Header with stats inline */}
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center">
          <div className="inline-flex items-center justify-center w-7 h-7 bg-white/20 rounded-full mr-2 animate-rocket-pulse">
            <span className="text-sm">🚀</span>
          </div>
          <div>
            <h2 className="text-lg font-bold text-white">{rocketName}</h2>
            <p className="text-xs text-white/70">Optimized configuration</p>
          </div>
        </div>

        {/* Key Stats - Inline to the right */}
        <div className="flex gap-3">
          <div className="glass rounded-lg px-3 py-2 text-center min-w-[80px]">
            <div className="text-xs font-medium text-blue-300">Stages</div>
            <div className="text-lg font-bold text-white">{rocket.totalStages}</div>
          </div>
          <div className="glass rounded-lg px-3 py-2 text-center min-w-[90px]">
            <div className="text-xs font-medium text-green-300">ΔV</div>
            <div className="text-sm font-bold text-white">
              {rocket.totalDeltaV.toFixed(0)}
              <span className="text-xs font-normal block text-white/70">m/s</span>
            </div>
          </div>
          <div className="glass rounded-lg px-3 py-2 text-center min-w-[90px]">
            <div className="text-xs font-medium text-purple-300">Payload</div>
            <div className="text-sm font-bold text-white">
              {rocket.payload.toFixed(0)}
              <span className="text-xs font-normal block text-white/70">kg</span>
            </div>
          </div>
        </div>
      </div>

      {/* Stage Details Table - More compact */}
      <div className="overflow-x-auto">
        <table className="w-full text-sm">
          <thead>
            <tr className="border-b border-white/20">
              <th className="text-left py-2 px-2 font-medium text-white text-xs">Property</th>
              {rocket.stages.map((stage, index) => (
                <th
                  key={index}
                  className="text-center py-2 px-1 font-medium text-white bg-white/5 rounded-t-lg text-xs"
                >
                  <div className="flex items-center justify-center">
                    <span className="w-1 h-1 bg-blue-400 rounded-full mr-1"></span>S{index + 1}
                  </div>
                </th>
              ))}
            </tr>
          </thead>
          <tbody className="divide-y divide-white/10">
            {properties.map((property, index) => (
              <DisplayProperty
                key={index}
                name={property.name}
                type={property.type}
                stages={rocket.stages}
              />
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
};

export default Rocket;
