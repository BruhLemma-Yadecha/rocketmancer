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
    <div className="space-y-6 animate-fade-in">
      {/* Header Card */}
      <div className="card-glass">
        <div className="text-center mb-6">
          <div className="inline-flex items-center justify-center w-12 h-12 bg-white/20 rounded-full mb-3 animate-rocket-pulse">
            <span className="text-2xl">🚀</span>
          </div>
          <h2 className="text-2xl font-bold text-white mb-2">{rocketName}</h2>
          <p className="text-white/70">Optimized rocket configuration</p>
        </div>
        
        {/* Key Stats */}
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
          <div className="glass rounded-lg p-4 text-center hover:glass-strong transition-all duration-300">
            <div className="text-sm font-medium text-blue-300 mb-1">Stages</div>
            <div className="text-2xl font-bold text-white">{rocket.totalStages}</div>
          </div>
          <div className="glass rounded-lg p-4 text-center hover:glass-strong transition-all duration-300">
            <div className="text-sm font-medium text-green-300 mb-1">Total ΔV</div>
            <div className="text-xl font-bold text-white">
              {rocket.totalDeltaV.toFixed(DECIMAL_PLACES)}
              <span className="text-sm font-normal ml-1 text-white/70">{VELOCITY_UNIT}</span>
            </div>
          </div>
          <div className="glass rounded-lg p-4 text-center hover:glass-strong transition-all duration-300">
            <div className="text-sm font-medium text-purple-300 mb-1">Payload</div>
            <div className="text-xl font-bold text-white">
              {rocket.payload.toFixed(DECIMAL_PLACES)}
              <span className="text-sm font-normal ml-1 text-white/70">kg</span>
            </div>
          </div>
          <div className="glass rounded-lg p-4 text-center hover:glass-strong transition-all duration-300">
            <div className="text-sm font-medium text-orange-300 mb-1">Total Mass</div>
            <div className="text-xl font-bold text-white">
              {rocket.totalMass.toFixed(DECIMAL_PLACES)}
              <span className="text-sm font-normal ml-1 text-white/70">kg</span>
            </div>
          </div>
        </div>
      </div>

      {/* Detailed Properties Table */}
      <div className="card">
        <h3 className="text-lg font-semibold text-white mb-6 flex items-center">
          <div className="w-2 h-2 bg-blue-400 rounded-full mr-3 animate-pulse"></div>
          Stage Details
        </h3>
        
        <div className="overflow-x-auto">
          <table className="w-full">
            <thead>
              <tr className="border-b border-white/20">
                <th className="text-left py-3 px-4 font-semibold text-white">Property</th>
                {rocket.stages.map((stage, index) => (
                  <th key={index} className="text-center py-3 px-4 font-semibold text-white bg-white/5 rounded-t-lg">
                    <div className="flex items-center justify-center">
                      <span className="w-2 h-2 bg-blue-400 rounded-full mr-2"></span>
                      Stage {index + 1}
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
    </div>
  );
};

export default Rocket;
