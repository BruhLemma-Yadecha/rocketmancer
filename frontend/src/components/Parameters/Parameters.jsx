import { useState, useEffect } from 'react';
import RocketService from '../../services/rocketService';
import ParametersStage from './ParametersStage';

const Parameters = ({ setRocket, rocketName, setRocketName }) => {
  const [config, setConfig] = useState(undefined);
  const [refresh, setRefresh] = useState(false);
  const [isOptimizing, setIsOptimizing] = useState(false);

  const loadDefaultConfiguration = () => {
    const defaultConfig = {
      name: "Billy Jean",
      totalStages: 2,
      totalDeltaV: 4760.08,
      payload: 1000.0,
      stages: [
        {
          specificImpulse: 307.36,
          propellantMassFraction: 0.83
        },
        {
          specificImpulse: 348.81,
          propellantMassFraction: 0.87
        }
      ]
    };
    setConfig(defaultConfig);
    setRocketName(defaultConfig.name);
  };

  useEffect(() => {
    // Load default configuration on component mount
    if (!config) {
      loadDefaultConfiguration();
      return;
    }

    // Use RocketService to optimize the rocket
    const optimizeRocket = async () => {
      try {
        setIsOptimizing(true);
        const optimizedRocket = await RocketService.optimize(config);
        setRocket(optimizedRocket);
      } catch (error) {
        console.error('Failed to optimize rocket:', error);
      } finally {
        setIsOptimizing(false);
      }
    };

    optimizeRocket();
  }, [refresh, config]);

  const setTotalStages = totalStages => {
    totalStages = parseInt(totalStages);
    if (totalStages < 1) return;
    if (totalStages < config.totalStages) {
      const newStages = config.stages.slice(0, totalStages);
      setConfig({ ...config, totalStages, stages: newStages });
    } else if (totalStages > config.totalStages) {
      // add new dummy rows
      const dummyStage = { specificImpulse: 300.0, propellantMassFraction: 0.9 };
      const newStages = [...config.stages];
      for (let i = config.totalStages; i < totalStages; i++) {
        newStages.push(dummyStage);
      }
      setConfig({ ...config, totalStages, stages: newStages, name: rocketName });
    }
  };

  const setTotalDeltaV = totalDeltaV => {
    totalDeltaV = parseFloat(totalDeltaV);
    if (totalDeltaV < 0) return;
    setConfig({ ...config, totalDeltaV, name: rocketName });
  };

  const setStages = stages => {
    setConfig({ ...config, stages, name: rocketName });
  };

  const setPayload = payload => {
    payload = parseFloat(payload);
    if (payload <= 0) return;
    setConfig({ ...config, payload, name: rocketName });
  };

  const setName = editedName => {
    setRocketName(editedName);
  };

  if (!config) {
    return (
      <div className="card animate-fade-in">
        <div className="flex items-center justify-center h-32">
          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600"></div>
        </div>
      </div>
    );
  }

  return (
    <div className="space-y-6 animate-fade-in">
      {/* Header */}
      <div className="card-glass">
        <h2 className="text-2xl font-bold text-white mb-6 flex items-center">
          <div className="w-2 h-2 bg-blue-400 rounded-full mr-3 animate-pulse"></div>
          Rocket Parameters
        </h2>
        
        {/* Main Configuration */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {/* Rocket Name */}
          <div>
            <label className="block text-sm font-medium text-white mb-2 flex items-center">
              <span className="mr-2">🚀</span>
              Rocket Name
            </label>
            <input
              type="text"
              value={rocketName}
              onChange={e => setName(e.target.value)}
              className="form-input"
              placeholder="Enter rocket name..."
            />
          </div>

          {/* Total Stages */}
          <div>
            <label className="block text-sm font-medium text-white mb-2 flex items-center">
              <span className="mr-2">🔢</span>
              Total Stages
            </label>
            <input
              type="number"
              min="1"
              max="10"
              value={config.totalStages}
              onChange={e => setTotalStages(e.target.value)}
              className="form-input"
            />
          </div>

          {/* Total Delta-V */}
          <div>
            <label className="block text-sm font-medium text-white mb-2 flex items-center">
              <span className="mr-2">⚡</span>
              Total Delta-V (m/s)
            </label>
            <input
              type="number"
              min="0"
              step="0.01"
              value={config.totalDeltaV}
              onChange={e => setTotalDeltaV(e.target.value)}
              className="form-input"
            />
          </div>

          {/* Payload */}
          <div>
            <label className="block text-sm font-medium text-white mb-2 flex items-center">
              <span className="mr-2">📦</span>
              Payload (kg)
            </label>
            <input
              type="number"
              min="0"
              step="0.01"
              value={config.payload}
              onChange={e => setPayload(e.target.value)}
              className="form-input"
            />
          </div>
        </div>
      </div>

      {/* Stages Configuration */}
      <div className="card">
        <h3 className="text-lg font-semibold text-white mb-6 flex items-center">
          <div className="w-2 h-2 bg-green-400 rounded-full mr-3 animate-pulse"></div>
          Stage Configuration
        </h3>
        
        <div className="overflow-x-auto">
          <table className="w-full">
            <thead>
              <tr className="border-b border-white/20">
                <th className="text-left py-3 px-4 font-semibold text-white">Stage</th>
                <th className="text-center py-3 px-4 font-semibold text-white">Specific Impulse (s)</th>
                <th className="text-center py-3 px-4 font-semibold text-white">Propellant Mass Fraction</th>
              </tr>
            </thead>
            <tbody className="divide-y divide-white/10">
              {config.stages.map((_, index) => {
                return (
                  <ParametersStage
                    key={index}
                    index={index}
                    stages={config.stages}
                    setStages={setStages}
                  />
                );
              })}
            </tbody>
          </table>
        </div>
      </div>

      {/* Actions */}
      <div className="acrylic rounded-xl p-6">
        <div className="flex flex-col sm:flex-row gap-4 items-center justify-between">
          <div className="flex items-center text-sm text-white/80">
            {isOptimizing ? (
              <>
                <div className="animate-spin rounded-full h-4 w-4 border-b-2 border-blue-400 mr-2"></div>
                Optimizing rocket configuration...
              </>
            ) : (
              <>
                <div className="w-2 h-2 bg-green-400 rounded-full mr-2 animate-pulse"></div>
                Configuration ready
              </>
            )}
          </div>
          
          <button
            onClick={() => setRefresh(!refresh)}
            disabled={isOptimizing}
            className="btn btn-primary disabled:opacity-50 disabled:cursor-not-allowed flex items-center"
          >
            <svg className="w-4 h-4 mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15" />
            </svg>
            {isOptimizing ? 'Optimizing...' : 'Recalculate'}
          </button>
        </div>
      </div>
    </div>
  );
};

export default Parameters;
