import { useState, useEffect } from 'react';
import RocketService from '../../services/rocketService';
import ParametersStage from './ParametersStage';

const INPUT_DEBOUNCE_MS = 500;

const Parameters = ({ setRocket, rocketName, setRocketName }) => {
  const [config, setConfig] = useState(undefined);
  const [refresh, setRefresh] = useState(false);
  const [isOptimizing, setIsOptimizing] = useState(false);

  const loadDefaultConfiguration = () => {
    const defaultConfig = {
      name: 'Billy Jean',
      totalStages: 2,
      totalDeltaV: 4760.08,
      payload: 1000.0,
      stages: [
        {
          specificImpulse: 307.36,
          propellantMassFraction: 0.83,
        },
        {
          specificImpulse: 348.81,
          propellantMassFraction: 0.87,
        },
      ],
    };
    setConfig(defaultConfig);
    setRocketName(defaultConfig.name);
  };

  useEffect(() => {
    if (!config) {
      loadDefaultConfiguration();
      return;
    }

    const debounceTimeout = setTimeout(() => {
      const optimizeRocket = async () => {
        if (!config) return;

        setIsOptimizing(true);
        try {
          const validation = RocketService.validateConfiguration(config);
          if (!validation.isValid) {
            console.warn('Configuration validation failed:', validation.errors);
          }

          const optimizedRocket = await RocketService.optimize(config);
          setRocket(optimizedRocket);
        } catch (error) {
          console.error('Optimization failed:', error);
        } finally {
          setIsOptimizing(false);
        }
      };
      optimizeRocket();
    }, INPUT_DEBOUNCE_MS);

    return () => clearTimeout(debounceTimeout);
  }, [config, refresh]);

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
    <>
      {/* Card 1: Basic Parameters */}
      <div className="card-glass animate-fade-in">
        <h2 className="text-base font-bold text-white mb-2 flex items-center">
          <div className="w-1.5 h-1.5 bg-blue-400 rounded-full mr-2 animate-pulse"></div>
          Parameters
        </h2>

        <div className="space-y-2">
          {/* Rocket Name - Full Width */}
          <div>
            <label className="block text-xs font-medium text-white/70 mb-0.5 flex items-center">
              <span className="mr-1 text-xs">🚀</span>
              Name
            </label>
            <input
              type="text"
              value={rocketName}
              onChange={e => setName(e.target.value)}
              className="form-input text-xs py-1 px-2 w-full"
              placeholder="Enter rocket name..."
            />
          </div>

          {/* Mini Cards Row - Stages, Delta-V, Payload */}
          <div className="grid grid-cols-3 gap-2">
            {/* Total Stages */}
            <div className="bg-white/5 rounded-lg p-2 backdrop-blur-sm border border-white/10">
              <label className="block text-xs font-medium text-white/70 mb-1 flex items-center">
                <span className="mr-1 text-xs">🔢</span>
                Stages
              </label>
              <input
                type="number"
                min="1"
                max="10"
                value={config.totalStages}
                onChange={e => setTotalStages(e.target.value)}
                className="form-input text-xs py-1 px-1 w-full text-center"
              />
            </div>

            {/* Total Delta-V */}
            <div className="bg-white/5 rounded-lg p-2 backdrop-blur-sm border border-white/10">
              <label className="block text-xs font-medium text-white/70 mb-1 flex items-center">
                <span className="mr-1 text-xs">⚡</span>
                Δv
              </label>
              <input
                type="number"
                min="0"
                step="0.01"
                value={config.totalDeltaV}
                onChange={e => setTotalDeltaV(e.target.value)}
                className="form-input text-xs py-1 px-1 w-full text-center"
                placeholder="4760"
              />
            </div>

            {/* Payload */}
            <div className="bg-white/5 rounded-lg p-2 backdrop-blur-sm border border-white/10">
              <label className="block text-xs font-medium text-white/70 mb-1 flex items-center">
                <span className="mr-1 text-xs">📦</span>
                Mass
              </label>
              <input
                type="number"
                min="0"
                step="0.01"
                value={config.payload}
                onChange={e => setPayload(e.target.value)}
                className="form-input text-xs py-1 px-1 w-full text-center"
                placeholder="1000"
              />
            </div>
          </div>
        </div>
      </div>

      {/* Card 2: Stage Configuration */}
      <div className="card-glass animate-fade-in pb-4">
        <div className="flex items-center justify-between mb-2">
          <h3 className="text-base font-bold text-white flex items-center">
            <div className="w-1.5 h-1.5 bg-green-400 rounded-full mr-2 animate-pulse"></div>
            Stages
            {isOptimizing && (
              <div className="ml-2 flex items-center text-xs text-white/60">
                <div className="animate-spin rounded-full h-1.5 w-1.5 border border-white/40 border-t-white mr-1"></div>
                optimizing
              </div>
            )}
          </h3>
          <button
            onClick={() => setRefresh(!refresh)}
            disabled={isOptimizing}
            className="glass p-1 rounded hover:glass-strong transition-all duration-200 disabled:opacity-50 disabled:cursor-not-allowed group"
            title="Force recalculate"
          >
            <svg
              className={`w-3 h-3 text-white/80 group-hover:text-white transition-colors duration-200 ${isOptimizing ? 'animate-spin' : ''}`}
              fill="none"
              stroke="currentColor"
              viewBox="0 0 24 24"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15"
              />
            </svg>
          </button>
        </div>

        <div className="overflow-x-auto">
          <table className="w-full">
            <thead>
              <tr className="border-b border-white/20">
                <th className="text-left py-1 px-2 font-medium text-white text-xs">#</th>
                <th className="text-center py-1 px-2 font-medium text-white text-xs">
                  Specific Impulse (s)
                </th>
                <th className="text-center py-1 px-2 font-medium text-white text-xs">
                  Mass Fraction
                </th>
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
    </>
  );
};

export default Parameters;
