import { useState, useEffect } from 'react';
import axios from 'axios';
import ParametersStage from './ParametersStage';
import '../../styles/Parameters.css';

const Parameters = ({ setRocket, rocketName, setRocketName }) => {
  const heading = ['Stage', 'Specific Impulse (s)', 'Propellant Mass Fraction'];

  const [config, setConfig] = useState(undefined);
  const [refresh, setRefresh] = useState(false);

  useEffect(() => {
    axios.get('/billyJean.json').then(response => {
      setConfig(response.data);
      setRocketName(response.data.name);
    });
  }, []);

  useEffect(() => {
    if (!config) return;
    axios.post('/api/optimize/', config).then(response => {
      setRocket(response.data.result);
    });
  }, [refresh]);

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

  if (!config) return <div>Loading...</div>;

  return (
    <div className={'parameters-display'}>
      <h2>Parameters</h2>
      <div>
        <label>Name: </label>
        <input
          className={'parameters-input'}
          type="text"
          value={rocketName}
          onChange={e => setName(e.target.value)}
        />{' '}
        <br />
      </div>
      <div>
        <label>Total Stages: </label>
        <input
          className={'parameters-input'}
          type="number"
          value={config.totalStages}
          onChange={e => setTotalStages(e.target.value)}
        />{' '}
        <br />
      </div>
      <div>
        <label>Total Delta-V: </label>
        <input
          className={'parameters-input'}
          type="number"
          value={config.totalDeltaV}
          onChange={e => setTotalDeltaV(e.target.value)}
        />{' '}
        <br />
      </div>
      <div>
        <label>Payload: </label>
        <input
          className={'parameters-input'}
          type="number"
          value={config.payload}
          onChange={e => setPayload(e.target.value)}
        />{' '}
        <br />
      </div>
      <h3>Stages</h3>
      <table>
        <thead>
          <tr>
            {heading.map(head => (
              <th key={head}>{head}</th>
            ))}
          </tr>
        </thead>
        <tbody>
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
      <button className={'parameters-refresh'} onClick={e => setRefresh(!refresh)}>
        Refresh
      </button>
    </div>
  );
};

export default Parameters;
