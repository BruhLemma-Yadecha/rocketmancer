import { useState, useEffect, useCallback } from 'react';
import { optimize } from './solver';
import Header from './components/Header';
import Footer from './components/Footer';
import Parameters from './components/Parameters';
import StagesInput from './components/StagesInput';
import { SettingsToggle, SettingsPanel } from './components/Settings';
import Results from './components/Results';
import './styles/index.css';

const DEFAULTS = {
  name: 'Billy Jean',
  totalDeltaV: 4760.08,
  payload: 1000,
  stages: [
    { specificImpulse: 307.36, propellantMassFraction: 0.83 },
    { specificImpulse: 348.81, propellantMassFraction: 0.87 },
  ],
};

function App() {
  const [name, setName] = useState(DEFAULTS.name);
  const [payload, setPayload] = useState(DEFAULTS.payload);
  const [totalDeltaV, setTotalDeltaV] = useState(DEFAULTS.totalDeltaV);
  const [stages, setStages] = useState(DEFAULTS.stages);
  const [result, setResult] = useState(null);
  const [minContribution, setMinContribution] = useState(0);
  const [settingsOpen, setSettingsOpen] = useState(false);

  const runOptimizer = useCallback(() => {
    if (payload <= 0 || totalDeltaV <= 0) return;
    if (
      stages.some(
        s =>
          s.specificImpulse <= 0 || s.propellantMassFraction <= 0 || s.propellantMassFraction >= 1
      )
    )
      return;

    try {
      setResult(optimize(payload, totalDeltaV, stages, minContribution));
    } catch {
      setResult(null);
    }
  }, [payload, totalDeltaV, stages, minContribution]);

  useEffect(() => {
    const t = setTimeout(runOptimizer, 300);
    return () => clearTimeout(t);
  }, [runOptimizer]);

  const updateStage = (index, field, value) => {
    setStages(stages.map((s, i) => (i === index ? { ...s, [field]: value } : s)));
  };

  const setStageCount = count => {
    count = Math.max(1, Math.min(10, parseInt(count) || 1));
    if (count < stages.length) {
      setStages(stages.slice(0, count));
    } else {
      const extra = Array.from({ length: count - stages.length }, () => ({
        specificImpulse: 300,
        propellantMassFraction: 0.9,
      }));
      setStages([...stages, ...extra]);
    }
  };

  return (
    <>
      <Header />
      <main className="container">
        <div className="layout">
          <Results result={result} name={name} />
          <Parameters
            name={name}
            setName={setName}
            payload={payload}
            setPayload={setPayload}
            totalDeltaV={totalDeltaV}
            setTotalDeltaV={setTotalDeltaV}
            stageCount={stages.length}
            setStageCount={setStageCount}
          />
          <StagesInput stages={stages} updateStage={updateStage} />
          {settingsOpen && (
            <SettingsPanel
              minContribution={minContribution}
              setMinContribution={setMinContribution}
              stageCount={stages.length}
            />
          )}
        </div>
      </main>
      <Footer />
      <SettingsToggle open={settingsOpen} setOpen={setSettingsOpen} />
    </>
  );
}

export default App;
