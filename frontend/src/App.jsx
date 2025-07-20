import { useState } from 'react';
import Rocket from './components/Rocket';
import Title from './components/Title';
import Parameters from './components/Parameters/Parameters';
import './styles/App.css';

function App() {
  const [rocket, setRocket] = useState(undefined);
  const [rocketName, setRocketName] = useState(undefined);
  return (
    <div className={'site-body'}>
      <Title />
      <Rocket rocket={rocket} rocketName={rocketName} />
      <Parameters rocketName={rocketName} setRocketName={setRocketName} setRocket={setRocket} />
    </div>
  );
}

export default App;
