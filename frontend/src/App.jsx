import { useState } from 'react';
import Rocket from './components/Rocket';
import Title from './components/Title';
import Parameters from './components/Parameters/Parameters';
import Footer from './Footer';
import './styles/index.css';

function App() {
  const [rocket, setRocket] = useState(undefined);
  const [rocketName, setRocketName] = useState(undefined);

  return (
    <div className="min-h-screen flex flex-col">
      <Title />

      <main className="flex-1 container mx-auto px-4 py-3">
        <div className="grid grid-cols-1 lg:grid-cols-5 gap-3 max-w-7xl mx-auto">
          {/* Left Column - Compact configuration cards */}
          <div className="lg:col-span-2 space-y-3 order-2 lg:order-1">
            <Parameters
              rocketName={rocketName}
              setRocketName={setRocketName}
              setRocket={setRocket}
            />
          </div>

          {/* Right Column - Rocket Display (takes 3/5 of the space) */}
          <div className="lg:col-span-3 order-1 lg:order-2">
            <Rocket rocket={rocket} rocketName={rocketName} />
          </div>
        </div>
      </main>

      <Footer />
    </div>
  );
}

export default App;
