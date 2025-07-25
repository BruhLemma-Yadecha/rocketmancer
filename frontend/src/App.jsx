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
    <div className="min-h-screen bg-gradient-to-br from-indigo-600 via-purple-600 to-purple-700">
      {/* Background overlay */}
      <div className="absolute inset-0 bg-black/5"></div>
      
      <div className="relative z-10 flex flex-col min-h-screen">
        <Title />
        
        <main className="flex-1 container mx-auto px-4 py-8">
          <div className="grid grid-cols-1 xl:grid-cols-2 gap-8 max-w-7xl mx-auto">
            {/* Parameters Panel */}
            <div className="order-2 xl:order-1">
              <Parameters 
                rocketName={rocketName} 
                setRocketName={setRocketName} 
                setRocket={setRocket} 
              />
            </div>
            
            {/* Rocket Display Panel */}
            <div className="order-1 xl:order-2">
              <Rocket rocket={rocket} rocketName={rocketName} />
            </div>
          </div>
        </main>
        
        <Footer />
      </div>
    </div>
  );
}

export default App;
