import { useState } from 'react';

const ParametersStage = ({ index, stages, setStages }) => {
  const [specificImpulse, setSpecificImpulse] = useState(stages[index].specificImpulse);
  const [propellantMassFraction, setPropellantMassFraction] = useState(
    stages[index].propellantMassFraction
  );

  const editSpecificImpulse = specificImpulse => {
    specificImpulse = parseFloat(specificImpulse);
    if (specificImpulse < 0) return;
    setSpecificImpulse(specificImpulse);

    const newStages = [...stages];
    newStages[index] = { ...newStages[index], specificImpulse };
    setStages(newStages);
  };

  const editPropellantMassFraction = propellantMassFraction => {
    propellantMassFraction = parseFloat(propellantMassFraction);
    if (propellantMassFraction < 0) return;
    setPropellantMassFraction(propellantMassFraction);

    const newStages = [...stages];
    newStages[index] = { ...newStages[index], propellantMassFraction };
    setStages(newStages);
  };

  const getStageColor = stageIndex => {
    const colors = [
      'bg-blue-50 text-blue-700',
      'bg-green-50 text-green-700',
      'bg-purple-50 text-purple-700',
      'bg-orange-50 text-orange-700',
      'bg-red-50 text-red-700',
      'bg-indigo-50 text-indigo-700',
      'bg-pink-50 text-pink-700',
      'bg-yellow-50 text-yellow-700',
    ];
    return colors[stageIndex % colors.length];
  };

  return (
    <tr className="hover:bg-white/5 transition-colors duration-200">
      <td className="py-1 px-2">
        <div
          className={`inline-flex items-center justify-center w-4 h-4 rounded-full text-xs font-semibold ${getStageColor(index)}`}
        >
          {index + 1}
        </div>
      </td>
      <td className="py-1 px-2">
        <div className="relative">
          <input
            type="number"
            min="0"
            step="0.01"
            value={specificImpulse}
            onChange={e => editSpecificImpulse(e.target.value)}
            className="form-input text-center text-xs py-1 px-1"
            placeholder="300"
          />
          <div className="absolute inset-y-0 right-0 flex items-center pr-1 pointer-events-none">
            <span className="text-white/60 text-xs">s</span>
          </div>
        </div>
      </td>
      <td className="py-1 px-2">
        <div className="relative">
          <input
            type="number"
            min="0"
            max="1"
            step="0.01"
            value={propellantMassFraction}
            onChange={e => editPropellantMassFraction(e.target.value)}
            className="form-input text-center text-xs py-1 px-1"
            placeholder="0.85"
          />
          <div className="absolute inset-y-0 right-0 flex items-center pr-1 pointer-events-none">
            <span className="text-white/60 text-xs">%</span>
          </div>
        </div>
      </td>
    </tr>
  );
};

export default ParametersStage;
