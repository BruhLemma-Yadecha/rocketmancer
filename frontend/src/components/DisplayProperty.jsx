const DECIMAL_PLACES = 4;
const MASS_UNIT = 't';
const VELOCITY_UNIT = 'm/s';
const SPECIFIC_IMPULSE_UNIT = 's';
const typeToUnit = {
  velocity: VELOCITY_UNIT,
  mass: MASS_UNIT,
  time: SPECIFIC_IMPULSE_UNIT,
};

function fixDecimals(number) {
  return number.toFixed(DECIMAL_PLACES);
}

const getTypeIcon = (type) => {
  switch (type) {
    case 'velocity': return '🚀';
    case 'mass': return '⚖️';
    case 'time': return '⏱️';
    case 'ratio': return '📊';
    case 'percentage': return '📈';
    default: return '📋';
  }
};

const getTypeColor = (type) => {
  switch (type) {
    case 'velocity': return 'text-blue-300';
    case 'mass': return 'text-green-300';
    case 'time': return 'text-purple-300';
    case 'ratio': return 'text-orange-300';
    case 'percentage': return 'text-red-300';
    default: return 'text-white/80';
  }
};

const DisplayProperty = ({ name, type, stages }) => {
  const convertName = () => {
    const firstLetter = name[0].toLowerCase();
    const withoutFirstLetter = name.slice(1);
    const withoutSpaces = withoutFirstLetter.replace(/ /g, '');
    return firstLetter + withoutSpaces;
  };
  
  const convertedName = convertName();
  const icon = getTypeIcon(type);
  const colorClass = getTypeColor(type);
  
  return (
    <tr className="hover:bg-white/5 transition-colors duration-200">
      <td className="py-3 px-4">
        <div className="flex items-center">
          <span className="mr-2 text-sm">{icon}</span>
          <div>
            <span className={`font-medium ${colorClass}`}>{name}</span>
            {(type !== 'ratio' && type !== 'percentage') && (
              <span className="text-sm text-white/60 ml-1">({typeToUnit[type]})</span>
            )}
          </div>
        </div>
      </td>
      {stages.map((stage, index) => (
        <td key={index} className="py-3 px-4 text-center">
          <span className="inline-flex items-center px-2.5 py-0.5 rounded-full text-sm font-medium bg-white/10 text-white backdrop-blur-md border border-white/20">
            {fixDecimals(stage[convertedName])}
          </span>
        </td>
      ))}
    </tr>
  );
};

export default DisplayProperty;
