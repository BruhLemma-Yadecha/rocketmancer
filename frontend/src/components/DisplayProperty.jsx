import '../styles/DisplayProperty.css';

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
  // take number to four decimal places
  return number.toFixed(DECIMAL_PLACES);
}

const DisplayProperty = ({ name, type, stages }) => {
  const convertName = () => {
    const firstLetter = name[0].toLowerCase();
    const withoutFirstLetter = name.slice(1);
    const withoutSpaces = withoutFirstLetter.replace(/ /g, '');
    return firstLetter + withoutSpaces;
  };
  const convertedName = convertName();
  return (
    <tr>
      <td>
        <b>
          {name} {type == 'ratio' || type == 'percentage' ? '' : `(${typeToUnit[type]})`}
        </b>{' '}
      </td>
      {stages.map((stage, index) => (
        <td key={index} className={'table-element'}>
          {fixDecimals(stage[convertedName])}
        </td>
      ))}
    </tr>
  );
};

export default DisplayProperty;
