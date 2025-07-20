import '../styles/Rocket.css';
import DisplayProperty from './DisplayProperty';

const MASS_UNIT = 't';
const VELOCITY_UNIT = 'm/s';
const SPECIFIC_IMPULSE_UNIT = 's';
const DECIMAL_PLACES = 4;

const Rocket = ({ rocket, rocketName }) => {
  if (!rocket) return <div>Loading...</div>;

  const properties = [
    { name: 'Delta V', type: 'velocity' },
    { name: 'Mass Ratio', type: 'ratio' },
    { name: 'Payload Mass', type: 'mass' },
    { name: 'Wet Mass', type: 'mass' },
    { name: 'Dry Mass', type: 'mass' },
    { name: 'Structural Mass', type: 'mass' },
    { name: 'Propellant Mass', type: 'mass' },
    { name: 'Exhaust Velocity', type: 'velocity' },
    { name: 'Specific Impulse', type: 'time' },
    { name: 'Propellant Mass Fraction', type: 'percentage' },
  ];
  return (
    <div className={'rocket-display'}>
      <h2 className={'rocket-display-name'}>{rocketName}</h2>
      <div>
        <b>Stages: </b> {rocket.totalStages} <br />
      </div>
      <div>
        <b>Total Delta V: </b> {rocket.totalDeltaV.toFixed(DECIMAL_PLACES)} {VELOCITY_UNIT} <br />
      </div>
      <div>
        <b>Payload ({MASS_UNIT}): </b> {rocket.payload.toFixed(DECIMAL_PLACES)} kg <br />
      </div>
      <div>
        <b>Total Mass ({MASS_UNIT}): </b> {rocket.totalMass.toFixed(DECIMAL_PLACES)} kg <br />
      </div>
      &nbsp;
      <table>
        <tbody>
          <tr>
            <td></td>
            {rocket.stages.map((stage, index) => (
              <td key={index} className={'rocket-display-header'}>
                <b>Stage {index + 1}</b>
              </td>
            ))}
          </tr>
          {properties.map((property, index) => (
            <DisplayProperty
              key={index}
              name={property.name}
              type={property.type}
              stages={rocket.stages}
            />
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default Rocket;
