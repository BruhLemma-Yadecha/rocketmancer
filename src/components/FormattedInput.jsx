import { useState } from 'react';
import { fmt } from '../palette';

export default function FormattedInput({ value, onChange, decimals = 2, ...props }) {
  const [focused, setFocused] = useState(false);

  return (
    <input
      {...props}
      type={focused ? 'number' : 'text'}
      value={focused ? value : fmt(value, decimals)}
      onFocus={e => {
        setFocused(true);
        requestAnimationFrame(() => e.target.select());
      }}
      onBlur={() => setFocused(false)}
      onChange={e => onChange(parseFloat(e.target.value) || 0)}
    />
  );
}
