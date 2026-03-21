import { describe, it, expect } from 'vitest';
import { render, screen } from '@testing-library/react';
import App from '../App';

describe('App', () => {
  it('renders without crashing', () => {
    render(<App />);
    expect(screen.getAllByText(/rocketmancer/i).length).toBeGreaterThan(0);
  });

  it('renders the main heading', () => {
    render(<App />);
    expect(screen.getByRole('heading', { level: 1, name: /rocketmancer/i })).toBeInTheDocument();
  });

  it('renders parameter inputs', () => {
    render(<App />);
    expect(screen.getByText('Parameters')).toBeInTheDocument();
    expect(screen.getAllByText('Stages').length).toBeGreaterThan(0);
  });
});
