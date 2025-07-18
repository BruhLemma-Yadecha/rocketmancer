/**
 * Utility functions for formatting values in the UI
 */

/**
 * Format mass values with appropriate units
 * @param {number} mass - Mass in kg
 * @returns {string} Formatted mass string
 */
export const formatMass = (mass) => {
  if (mass >= 1000) {
    return `${(mass / 1000).toFixed(2)} t`;
  }
  return `${mass.toFixed(2)} kg`;
};

/**
 * Format velocity values
 * @param {number} velocity - Velocity in m/s
 * @returns {string} Formatted velocity string
 */
export const formatVelocity = (velocity) => {
  return `${velocity.toFixed(0)} m/s`;
};

/**
 * Format specific impulse values
 * @param {number} isp - Specific impulse in seconds
 * @returns {string} Formatted Isp string
 */
export const formatSpecificImpulse = (isp) => {
  return `${isp.toFixed(0)} s`;
};

/**
 * Format mass ratio values
 * @param {number} ratio - Mass ratio
 * @returns {string} Formatted ratio string
 */
export const formatMassRatio = (ratio) => {
  return ratio.toFixed(3);
};

/**
 * Format percentage values
 * @param {number} fraction - Fraction (0-1)
 * @returns {string} Formatted percentage string
 */
export const formatPercentage = (fraction) => {
  return `${(fraction * 100).toFixed(1)}%`;
};

/**
 * Format numbers with thousands separators
 * @param {number} num - Number to format
 * @returns {string} Formatted number string
 */
export const formatNumber = (num) => {
  return new Intl.NumberFormat().format(num);
};