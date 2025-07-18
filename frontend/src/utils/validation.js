/**
 * Validation utilities for form inputs
 */

/**
 * Validate numeric input
 * @param {string|number} value - Value to validate
 * @param {Object} options - Validation options
 * @returns {Object} Validation result
 */
export const validateNumber = (value, options = {}) => {
  const { min, max, required = true } = options;
  const num = parseFloat(value);

  if (required && (value === '' || value === null || value === undefined)) {
    return { isValid: false, error: 'This field is required' };
  }

  if (isNaN(num)) {
    return { isValid: false, error: 'Must be a valid number' };
  }

  if (min !== undefined && num < min) {
    return { isValid: false, error: `Must be at least ${min}` };
  }

  if (max !== undefined && num > max) {
    return { isValid: false, error: `Must be at most ${max}` };
  }

  return { isValid: true, error: null };
};

/**
 * Validate mass fraction (0 < value < 1)
 * @param {string|number} value - Mass fraction value
 * @returns {Object} Validation result
 */
export const validateMassFraction = (value) => {
  const result = validateNumber(value, { min: 0.01, max: 0.99 });
  
  if (!result.isValid) {
    return result;
  }

  const num = parseFloat(value);
  if (num <= 0 || num >= 1) {
    return { isValid: false, error: 'Mass fraction must be between 0 and 1' };
  }

  return { isValid: true, error: null };
};

/**
 * Validate specific impulse
 * @param {string|number} value - Specific impulse value
 * @returns {Object} Validation result
 */
export const validateSpecificImpulse = (value) => {
  return validateNumber(value, { min: 50, max: 500 });
};

/**
 * Validate payload mass
 * @param {string|number} value - Payload mass value
 * @returns {Object} Validation result
 */
export const validatePayload = (value) => {
  return validateNumber(value, { min: 0.1, max: 1000000 });
};

/**
 * Validate delta-V
 * @param {string|number} value - Delta-V value
 * @returns {Object} Validation result
 */
export const validateDeltaV = (value) => {
  return validateNumber(value, { min: 100, max: 20000 });
};