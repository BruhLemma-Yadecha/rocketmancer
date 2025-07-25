import api from './api';

/**
 * Service for rocket optimization API calls
 */
export class RocketService {
  /**
   * Optimize rocket configuration
   * @param {Object} rocketConfig - Rocket configuration object
   * @param {number} rocketConfig.payload - Payload mass in kg
   * @param {number} rocketConfig.totalDeltaV - Total delta-V requirement in m/s
   * @param {number} rocketConfig.totalStages - Number of stages
   * @param {Array} rocketConfig.stages - Array of stage configurations
   * @returns {Promise<Object>} Optimized rocket configuration
   */
  static async optimize(rocketConfig) {
    try {
      const response = await api.post('/api/v1/optimize/', rocketConfig);
      return response.data.result;
    } catch (error) {
      console.error('Rocket optimization failed:', error);
      throw error;
    }
  }

  /**
   * Validate rocket configuration before optimization
   * @param {Object} rocketConfig - Rocket configuration to validate
   * @returns {Object} Validation result with isValid boolean and errors array
   */
  static validateConfiguration(rocketConfig) {
    const errors = [];

    // Validate basic parameters
    if (!rocketConfig.payload || rocketConfig.payload <= 0) {
      errors.push('Payload mass must be greater than 0');
    }

    if (!rocketConfig.totalDeltaV || rocketConfig.totalDeltaV <= 0) {
      errors.push('Total delta-V must be greater than 0');
    }

    if (!rocketConfig.totalStages || rocketConfig.totalStages < 1) {
      errors.push('Must have at least 1 stage');
    }

    if (!rocketConfig.stages || rocketConfig.stages.length !== rocketConfig.totalStages) {
      errors.push('Number of stages must match totalStages parameter');
    }

    // Validate each stage
    if (rocketConfig.stages) {
      rocketConfig.stages.forEach((stage, index) => {
        if (!stage.specificImpulse || stage.specificImpulse <= 0) {
          errors.push(`Stage ${index + 1}: Specific impulse must be greater than 0`);
        }

        if (
          !stage.propellantMassFraction ||
          stage.propellantMassFraction <= 0 ||
          stage.propellantMassFraction >= 1
        ) {
          errors.push(`Stage ${index + 1}: Propellant mass fraction must be between 0 and 1`);
        }
      });
    }

    return {
      isValid: errors.length === 0,
      errors,
    };
  }
}

export default RocketService;
