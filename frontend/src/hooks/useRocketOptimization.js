import { useState, useCallback } from 'react';
import { RocketService } from '../services/rocketService';

/**
 * Custom hook for managing rocket optimization state and operations
 */
export const useRocketOptimization = () => {
  const [rocket, setRocket] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const optimizeRocket = useCallback(async rocketConfig => {
    setLoading(true);
    setError(null);

    try {
      // Validate configuration first
      const validation = RocketService.validateConfiguration(rocketConfig);
      if (!validation.isValid) {
        throw new Error(validation.errors.join(', '));
      }

      // Perform optimization
      const optimizedRocket = await RocketService.optimize(rocketConfig);
      setRocket(optimizedRocket);
      return optimizedRocket;
    } catch (err) {
      const errorMessage = err.message || 'Optimization failed';
      setError(errorMessage);
      throw err;
    } finally {
      setLoading(false);
    }
  }, []);

  const clearRocket = useCallback(() => {
    setRocket(null);
    setError(null);
  }, []);

  const clearError = useCallback(() => {
    setError(null);
  }, []);

  return {
    rocket,
    loading,
    error,
    optimizeRocket,
    clearRocket,
    clearError,
  };
};
