import axios from 'axios';

// Create axios instance with base configuration
const api = axios.create({
  baseURL: import.meta.env.VITE_API_BASE_URL || 'http://localhost:8000',
  timeout: 30000, // 30 seconds for optimization calculations
  headers: {
    'Content-Type': 'application/json',
  },
});

// Request interceptor for logging
api.interceptors.request.use(
  (config) => {
    console.log(`API Request: ${config.method?.toUpperCase()} ${config.url}`);
    return config;
  },
  (error) => {
    console.error('API Request Error:', error);
    return Promise.reject(error);
  }
);

// Response interceptor for error handling
api.interceptors.response.use(
  (response) => {
    return response;
  },
  (error) => {
    console.error('API Response Error:', error);
    
    // Handle common error scenarios
    if (error.response?.status === 500) {
      throw new Error('Server error occurred during optimization');
    } else if (error.response?.status === 400) {
      throw new Error('Invalid rocket configuration provided');
    } else if (error.code === 'ECONNABORTED') {
      throw new Error('Optimization request timed out');
    } else if (!error.response) {
      throw new Error('Unable to connect to optimization service');
    }
    
    return Promise.reject(error);
  }
);

export default api;