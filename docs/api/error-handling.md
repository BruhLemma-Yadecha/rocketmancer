# Error Handling

This document describes how errors are handled in the Rocketmancer Web API.

## HTTP Status Codes

The API uses standard HTTP status codes to indicate success or failure:

- `200 OK`: Request successful
- `400 Bad Request`: Invalid request data or malformed JSON
- `405 Method Not Allowed`: HTTP method not supported for endpoint
- `500 Internal Server Error`: Server-side error during processing

## Error Response Format

When an error occurs, the API returns a JSON response with error details:

```json
{
  "error": "Error message describing what went wrong",
  "details": "Additional details about the error (when available)"
}
```

## Common Errors

### Invalid JSON (400)

```json
{
  "error": "Invalid JSON in request body"
}
```

### Missing Parameters (400)

```json
{
  "error": "Missing required parameter: payload"
}
```

### Method Not Allowed (405)

```json
{
  "error": "Method 'GET' not allowed. Use POST instead."
}
```

### Calculation Error (500)

```json
{
  "error": "Optimization calculation failed",
  "details": "Invalid stage configuration or impossible mission requirements"
}
```

## Best Practices

1. Always check the HTTP status code before processing the response
2. Handle network errors and timeouts gracefully
3. Validate input data on the client side before sending requests
4. Log errors for debugging purposes
5. Provide user-friendly error messages in your application

## Debugging Tips

- Use browser developer tools to inspect network requests
- Check the request payload format matches the expected structure
- Verify all required parameters are included
- Ensure numeric values are within reasonable ranges