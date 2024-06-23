function dVdt = elementAcceleration(x, y, xi, eta, xi_w, eta_w, element, N)
%ELEMENTACCELERATION Computes the time derivative of velocity integrated across a 2D element.
%   dVdt = ELEMENTACCELERATION(X, Y, XI, ETA, XI_W, ETA_W, ELEMENT, N)
%   returns the matrix representing the time derivative of velocity integrated
%   across the given element using the provided quadrature points and weights.
%
%   This function evaluates the integral:
%   \frac{dV}{dt} ≡ ∫∫_D N^T N dx dy
%
%   Syntax:
%   dVdt = ELEMENTACCELERATION(X, Y, XI, ETA, XI_W, ETA_W, ELEMENT, N)
%
%   Inputs:
%   X - X-coordinates of the element nodes (numeric vector)
%   Y - Y-coordinates of the element nodes (numeric vector)
%   XI - Quadrature points in the xi direction (numeric vector)
%   ETA - Quadrature points in the eta direction (numeric vector)
%   XI_W - Quadrature weights in the xi direction (numeric vector)
%   ETA_W - Quadrature weights in the eta direction (numeric vector)
%   ELEMENT - Handle to the element function that returns shape functions and their derivatives
%   N - Number of shape functions (integer)
%
%   Outputs:
%   dVdt - Matrix representing the time derivative of velocity integrated across the element (N x N numeric matrix)
%
%   Example:
%   x = [0, 1, 1, 0];
%   y = [0, 0, 1, 1];
%   xi = [-1, 1];
%   eta = [-1, 1];
%   xi_w = [1, 1];
%   eta_w = [1, 1];
%   element = @linearQuadrilateral;
%   N = 4;
%   dVdt = elementAcceleration(x, y, xi, eta, xi_w, eta_w, element, N);
%
%   See also: ELEMENTDIFFUSIONMATRIX, LINEARQUADRILATERAL, LINEARTRIANGLE, QUADRATICQUADRILATERAL

    % Initialize the matrix for the time derivative of velocity
    dVdt = zeros(N, N);

    % Loop over quadrature points
    for i = 1:numel(xi)
        for j = 1:numel(eta)
            [N_vals, ~, J] = element(x, y, xi(i), eta(j));
            dVdt = dVdt + det(J) * xi_w(i) * eta_w(j) * (transpose(N_vals) .* N_vals);
        end
    end
end