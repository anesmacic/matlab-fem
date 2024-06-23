function [N, dNdX, J] = quadraticTriangle(x, y, xi, eta)
%QUADRATICTRIANGLE Computes the shape functions and their derivatives for a 6-node quadratic triangular element.
%   [N, dNdX, J] = QUADRATICTRIANGLE(X, Y, XI, ETA) returns the shape functions N, 
%   their derivatives with respect to global coordinates dNdX, and the Jacobian matrix J 
%   for the given element at the specified local coordinates (XI, ETA).
%
%   Syntax:
%   [N, dNdX, J] = QUADRATICTRIANGLE(X, Y, XI, ETA)
%
%   Inputs:
%   X - X-coordinates of the element nodes (numeric vector of length 6)
%   Y - Y-coordinates of the element nodes (numeric vector of length 6)
%   XI - Local coordinate in the xi direction (scalar)
%   ETA - Local coordinate in the eta direction (scalar)
%
%   Outputs:
%   N - Shape functions evaluated at (XI, ETA) (1x6 numeric vector)
%   dNdX - Derivatives of the shape functions with respect to global coordinates (2x6 numeric matrix)
%   J - Jacobian matrix evaluated at (XI, ETA) (2x2 numeric matrix)
%
%   Example:
%   x = [0, 1, 0, 0.5, 0.5, 0];
%   y = [0, 0, 1, 0, 0.5, 0.5];
%   xi = 0.3;
%   eta = 0.3;
%   [N, dNdX, J] = quadraticTriangle(x, y, xi, eta);
%
%   See also: ELEMENTDIFFUSIONMATRIX, LINEARQUADRILATERAL, LINEARTRIANGLE, QUADRATICQUADRILATERAL

    % Define the shape functions for a 6-node quadratic triangle
    N = [(1 - xi - eta) * (1 - 2 * xi - 2 * eta), ...
         xi * (2 * xi - 1), ...
         eta * (2 * eta - 1), ...
         4 * xi * (1 - xi - eta), ...
         4 * xi * eta, ...
         4 * eta * (1 - xi - eta)];
    
    % Derivatives of the shape functions with respect to local coordinates
    dNdxi = [-3 + 4 * xi + 4 * eta, 4 * xi - 1, 0, 4 * (1 - 2 * xi - eta), 4 * eta, -4 * eta;
             -3 + 4 * xi + 4 * eta, 0, 4 * eta - 1, -4 * xi, 4 * xi, 4 * (1 - xi - 2 * eta)];
    
    % Initialize the Jacobian matrix
    J = zeros(2, 2);
    
    % Calculate the Jacobian matrix
    for i = 1:6
        J(1, :) = J(1, :) + dNdxi(1, i) * [x(i), y(i)];
        J(2, :) = J(2, :) + dNdxi(2, i) * [x(i), y(i)];
    end
    
    % Compute derivatives of N with respect to global coordinates
    dNdX = J \ dNdxi;
end