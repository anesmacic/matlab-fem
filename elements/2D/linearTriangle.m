function [N, dNdX, J] = linearTriangle(x, y, xi, eta)
%LINEARTRIANGLE Computes the shape functions and their derivatives for a 3-node linear triangular element.
%   [N, dNdX, J] = LINEARTRIANGLE(X, Y, XI, ETA) returns the shape functions N, 
%   their derivatives with respect to global coordinates dNdX, and the Jacobian matrix J 
%   for the given element at the specified local coordinates (XI, ETA).
%
%   Syntax:
%   [N, dNdX, J] = LINEARTRIANGLE(X, Y, XI, ETA)
%
%   Inputs:
%   X - X-coordinates of the element nodes (numeric vector of length 3)
%   Y - Y-coordinates of the element nodes (numeric vector of length 3)
%   XI - Local coordinate in the xi direction (scalar)
%   ETA - Local coordinate in the eta direction (scalar)
%
%   Outputs:
%   N - Shape functions evaluated at (XI, ETA) (1x3 numeric vector)
%   dNdX - Derivatives of the shape functions with respect to global coordinates (2x3 numeric matrix)
%   J - Jacobian matrix evaluated at (XI, ETA) (2x2 numeric matrix)
%
%   Example:
%   x = [0, 1, 0];
%   y = [0, 0, 1];
%   xi = 0.5;
%   eta = 0.2;
%   [N, dNdX, J] = linearTriangle(x, y, xi, eta);
%
%   See also: ELEMENTDIFFUSIONMATRIX

    % Define the shape functions for a 3-node linear triangle
    N = [1 - xi - eta, xi, eta];
    
    % Derivatives of the shape functions with respect to local coordinates
    dNdxi = [-1, 1, 0;
             -1, 0, 1];
    
    % Initialize the Jacobian matrix
    J = zeros(2, 2);
    
    % Calculate the Jacobian matrix
    for i = 1:3
        J(1, :) = J(1, :) + dNdxi(1, i)*[x(i), y(i)];
        J(2, :) = J(2, :) + dNdxi(2, i)*[x(i), y(i)];
    end
    
    % Compute derivatives of N with respect to global coordinates
    dNdX = J \ dNdxi;
end