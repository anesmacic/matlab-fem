function [N, dNdX, J] = linearQuadrilateral(x, y, xi, eta)
%LINEARQUADRILATERAL Computes the shape functions and their derivatives for a 4-node linear quadrilateral element.
%   [N, dNdX, J] = LINEARQUADRILATERAL(X, Y, XI, ETA) returns the shape functions N, 
%   their derivatives with respect to global coordinates dNdX, and the Jacobian matrix J 
%   for the given element at the specified local coordinates (XI, ETA).
%
%   Syntax:
%   [N, dNdX, J] = LINEARQUADRILATERAL(X, Y, XI, ETA)
%
%   Inputs:
%   X - X-coordinates of the element nodes (numeric vector of length 4)
%   Y - Y-coordinates of the element nodes (numeric vector of length 4)
%   XI - Local coordinate in the xi direction (scalar)
%   ETA - Local coordinate in the eta direction (scalar)
%
%   Outputs:
%   N - Shape functions evaluated at (XI, ETA) (1x4 numeric vector)
%   dNdX - Derivatives of the shape functions with respect to global coordinates (2x4 numeric matrix)
%   J - Jacobian matrix evaluated at (XI, ETA) (2x2 numeric matrix)
%
%   Example:
%   x = [0, 1, 1, 0];
%   y = [0, 0, 1, 1];
%   xi = 0.5;
%   eta = -0.5;
%   [N, dNdX, J] = linearQuadrilateral(x, y, xi, eta);
%
%   See also: ELEMENTDIFFUSIONMATRIX

    % Define the shape functions for a 4-node linear quadrilateral
    N = [(1-xi)*(1-eta)/4,  ...
         (1+xi)*(1-eta)/4, ...
         (1+xi)*(1+eta)/4, ...
         (1-xi)*(1+eta)/4];
    
    % Derivatives of the shape functions with respect to local coordinates
    dNdxi = [-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4;
             -(1-xi)/4, -(1+xi)/4, (1+xi)/4, (1-xi)/4];
    
    % Initialize the Jacobian matrix
    J = zeros(2, 2);
    
    % Calculate the Jacobian matrix
    for i = 1:4
        J(1, :) = J(1, :) + dNdxi(1, i)*[x(i), y(i)];
        J(2, :) = J(2, :) + dNdxi(2, i)*[x(i), y(i)];
    end
    
    % Compute derivatives of N with respect to global coordinates
    dNdX = J \ dNdxi;
end