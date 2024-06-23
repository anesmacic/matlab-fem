function [N, dNdX, J] = quadraticQuadrilateral(x, y, xi, eta)
%QUADRATICQUADRILATERAL Computes the shape functions and their derivatives for a 9-node Lagrange quadratic quadrilateral element.
%   [N, dNdX, J] = QUADRATICQUADRILATERAL(X, Y, XI, ETA) returns the shape functions N, 
%   their derivatives with respect to global coordinates dNdX, and the Jacobian matrix J 
%   for the given element at the specified local coordinates (XI, ETA).
%
%   Syntax:
%   [N, dNdX, J] = QUADRATICQUADRILATERAL(X, Y, XI, ETA)
%
%   Inputs:
%   X - X-coordinates of the element nodes (numeric vector of length 9)
%   Y - Y-coordinates of the element nodes (numeric vector of length 9)
%   XI - Local coordinate in the xi direction (scalar)
%   ETA - Local coordinate in the eta direction (scalar)
%
%   Outputs:
%   N - Shape functions evaluated at (XI, ETA) (1x9 numeric vector)
%   dNdX - Derivatives of the shape functions with respect to global coordinates (2x9 numeric matrix)
%   J - Jacobian matrix evaluated at (XI, ETA) (2x2 numeric matrix)
%
%   Example:
%   x = [0, 1, 1, 0, 0.5, 1, 0.5, 0, 0.5];
%   y = [0, 0, 1, 1, 0, 0.5, 1, 0.5, 0.5];
%   xi = 0.5;
%   eta = -0.5;
%   [N, dNdX, J] = quadraticQuadrilateral(x, y, xi, eta);
%
%   See also: ELEMENTDIFFUSIONMATRIX, LINEARQUADRILATERAL, LINEARTRIANGLE

    % Define the shape functions for a 9-node Lagrange quadratic quadrilateral
    N = [xi*eta*(xi-1)*(eta-1)/4, xi*eta*(xi+1)*(eta-1)/4, ...
         xi*eta*(xi+1)*(eta+1)/4, xi*eta*(xi-1)*(eta+1)/4, ...
         (1-xi^2)*eta*(eta-1)/2, xi*(xi+1)*(1-eta^2)/2, ...
         (1-xi^2)*eta*(eta+1)/2, xi*(xi-1)*(1-eta^2)/2, ...
         (1-xi^2)*(1-eta^2)];
    
    % Derivatives of the shape functions with respect to local coordinates (xi, eta)
    dNdxi = [eta*(2*xi - 1)*(eta - 1)/4, eta*(2*xi + 1)*(eta - 1)/4, ...
             eta*(2*xi + 1)*(eta + 1)/4, eta*(2*xi - 1)*(eta + 1)/4, ...
             -xi*eta*(eta - 1), (1 - eta^2)*(2*xi + 1)/2, ...
             -xi*eta*(eta + 1), (1 - eta^2)*(2*xi - 1)/2, ...
             -2*xi*(1 - eta^2)];
         
    dNdeta = [xi*(xi - 1)*(2*eta - 1)/4, xi*(xi + 1)*(2*eta - 1)/4, ...
              xi*(xi + 1)*(2*eta + 1)/4, xi*(xi - 1)*(2*eta + 1)/4, ...
              (1 - xi^2)*(2*eta - 1)/2, -xi*(xi + 1)*eta, ...
              (1 - xi^2)*(2*eta + 1)/2, -xi*(xi - 1)*eta, ...
              -2*(1 - xi^2)*eta];
    
    % Initialize the Jacobian matrix
    J = zeros(2, 2);
    
    % Calculate the Jacobian matrix
    for i = 1:9
        J(1, :) = J(1, :) + dNdxi(i)*[x(i), y(i)];
        J(2, :) = J(2, :) + dNdeta(i)*[x(i), y(i)];
    end
    
    % Compute derivatives of N with respect to global coordinates
    dNdX = J \ [dNdxi; dNdeta];
    
end