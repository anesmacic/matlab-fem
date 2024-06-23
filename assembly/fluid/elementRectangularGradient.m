function [Dx, Dy] = elementRectangularGradient(x, y, xi, eta, xi_w, eta_w, elementu, elementp, Nu, Np)
%ELEMENTRECTANGULARGRADIENT Computes the discrete differential operators Dx and Dy for a 2D element.
%   [Dx, Dy] = ELEMENTRECTANGULARGRADIENT(X, Y, XI, ETA, XI_W, ETA_W, ELEMENTU, ELEMENTP, NU, NP)
%   returns the discrete differential operators Dx and Dy for the given element
%   using the provided quadrature points and weights.
%
%   This function evaluates the integral:
%        D_ij ≡ ∫∫_D φ^p_i ∇φ^u_j dx dy
%
%   Syntax:
%   [Dx, Dy] = ELEMENTRECTANGULARGRADIENT(X, Y, XI, ETA, XI_W, ETA_W, ELEMENTU, ELEMENTP, NU, NP)
%
%   Inputs:
%   X - X-coordinates of the element nodes (numeric vector)
%   Y - Y-coordinates of the element nodes (numeric vector)
%   XI - Quadrature points in the xi direction (numeric vector)
%   ETA - Quadrature points in the eta direction (numeric vector)
%   XI_W - Quadrature weights in the xi direction (numeric vector)
%   ETA_W - Quadrature weights in the eta direction (numeric vector)
%   ELEMENTU - Handle to the element function for velocity that returns shape functions and their derivatives
%   ELEMENTP - Handle to the element function for pressure that returns shape functions and their derivatives
%   NU - Number of shape functions for velocity (integer)
%   NP - Number of shape functions for pressure (integer)
%
%   Outputs:
%   Dx - Discrete differential operator in the x direction (Nu x Np numeric matrix)
%   Dy - Discrete differential operator in the y direction (Nu x Np numeric matrix)
%
%   Example:
%   x = [0, 1, 1, 0, 0.5, 1, 0.5, 0];
%   y = [0, 0, 1, 1, 0, 0.5, 1, 0.5];
%   xi = [-1, 1];
%   eta = [-1, 1];
%   xi_w = [1, 1];
%   eta_w = [1, 1];
%   elementu = @serendipityQuadrilateral; % Quadratic element function for u
%   elementp = @linearQuadrilateral; % Linear element function for p
%   Nu = 8; % Number of shape functions for quadratic element
%   Np = 4; % Number of shape functions for linear element
%   [Dx, Dy] = elementRectangularGradient(x, y, xi, eta, xi_w, eta_w, elementu, elementp, Nu, Np);
%
%   See also: ELEMENTDIFFUSIONMATRIX, LINEARQUADRILATERAL, LINEARTRIANGLE, QUADRATICQUADRILATERAL

    % Initialize the differential operators
    Dx = zeros(Nu, Np);
    Dy = zeros(Nu, Np);

    % Loop over quadrature points
    for i = 1:numel(xi)
        for j = 1:numel(eta)
            [~, dNdXu, J] = elementu(x(1:Nu), y(1:Nu), xi(i), eta(j));
            [N, ~, ~] = elementp(x(1:Np), y(1:Np), xi(i), eta(j));
            Dx = Dx + det(J) * xi_w(i) * eta_w(j) * transpose(dNdXu(1,:)) * N;
            Dy = Dy + det(J) * xi_w(i) * eta_w(j) * transpose(dNdXu(2,:)) * N;
        end
    end
end