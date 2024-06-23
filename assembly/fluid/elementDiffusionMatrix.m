function D = elementDiffusionMatrix(x, y, xi, eta, xi_w, eta_w, element, N)
% ELEMENTDIFFUSIONMATRIX 
% 
%   Computes the element diffusion matrix. D = ∇φ · ∇φ
%   
%   D = ELEMENTDIFFUSIONMATRIX(X, Y, XI, ETA, XI_W, ETA_W, ELEMENT, N)
%   returns the diffusion matrix D for the given element using the provided
%   quadrature points and weights.
%
%   Syntax:
%   D = ELEMENTDIFFUSIONMATRIX(X, Y, XI, ETA, XI_W, ETA_W, ELEMENT, N)
%
%   Inputs:
%   X - X-coordinates of the element nodes (numeric vector)
%   Y - Y-coordinates of the element nodes (numeric vector)
%   XI - Quadrature points in the xi direction (numeric vector)
%   ETA - Quadrature points in the eta direction (numeric vector)
%   XI_W - Quadrature weights in the xi direction (numeric vector)
%   ETA_W - Quadrature weights in the eta direction (numeric vector)
%   ELEMENT - Handle to the element function that returns shape functions
%   and their derivatives (see elements/*/*)
%   N - Number of shape functions (integer)
%
%   Outputs:
%   D - Element diffusion matrix (N x N numeric matrix)
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
%   D = elementDiffusionMatrix(x, y, xi, eta, xi_w, eta_w, element, N);
%
%   See also: linearQuadrilateral

    % Initialize the diffusion matrix
    D = zeros(N, N);

    % Loop over quadrature points
    for i = 1:numel(xi)
        for j = 1:numel(eta)
            [~, dNdX, J] = element(x, y, xi(i), eta(j));
            D = D + det(J) * xi_w(i) * eta_w(j) * ( ...
                  transpose(dNdX(1,:)) * dNdX(1,:) ...
                + transpose(dNdX(2,:)) * dNdX(2,:));
        end
    end
end