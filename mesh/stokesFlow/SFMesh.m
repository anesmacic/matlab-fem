classdef SFMesh < handle
%SFMESH Class for managing a structured finite element mesh.
%   The SFMESH class handles the creation and management of a structured finite element mesh,
%   including elements for velocity and pressure, boundary conditions, and mesh connectivity.
%
%   Properties:
%       velocityElement  - Handle to the velocity element function (function_handle)
%       pressureElement  - Handle to the pressure element function (function_handle)
%       nVE              - Number of nodes in the velocity element (integer)
%       nPE              - Number of nodes in the pressure element (integer)
%       bc               - Array of boundary conditions (boundaryConditions array)
%
%       nElements        - Number of elements in the mesh (integer, read-only)
%       nNodes           - Number of nodes in the mesh (integer, read-only)
%       connectivity     - Connectivity matrix (nElements x max(nVE, nPE) integer matrix, read-only)
%       nodes            - Node coordinates (nNodes x 2 numeric matrix, read-only)
%
%   Methods:
%       SFMesh           - Constructor to create an SFMesh object
%       removeUnusedMesh - Removes unused nodes from the mesh
%
%   Syntax:
%   mesh = SFMesh(nodes, connectivity, boundaries);
%   mesh = SFMesh(nodes, connectivity, boundaries, 'velocityElement', @quadraticQuadrilateral, 'pressureElement', @linearTriangle, 'nVE', 9, 'nPE', 3);
%
%   Example:
%   nodes = [0, 0; 1, 0; 1, 1; 0, 1];
%   connectivity = [1, 2, 3, 4];
%   boundaries = {100, [1 2 3 4]}; 
%   mesh = SFMesh(nodes, connectivity, boundaries);
%
%   See also: BOUNDARYCONDITIONS, BOUNDARYCONDITIONTYPE

    properties
        velocityElement  % Handle to the velocity element function (function_handle)
        pressureElement  % Handle to the pressure element function (function_handle)
        nVE % Number of nodes in velocity element (integer)
        nPE % Number of nodes in pressure element (integer)
        bc % Array of boundary conditions (boundaryCondition array)
    end

    properties (SetAccess = private)
        nElements % Number of elements in the mesh (integer, read-only)
        nNodes % Number of nodes in the mesh (integer, read-only)
        connectivity % Connectivity matrix (nElements x max(nVE, nPE) integer matrix, read-only)
        nodes % Node coordinates (nNodes x 2 numeric matrix, read-only)
        idof
        dof2node
        ndof
    end

    methods
        function obj = SFMesh(nodes, connectivity, boundaries, options)
            %SFMESH Constructor to create an SFMesh object.
            %   OBJ = SFMESH(NODES, CONNECTIVITY, BOUNDARIES, OPTIONS) creates an SFMesh object
            %   with the specified nodes, connectivity matrix, and boundary
            %   conditions. Boundary conditions should be provided as a
            %   cell array where each entry is a tuple whose first element
            %   is the boundary identifier, and the second element is array
            %   of node IDs that lie on the boundary
            %   OPTIONS can include:
            %       - 'velocityElement': Handle to the velocity element function (default: @serendipityQuadrilateral)
            %       - 'pressureElement': Handle to the pressure element function (default: @linearQuadrilateral)
            %       - 'nVE': Number of nodes in the velocity element (default: 8)
            %       - 'nPE': Number of nodes in the pressure element (default: 8)
            arguments
                nodes
                connectivity
                boundaries
                options.velocityElement = @serendipityQuadrilateral
                options.pressureElement = @linearQuadrilateral
                options.nVE = 8
                options.nPE = 4
            end
            obj.nodes = nodes;
            obj.connectivity = connectivity;
            obj.nNodes = size(obj.nodes, 1);
            obj.nElements = size(obj.connectivity, 1);
            obj.bc = boundaryConditions.empty;
            for i = 1:size(boundaries, 1)
                obj.bc(i) = boundaryConditions(boundaries{i, 1}, boundaries{i, 2});
            end
            obj.velocityElement = options.velocityElement;
            obj.pressureElement = options.pressureElement;
            obj.nVE = options.nVE;
            obj.nPE = options.nPE;
            obj.removeUnusedMesh;
            pconnectivity = obj.connectivity(:,1:obj.nPE);
            np = numel(unique(pconnectivity(:)));
            nu = obj.nNodes;
            obj.ndof = nu*2 + np;
            obj.idof = zeros(obj.ndof, 3,'logical');
            obj.idof(1     :nu,   1) = 1;
            obj.idof(nu+1  :2*nu, 2) = 1;
            obj.idof(2*nu+1:end,  3) = 1;
            obj.dof2node = zeros(ndof, 1);
            obj.dof2node(1:nu)      = 1:obj.nNodes;
            obj.dof2node(1+nu:2*nu) = 1:obj.nNodes;
            obj.dof2node(2*nu+1:end)= sort(unique(pconnectivity(:)));
        end

        function removeUnusedMesh(obj)
            %REMOVEUNUSEDMESH Removes unused nodes from the mesh.
            %   This method updates the connectivity matrix and node coordinates
            %   to remove any nodes that are not referenced in the connectivity matrix.
            numberOfUsedNodes = numel(unique(obj.connectivity(:)));
            if numberOfUsedNodes > obj.nNodes
                error('SFMesh:InvalidInput', 'Number of nodes in the referenced connectivity matrix exceeds the number of provided nodes.')
            end
            if numberOfUsedNodes < obj.nNodes
                nodeIsReferenced = ismember(1:obj.nNodes, obj.connectivity(:));
                mapping = zeros(obj.nNodes, 1);
                mapping(nodeIsReferenced) = 1:sum(nodeIsReferenced);
                obj.connectivity = mapping(obj.connectivity);
                obj.nodes = obj.nodes(nodeIsReferenced, :);
                obj.nNodes = sum(nodeIsReferenced);
                for i = 1:numel(obj.bc)
                    obj.bc(i).nodes = mapping(obj.bc(i).nodes);
                end
            end
        end

    end
end