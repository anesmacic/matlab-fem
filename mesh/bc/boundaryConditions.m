classdef boundaryConditions < handle
    properties
        id
        nodes
        u boundaryConditionParameters
        v boundaryConditionParameters
        p boundaryConditionParameters
    end

    methods
        function obj = boundaryConditions(id, nodes)
            obj.id = id;
            obj.nodes = nodes;
            obj.u = boundaryConditionParameters;
            obj.v = boundaryConditionParameters;
            obj.p = boundaryConditionParameters;
        end
        
        function prescribeConstantBC(obj, dof, value)
            arguments
                obj boundaryConditions
                dof {mustBeMember(dof, {'u','v','p'})}
                value double
            end
            obj.(dof) = boundaryConditionParameters(constant = value);
        end

        function prescribeFunctionBC(obj, dof, fcn)
            arguments
                obj boundaryConditions
                dof {mustBeMember(dof, {'u','v','p'})}
                fcn function_handle
            end
            obj.(dof) = boundaryConditionParameters(fcn = fcn);
        end

        function valid = isValid(obj)
            valid = sum([obj.u.isempty, obj.v.isempty, obj.p.isempty]) == 1;
        end


    end

end