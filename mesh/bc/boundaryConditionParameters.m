classdef boundaryConditionParameters < handle
    properties
        type boundaryConditionType
        value 
        fcn
    end

    methods
        function obj = boundaryConditionParameters(options)
            arguments
                options.constant
                options.fcn
            end
            if isempty(fieldnames(options))
                obj.type = boundaryConditionType.None;
            else
                fn = fieldnames(options);
                switch fn{1}
                    case 'constant'
                        assert(isa(options.constant,'double'));
                        obj.type = boundaryConditionType.Constant;
                        obj.value = options.constant;
                    case 'fcn'
                        assert(isa(options.fcn, 'function_handle'));
                        obj.type = boundaryConditionType.Func;
                        obj.fcn = options.fcn;
                end
            end

        end
    
        function result = isempty(obj)
            result = obj.type == boundaryConditionType.None;
        end
    end
end