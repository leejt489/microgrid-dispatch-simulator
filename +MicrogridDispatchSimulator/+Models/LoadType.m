classdef LoadType < handle
    properties (SetAccess = immutable)
        Name
        OnW
    end
    
    methods
        function self = LoadType(onW,name)
            self.Name = name;
            self.OnW = onW;
        end
    end
end