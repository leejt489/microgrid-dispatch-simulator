classdef Load < handle
    
    properties (Access = public)
        Connected = false % Whether or not the load is connected to supply, i.e. switched on
        Bus % The source of the load. Has a 'V' voltage property.
    end
    
    properties (SetAccess = private)
       OnW
       Type
    end
    
    properties (Dependent)
        On
        P
    end
    
    methods
        function self = Load(loadType)
            self.OnW = loadType.OnW;
            self.Type = loadType;
        end
        
        function y = get.P(self)
            %y = self.On*self.OnW; % Optimizing to remove a function call
            y = ((self.Bus.V > 0) && self.Connected)*self.OnW;
        end
        
        function y = get.On(self)
            y = (self.Bus.V > 0) && self.Connected;
        end
    end
    
end