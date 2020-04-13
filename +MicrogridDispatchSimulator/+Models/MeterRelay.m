classdef MeterRelay < handle
    properties
        SourceBus
        LoadBus
        
        % State variables
        RelayClosed = true
        EnergyLimitAmount = Inf % Energy limit in Wh
        EnergyLimitTime = 0 % Time remainting on energy limit in seconds
        PowerLimit = Inf % Power limit in W
    end
    
%     properties (Dependent)
%         Pl
%     end
        
    
    methods
        function self = MeterRelay()
            % Create buses for source and load. Can add constructor to take
            % existing buses as an argument, but don't need that.
            import MicrogridDispatchSimulator.Models.Bus
            
            self.LoadBus = Bus();
            self.SourceBus = Bus();
            self.SourceBus.ChildBuses = self.LoadBus;
        end
        
        function update(self,t,deltaT)
            if ~isinf(self.EnergyLimitAmount) || ~isinf(self.PowerLimit)
                Pl = self.LoadBus.Pl;
                self.EnergyLimitAmount = max(self.EnergyLimitAmount - deltaT*Pl/3600,0);
                self.EnergyLimitTime = max(self.EnergyLimitTime - deltaT,0);
                if ((self.EnergyLimitTime > 0 && self.EnergyLimitAmount == 0) || (Pl > self.PowerLimit)) && self.RelayClosed
                    self.RelayClosed = false;
                elseif (~self.RelayClosed)
                    self.RelayClosed = true;
                end
            elseif (~self.RelayClosed)
                self.RelayClosed = false;
            end
        end
        
        function set.RelayClosed(self,x)
            self.RelayClosed = x;
            self.LoadBus.V = x*1;
        end
    end
end