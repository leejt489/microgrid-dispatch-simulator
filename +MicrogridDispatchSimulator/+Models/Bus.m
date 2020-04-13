classdef Bus < handle
    properties
        V = 1
        Loads
        ChildBuses = [];
        PVCapacity
        BatteryEnergyCapacity % Storage capacity connected to this bus, Wh
        %BatteryPowerCapacity % Maximum charge/discharge power, positive, W. Minimum assumed to be negative of maximum.
        InjectionMax = inf % W Maximum power injection, positive, W. Minimum assumed to be negative of maximum for now
        LoadMax % Maximum load power (rating of load), W.
    end
    
    properties (Dependent)
        Pl % Power injected
    end
    
    methods
        function Pl = get.Pl(self)
            % Compute self load power
            Pl = 0;
            for l = 1:length(self.Loads)
                Pl = Pl+self.Loads(l).P;
            end
            % Plus children load power
            for c = 1:length(self.ChildBuses)
                Pl = Pl+self.ChildBuses(c).getPl();
            end
        end
        
        function Pl = getPl(self)
            % Same as get.Pl above, but defining it separately so it can be
            % called from get.Pl.
            Pl = 0;
            for l = 1:length(self.Loads)
                Pl = Pl+self.Loads(l).P;
            end
            for c = 1:length(self.ChildBuses)
                Pl = Pl+self.ChildBuses(c).getPl();
            end
        end
        
        function set.Loads(self,loads)
            self.Loads = loads;
            for i = 1:length(self.Loads)
                self.Loads(i).Bus = self;
            end
        end
        
        function set.V(self,v)
            self.V = v;
            for c = 1:length(self.ChildBuses)
                self.ChildBuses(c).setV(v);
            end
        end
        
        function setV(self,v)
            self.V = v;
        end
    end
end