classdef MicrogridController < handle
    properties (Constant)
        MinSoC = 0.1 % Constrain controller to attempt to treat 10% as the minimum allowable SoC
    end
    properties (SetAccess = private)
        RequiresSolarForecast = false
        RequiresLoadForecast = false
        HasMemory = false
    end
    
    properties (Access = private)
        Memory = struct % Use for learning-based controllers
        Microgrid
        LoadControlFun
        Params
    end
    
    properties (Dependent, Access = private)
        E % Energy in storage, adjusted for MinSoC
        Emax % Max energy capacity, adjusted for MinSoC
        E_hor % Energy in storage, adjusted for MinSoC, converted from Wh to W*(deltaT_hor)
        Emax_hor % Max energy capacity, adjusted for MinSoC, converted from Wh to W*(deltaT_hor)
        E_cntr % Energy in storage, adjusted for MinSoC, converted from Wh to W*(deltaT_cntr)
        Emax_cntr % Max energy capacity, adjusted for MinSoc, converted from Wh to W*(deltaT_cntr)
    end
    
    methods
        function self = MicrogridController(microgrid,loadControllerName,params)
            self.Microgrid = microgrid;
            self.Params = params;
            switch loadControllerName
                case 'NoControl'
                    self.LoadControlFun = @self.NoControl;
                case 'Reactive' % Limits defined assuming Pl_max is 10 kW.
                    self.LoadControlFun = @self.Reactive;
                case 'Deterministic'
                    self.LoadControlFun = @self.Deterministic;
                    self.RequiresLoadForecast = true;
                    self.RequiresSolarForecast = true;
                case 'StochasticTrajectory'
                    self.LoadControlFun = @self.StochasticTrajectory;
                    self.RequiresLoadForecast = true;
                    self.RequiresSolarForecast = true;
                case 'StochasticDP'
                    self.LoadControlFun = @self.StochasticDP;
                    self.RequiresLoadForecast = true;
                    self.RequiresSolarForecast = true;
                otherwise
                    error('Unknown algorithm name');
            end
        end
        
        function y = update(self,t,u)
            % u can include past outputs from the environment and also forecast data
            
            if (self.RequiresLoadForecast || self.RequiresSolarForecast)
                [l,f] = self.LoadControlFun(u.W);
            else
                [l,f] = self.LoadControlFun();
            end
            
            % 'l' will be average power over the time horizon in W
            energyLimit = l*self.Params.deltaT_cntr/3600; % Convert average power over period deltaT to Wh, deltaT in seconds
            loadLimit.T = self.Params.deltaT_cntr;
            
            for n = 1:length(self.Microgrid.Users)
                % Send limit to user so they can adjust activities
                loadLimit.energyLimit = energyLimit(n);
                self.Microgrid.Users(n).sendLoadLimit(t,loadLimit);
                % Send limit to the meter so it can enforce
                self.Microgrid.Users(n).Meter.EnergyLimitAmount = loadLimit.energyLimit;
                self.Microgrid.Users(n).Meter.EnergyLimitTime = loadLimit.T;
            end
            
            % Compute power injection (can move this to another function if
            % more computation desired.
            netSOC = sum(self.E)/sum(self.Emax); % Average SOC of entire system
            Pinj = (self.E-netSOC*self.Emax)*self.Params.K/self.Params.deltaT_cntr*3600;
            self.Microgrid.Pinj = Pinj; % Update it on the microgrid
            
            y.Pinj = Pinj;
            loadLimit.energyLimit = energyLimit;
            y.loadLimit = loadLimit;
            y.objPred = f;
        end
        
        function E = get.E(self)
            E = max(self.Microgrid.E-self.MinSoC*self.Microgrid.Emax,0); % Define SOC relative to 10%
        end
        
        function Emax = get.Emax(self)
            Emax = (1-self.MinSoC)*self.Microgrid.Emax; % Define max relative to 10% 
        end
        
        function E = get.E_hor(self)
            E = self.E/self.Params.deltaT_hor*3600;
        end
        
        function Emax = get.Emax_hor(self)
            Emax = self.Emax/self.Params.deltaT_hor*3600;
        end
        
        function E = get.E_cntr(self)
            E = self.E/self.Params.deltaT_cntr*3600;
        end
        
        function Emax = get.Emax_cntr(self)
            Emax = self.Emax/self.Params.deltaT_cntr*3600;
        end
    end
    
    methods (Access = private)
        function [l,f] = NoControl(self)
            l = inf(length(self.Microgrid.E),1);
            f = nan;
        end
        
        function [l,f] = Reactive(self)
            netSOC = sum(self.E)/sum(self.Emax);
            if (netSOC >= 0.3)
                l = inf(length(self.E),1);
            elseif (netSOC >= 0.2)
                l = 0.1*self.Microgrid.Plmax;
            elseif (netSOC >= 0.1)
                l = 0.05*self.Microgrid.Plmax;
            else
                l = 0.01*self.Microgrid.Plmax;
            end
            f = nan;
        end
        
        function [l,f] = Deterministic(self,forecast)
            [l,f] = MicrogridDispatchSimulator.DispatchControllers.scenarioFormulation2Stage1LookGurobi(...
                self.Params.gamma,mean(forecast.Pg,3),mean(forecast.Pl,3),1,...
                self.E_hor/1000,self.Emax_hor/1000,self.Microgrid.Pmax/1000,self.Microgrid.Plmax/1000,...
                self.Microgrid.Pcmax/1000 ...
            );
            l = l*1000; % convert to W
        end
        
        function [l,f] = StochasticTrajectory(self,forecast)
            [l,f] = MicrogridDispatchSimulator.DispatchControllers.scenarioFormulation2Stage1LookGurobi(...
                self.Params.gamma,forecast.Pg,forecast.Pl,forecast.ps,...
                self.E_hor/1000,self.Emax_hor/1000,self.Microgrid.Pmax/1000,self.Microgrid.Plmax/1000,...
                self.Microgrid.Pcmax/1000 ...
            );
            l = l*1000; % convert to W
        end
        
        function [l,f] = StochasticDP(self,forecast)
            [l,f] = MicrogridDispatchSimulator.DispatchControllers.scenarioFormulationBackRecGurobi(...
                self.Params.gamma,forecast.Pg,forecast.Pl,forecast.ps,...
                self.E_hor/1000,self.Emax_hor/1000,self.Microgrid.Pmax/1000,self.Microgrid.Plmax/1000,...
                self.Microgrid.Pcmax/1000 ...
            );
            l = l*1000;
        end
    end
end