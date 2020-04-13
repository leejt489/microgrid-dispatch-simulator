classdef Activity < handle
    properties (SetAccess = private)
        % State variables
        Status % Discrete variable: 0 = not started; 1 = in progress; 2 = completed; 3 = interrupted; 4 = canceled
        RemainingTime
    end
    
    properties (SetAccess = immutable)
        Type % Link to activity type
        % Parameters
        CompletionTime
        CompletionValue
        InterruptionCost
        StartTime
    end
    
    properties (Access = private)
        Load
    end
    
    properties (Dependent)
        OnW
    end
    
    methods
        function self = Activity(startTime,timeToComplete,completionValue,interruptionCost,load,type)
            % Set parameters
            self.StartTime = startTime;
            self.Load = load;
            self.CompletionValue = completionValue;
            self.InterruptionCost = interruptionCost;
            self.CompletionTime = timeToComplete;
            % Link to activity type
            self.Type = type;
            % Initialize state
            self.initializeState();
        end
        
        function initializeState(self)
            self.Status = 0;
            self.RemainingTime = self.CompletionTime;
        end
        
        function update(self,t,deltaT)
            if ((self.Status == 0) && (self.StartTime >= t) && (self.StartTime < t + deltaT)) % Start condition
                self.Load.Connected = true; % Connect the load
                if (self.Load.On)
                    self.Status = 1; % Move to in progress
                else
                    % Move to cancelled
                    self.Load.Connected = false;
                    self.Status = 4;
                end
            end
            if self.Status == 1
                if ~self.Load.On
                    self.cancel(); % Load interrupted, so cancel the activity
                elseif self.RemainingTime == 0
                        % Activity was completed, move to that status
                        self.Status = 2;
                        self.Load.Connected = false;
                else
                    % Decrement remaining time. If remaining time reaches
                    % 0, activity will be moved to completed in next time
                    % step as long as the load remains on.
                    self.RemainingTime = max(self.RemainingTime - deltaT - max((self.StartTime - t),0),0);
                end
            end
                
        end
        
        function cancel(self)
            switch self.Status
                case 0
                    self.Status = 4;
                case 1
                    self.Status = 3;
                    self.Load.Connected = false;
                case 2
                case 3
                case 4
                    error('Cannot cancel a completed, interrupted, or cancelled task');
            end
        end
        
        function y = get.OnW(self)
            y = self.Load.OnW;
        end
    end
end