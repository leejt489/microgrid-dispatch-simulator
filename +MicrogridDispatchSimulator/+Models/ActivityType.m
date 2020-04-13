classdef ActivityType < handle
    properties (SetAccess = immutable)
        %TimeToComplete
        LoadType
        CompletionValue
        InterruptionCost
        TimeToCompleteMin % seconds
        TimeToCompleteMax % seconds
        StartHourProbability
        %TODO: add time dependent value
        %DeltaT
    end
    
    methods
        function self = ActivityType(loadType,completionValue,interruptionCost,timeToCompleteMin,timeToCompleteMax,startHourProbability)
            self.LoadType = loadType;
            self.CompletionValue = completionValue;
            self.InterruptionCost = interruptionCost;
            self.TimeToCompleteMin = timeToCompleteMin;
            self.TimeToCompleteMax = timeToCompleteMax;
            self.StartHourProbability = startHourProbability;
        end
        
        function activities = createActivities(self,load,T,deltaT)
            import MicrogridDispatchSimulator.Models.Activity
            
            % Creates a vector of activities, non-overlapping
            if (self.TimeToCompleteMin < deltaT)
                warning('ActivityType for LoadType %s has a TimeToCompleteMin of %i, which is less than the deltaT of %i',...
                    self.LoadType.Name, self.TimeToCompleteMin, deltaT);
            end
            actInd = 1;
            endTime = 0; % State, end time
            r = rand(ceil(T/deltaT),1);
            for i = 0:floor(T/deltaT)-1
                t = deltaT*i; % Time
                if (t >= endTime && r(i+1) <= self.startProbability(t,deltaT))
                    % First clause checks that we are not overlapping with an existing activity,
                    % second is random start. If both, create a new
                    % activity starting at that time.
                    
                    % Random time to complete in interval. Round up to be a
                    % multiple of deltaT.
                    timeToComplete = self.TimeToCompleteMin + (self.TimeToCompleteMax-self.TimeToCompleteMin)*rand;
                    timeToComplete = ceil(timeToComplete/deltaT)*deltaT;
                    
                    activities(actInd,1) = Activity(t,timeToComplete,...
                        self.CompletionValue,self.InterruptionCost,load,...
                        self);
                    actInd = actInd + 1;
                    
                    endTime = t + timeToComplete;
                end
            end
            if ~exist('activities','var')
                activities = [];
            end
        end
        
        function p = startProbability(self,t,deltaT)
            % Get day type
            dayOfWeek = ceil(mod(t,3600*24));
            if dayOfWeek == 1
                d = 3; % Sunday
            elseif dayOfWeek == 7
                d = 2; % Saturday
            else
                d = 1; % Weekday
            end
            h = fix(mod(t/3600,24))+1; % Hour of day, 1 to 24
            ph = self.StartHourProbability(h,d);
            p = 1-(1-ph)^(deltaT/3600);
        end
    end
end