classdef User < handle
    
    properties
        Loads
        Meter
    end
    
    properties (Dependent)
        %PVCapacity
        Pl
        CompletedValue % Current cumulative completed value of all activities; state dependent
        InterruptedCost % Current interruption cost of all activities; state dependent
        TotalUtility % Current utility from activities. Equal to completed value - interrupted cost
        ActivityInterruptions % Number of activity interruptions
        LoadNames % Names of loads
        Pa % Power by load
    end
    
    properties (SetAccess = immutable)
        NBatteries % Number of batteries
        NPV % Number of PV modules
        BatteryUnitCapacity % Capacity of each battery, Wh
        PVUnitCapacity % Capacity of each unit, W
    end
    
    properties (SetAccess = private)
        %x
        %ThermalLoads
        Activities = [];
        
        MILPOpts = optimoptions('intlinprog','Display','off');
    end
    
    properties (Access = private)
        Type
        % For performance, keep references for activities that can be updated, sorted by start time.
        queuedActivities
        inProgressActivities = []
    end
    
    methods
        function self = User(userType,nBatteries,batteryUnitCapacity,nPV,pvUnitCapacity)
            import MicrogridDispatchSimulator.Models.Load
            self.Type = userType;
            for i = 1:length(self.Type.ActivityTypes)
                loads(i,1) = Load(self.Type.ActivityTypes(i).LoadType);
            end
            self.Loads = loads;
            
            if (nargin > 1)
                self.NBatteries = nBatteries;
                self.BatteryUnitCapacity = batteryUnitCapacity;
                self.NPV = nPV;
                self.PVUnitCapacity = pvUnitCapacity;
            end
        end
        
        function initializeState(self)
            % Initialize all activities
            for a = self.Activities'
                a.initializeState();
            end
        end
        
        function createActivities(self,T,deltaT)
            % TODO: make this accept a start time
            nextInd = 1;
            for i = 1:length(self.Type.ActivityTypes)
                newActs = self.Type.ActivityTypes(i).createActivities(self.Loads(i),T,deltaT);
                if (~isempty(newActs))
                    activities(nextInd:nextInd+length(newActs)-1,1) = newActs;
                    nextInd = length(activities)+1;
                end
            end
            if (nextInd > 1)
                self.Activities = activities;
                % Sort by start time
                [startTimes, sortedInd] = sort([activities.StartTime]');
                self.queuedActivities = [sortedInd,startTimes];
            else
                warning('No activities created');
            end
        end
        
        function Pl = get.Pl(self)
            %Pl = sum([self.Loads.P]);
            Pl = 0;
            for l = 1:length(self.Loads)
                Pl = Pl+self.Loads(l).P;
            end
        end
        
        function Pa = get.Pa(self)
            Pa = nan(length(self.Loads),1);
            for i = 1:length(self.Loads)
                Pa(i) = self.Loads(i).P;
            end
        end
        
        function y = get.LoadNames(self)
            y = cell(1,length(self.Loads));
            for i = 1:length(self.Loads)
                y{i} = self.Loads(i).Type.Name;
            end
        end
        
        function y = get.CompletedValue(self)
            y = sum([self.Activities([self.Activities.Status] == 2).CompletionValue]);
        end
        
        function y = get.InterruptedCost(self)
            y = sum([self.Activities([self.Activities.Status] == 3).InterruptionCost]);
        end
        
        function y = get.ActivityInterruptions(self)
            y = sum([self.Activities.Status] == 3);
        end
        
        function y = get.TotalUtility(self)
            y = self.CompletedValue - self.InterruptedCost;
        end
        
        function set.Loads(self,loads)
            self.Loads = loads;
            self.Meter.LoadBus.Loads = loads;
        end
        
        function set.Meter(self,meter)
            self.Meter = meter;
            self.Meter.LoadBus.Loads = self.Loads;
        end
        
        function sendLoadLimit(self,t,u)
            %If limit is below infinity, adjust activities
            if isfield(u,'energyLimit') && u.energyLimit < inf
                self.adjustToLimit(t,u.energyLimit,u.T);
            end
        end
        
        function update(self,t,deltaT)
            % Could return information about activity interruptions and/or
            % costs and value from activities
            
            % Update all activities. Using explicit loop instead of
            % arrayfun because it timed faster.
            % Simplest would be to loop over all activities, but optimize
            % to only loop over ones that can be updated
            
            % Update in progress activities
            
            % Get indices of activities that are in progress or scheduled
            % to start in the time window.
            removeInd = false(length(self.inProgressActivities),1);
            for i = 1:length(self.inProgressActivities)
                a = self.Activities(self.inProgressActivities(i));
                a.update(t,deltaT);
                if (a.Status ~= 1)
                    % Remove from list of in progress activities
                    removeInd(i) = true;
                end
            end
            if (any(removeInd))
                self.inProgressActivities = self.inProgressActivities(~removeInd);
            end
            
            
            % Update activities that could start
            
            L = size(self.queuedActivities,1);
            i = 1;
            % Queued activities are sorted by start time, so increment
            % until found the one that is higher.
            while (i <= L && self.queuedActivities(i,2) < t + deltaT)
                a = self.Activities(self.queuedActivities(i,1));
                a.update(t,deltaT);
                if (a.Status == 0)
                    error('this should not be happening');
                elseif (a.Status == 1)
                    self.inProgressActivities = [self.inProgressActivities;self.queuedActivities(i,1)];
                end
                i = i+1;
            end
            if (i > 1)
                self.queuedActivities = self.queuedActivities(i:end,:); % Only keep activities that are still queued
            end
            
        end
        
        function adjustToLimit(self,t,energyLimit,T)
            % t is global time in seconds
            % energyLimit: consumption limit in Wh over horizon
            % T: horizon in seconds
            
            if (isinf(energyLimit) || T == 0)
                % There is no limit
                return;
            end
            
            % Queued activities (activities that are not started but scheduled) are sorted by start time, so increment
            % until found the one that is higher.
            i = 1;
            while (i <= size(self.queuedActivities,1) && self.queuedActivities(i,2) < t + T)
                i = i+1;
            end
            a0Ind = self.queuedActivities(1:i-1,1); % Indices of activities that are scheduled to start within the horizon
            a1Ind = self.inProgressActivities; % Activities that are in progress
            
            if (isempty(a0Ind) && isempty(a1Ind))
                return
            end
            
            % Energy constraint
            if (~isempty(a0Ind))
                A0 = [self.Activities(a0Ind).OnW].*min([self.Activities(a0Ind).RemainingTime],T);
            else
                A0 = [];
            end
            if (~isempty(a1Ind))
                A1 = [self.Activities(a1Ind).OnW].*min([self.Activities(a1Ind).RemainingTime],T);
            else
                A1 = [];
            end
            A = [A0 A1]/3600; % Stack, and convert J to Wh
            
            b = energyLimit;

            % First check if constraint satisfied for no interruption
            if (sum(A) > b)
                % User must make some interruptions

                % Setup MILP
                % Components of objective (maximization form)
                f0 = [self.Activities(a0Ind).CompletionValue].*min(T./[self.Activities(a0Ind).RemainingTime],1);
                f1 = [self.Activities(a1Ind).CompletionValue].*min(T./[self.Activities(a1Ind).RemainingTime],1)...
                    +[self.Activities(a1Ind).InterruptionCost];

                % Compile components into standard form (minimization)
                f = -[f0';f1'];
                intcon=1:size(A,2);
                lb=zeros(size(A));
                ub=ones(size(A));

                % Compute which loads to interrupt.
                [x_opt,~,exitFlag]=intlinprog(f,intcon,A,b,[],[],lb,ub,zeros(size(A))',self.MILPOpts); % x returns activites and thermal loads to interrupt
                if (exitFlag ~= 1)
                    error('intlinprog exited with bad exit flag: %i',exitFlag);
                end

                % Coerce x_opt to binary
                x_opt = abs(x_opt);
                x_opt(x_opt < 1e-6) = 0; % Check for values very close to zero
                x_opt = logical(x_opt);
                
                removeIPInd = false(length(a1Ind),1); % Remove these activities from the stored inprogress
                for i = find(~x_opt(:)') % Iterate over activities to be cancelled
                    if i <= length(a0Ind)
                        self.Activities(a0Ind(i)).cancel();
                    else
                        self.Activities(a1Ind(i-length(a0Ind))).cancel();
                        removeIPInd(i-length(a0Ind)) = true;
                    end
                end
                self.inProgressActivities = self.inProgressActivities(~removeIPInd);
            end
        end
    end
end