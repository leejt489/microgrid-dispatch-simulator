classdef Microgrid < handle
    
    properties (Constant)
        
    end
    
    properties
        Users
        Pinj % Power injection setpoint, can be updated by controller
    end
    
    properties (SetAccess = private)
        E
    end
    
    properties (Access = private)
        Buses
        %Areas = []
        Blackout = false
        
        % Parameters
        ERestart % Fraction of SOC to reconnect after blackout
        BatteryChargeRate % Nominal maximum charge / discharge rate (hours)
        UserMaxLoad
        BusMaxInjection
        BetaPU % DER area tie-line flow frequency response stiffness, p.u. power per Hz
    end
    
    properties (Dependent)
        Emax % Nominal storage capacity per user bus, Wh, Nx1 vector
        Pmax % Nominal max power transfer at connection. Could be same as max load, or different if load meterd separately. Should be defined at network level.
        Plmax % Nominal maximum possible load of each user, W, Nx1 vector. Could define this on the user. Should be thought of as meter or load power supply rating
        Pcmax % Nominal maximum battery charge/rate ber user bus, W, Nx1 vector.
        Beta % DER area tie-line flow frequency response stiffness per user bus, W/Hz, Nx1 vector
    end
    
    methods
        function self = Microgrid(params)
            %self.Buses = Bus();
            %self.Buses.V = 1;
            self.UserMaxLoad = params.UserMaxLoad;
            self.ERestart = params.ERestart;
            self.BatteryChargeRate = params.BatteryChargeRate;
            self.BusMaxInjection = params.BusMaxInjection;
            self.BetaPU = params.Beta;
        end
        
        function initialize(self,initialSoC)
            % Set initial state.
            % initialSoC can be a scalar or Nx1 vector in interval [0,1] of
            % the initial state of charge of each DER battery
            
            if (nargin < 2)
                initialSoC = 0.7; % Default initial state of charge
            end
            self.E = [self.Emax].*initialSoC; % Initial SOC
            for i = 1:length(self.Buses)
                self.Buses(i).V = 1;
            end
        end
        
        function set.Users(self,users)
            import MicrogridDispatchSimulator.Models.MeterRelay
            
            for i = 1:length(users)
                user = users(i);
                
                % This keeps 1-1-1 bus, area, and user index. WARNING: need to change if
                % start using buses and areas differently
                
                % Create a meter with a source bus and load bus and assign
                % to user
                m = MeterRelay();
                b = m.SourceBus;
                user.Meter = m;

                % Associate user DER information directly to the bus
                b.BatteryEnergyCapacity = user.NBatteries*user.BatteryUnitCapacity;
                b.PVCapacity = user.NPV*user.PVUnitCapacity;
                b.LoadMax = self.UserMaxLoad;
                b.InjectionMax = self.BusMaxInjection; % Could define this differently for users
                buses(i,1) = b;
            end
            
            self.Users = users;
            self.Buses = buses;
        end
        
        function y = update(self,deltaT,u)

            pSet = self.Pinj; % Injection setpoint in W
            PVMaxGen = u.irradiance*[self.Buses.PVCapacity]'; % Solar generation in W
            Pl = [self.Buses.Pl]'; % W
            
            
            % Step 1
            % Calculate the maximum tie-line limits for this user's DER
            % system and clamp the control setpoint if needed
            E2 = max(self.E,0); % Treat negative charge as zero.
            battSOC=E2./self.Emax;
            
            % Calculate battery injection and withdrawal limits for each battery given current SOC
            maxBattOut=min(E2/deltaT*3600,self.Pcmax.*(...
                (battSOC/0.1-1).*(battSOC < 0.1)+1 ... % When battSOC < 10%, reduces max output to zero
                ));
            maxBattIn = min((self.Emax-E2)/deltaT*3600,self.Pcmax.*(...
                ((1-battSOC)/0.1-1).*(battSOC > 0.9)+1 ... % When battSOC > 90%, reduces max input to zero
                ));

            % maximum injection and withdrawal
            der_max_P=PVMaxGen+maxBattOut; 
            der_min_P=-maxBattIn;

            % correct setpoint
            pSet=max(der_min_P,min(der_max_P,pSet));
            uGBalance=sum(pSet-Pl); % net injection from all the users

            % Step 2
            % Calculate frequency deviation. If nan returned, no valid
            % frequency
            dF=calc_frequency(pSet,uGBalance,der_max_P,der_min_P,self.Beta); % frequency deviation and blackout condition. dF = nan is the implicit blackout condition
            
            if (isnan(dF))
                % This registers a blackout DURING the current time period
                % by being unable to provide power demand. Results in
                % clearing power. SoC does not change.
                self.Blackout = true; % This effectively disconnects all load
                der_P = zeros(size(self.E));
                PVMaxGen = zeros(size(self.E));
                pBatt = zeros(size(self.E));
                Pl = zeros(size(self.E));
            else
                der_P=max(der_min_P,min(der_max_P,pSet-dF*self.Beta)); % Compute DER output from frequency deviation
                pBatt=min(maxBattOut,max(-maxBattIn,der_P-PVMaxGen)); % results in the least PV curtailment; i.e. assume re-dispatched string to prioritize PV
                self.E = self.E - pBatt*deltaT/3600; % Next state of charge. Computation of pBatt should ensure it is within bounds
                
                % Update blackout state
                if (self.Blackout && sum(self.E)/sum(self.Emax) >= self.ERestart)
                    % Restart and leave blackout state.
                    % TODO: implement proper black start with load control
                    self.Blackout = false;
                end
            end
            
            y.Pl = Pl;
            y.P = y.Pl-der_P;
            y.Pg = der_P - pBatt;
            y.Pw = PVMaxGen - y.Pg;
            y.Pc = -pBatt;
            y.E = self.E;
            y.blackout = self.Blackout;
            y.deltaOmega = dF;
            
        end
        
        function y = get.Beta(self)
            y = self.BetaPU*(self.Pcmax+[self.Buses.PVCapacity]');
        end
        
        function y = get.Emax(self)
            y = [self.Buses.BatteryEnergyCapacity]';
        end
        
        function y = get.Pcmax(self)
            y = self.Emax*self.BatteryChargeRate;
        end
        
        function y = get.Plmax(self)
            y = [self.Buses.LoadMax]';
        end
        
        function y = get.Pmax(self)
            y = [self.Buses.InjectionMax]';
        end
        
        function set.Blackout(self,val)
            self.Blackout = logical(val);
            for n = 1:length(self.Buses)
                self.Buses(n).V = (~self.Blackout)*1;
            end
        end
                
    end
end

function dF = calc_frequency(pSet,uGBalance,der_max_P,der_min_P,beta)
% Returns frequency deviation dF, or nan if no solution given power bounds.
% The deviation is computed by finding the steady state equilibrium of
% frequency and power injection given the setpoints, power imbalance,
% minimum and maximum outputs of each DER, and the stiffness

if uGBalance == 0
    % No imbalance so no frequency deviation
    dF = 0;
else
    % Set up variables for search
    n = 1;
    dP_curr = 0;
    dF_curr = 0;
    dF = nan;
    
    if uGBalance > 0
        % Injections exceed withdrawals -> need to increase the frequency
        
        % First check if a solution can't be found
        if (-sum(der_min_P) < uGBalance)
            return;
        end
        
        % Upward frequency response curve, used to reduce injections from
        % DER
        der_max_dF=(pSet-der_min_P)./beta; % When the frequency increases, the injection from DER decreases. Calculate the frequency at which each DER can no longer increases output
        can_reduce = der_max_dF > 0; % Index of infinity systems that can reduce output
        der_max_dF=der_max_dF(can_reduce);

        beta_up=beta(can_reduce); % responses to overfrequency
        N=length(beta_up);

        % Search over units to find which are at their minimum output
        [max_dF,ind]=sort(der_max_dF);
        n = 1;
        dP_curr = 0;
        dF_curr = 0;
        while (true)
            if (n > N)
                % Didn't find solution; break
                break;
            end

            beta_int = sum(beta_up(ind(n:end)));
            dP_next = dP_curr - (max_dF(n)-dF_curr)*beta_int;
            if (dP_next <= -uGBalance)
                % Frequency and power will be in this interval
                dF = dF_curr - (-uGBalance - dP_curr)/beta_int;
                break;
            end
            % Skip iterations where max_dF is the same
            while (n < N && max_dF(n+1) == dF_curr)
                n = n+1;
            end

            % Update variables for next iteration
            dF_curr = max_dF(n);
            dP_curr = dP_next;
            n = n+1;
        end

    else
        % Withdrawals exceed injections -> need to decrease the frequency
        
        % First check if a solution can't be found
        if (-sum(der_max_P) > uGBalance)
            return;
        end
        
        % Downward frequency response curve, used to increase injections from DER
        der_min_dF=-(der_max_P-pSet)./beta; % % When the frequency decreases, the injection from DER increases. Calculate the frequency at which each DER can no longer increase output
        can_increase = der_min_dF < 0; % Index of DER systems that can increase output
        der_min_dF=der_min_dF(can_increase);

        beta_down=beta(can_increase); % responses to underfrequency
        N=length(beta_down);

        [min_dF,ind]=sort(-der_min_dF);
        min_dF = -min_dF;
        while (true)
            if (n > N)
                % Didn't find solution; break
                break;
            end

            beta_int = sum(beta_down(ind(n:end)));
            dP_next = dP_curr - (min_dF(n)-dF_curr)*beta_int;
            if (dP_next >= -uGBalance)
                % Frequency and power will be in this interval
                dF = dF_curr - (-uGBalance - dP_curr)/beta_int;
                break;
            end
            % Skip iterations where max_dF is the same
            while (n < N && min_dF(n+1) == dF_curr)
                n = n+1;
            end

            % Update variables for next iteration
            dF_curr = min_dF(n);
            dP_curr = dP_next;
            n = n+1;
        end
    end
end

end