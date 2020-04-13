classdef SimulationOutputs < handle
    properties (SetAccess = private)
        % Time series (simulation time)
        E % Energy stored
        P % Power injected; positive is into the grid
        Pc % Charge power into battery
        Pl % Load power consumed
        Pg % Generation power
        Pw % Solar power wasted / curtailed
        xi % Blackout
        %chi
        %ci
        %l2
        deltaOmega % Change in frequency
        timeSim
        timeCntr
        
        % Time series (control time)
        Pinj % Power injection setpoint
        l % Load limit (energy)
        T_l % Time horizon of load limit
        objPred % Predicted objective
        
        % Final state
        Users
        Microgrid
    end
    
    properties (SetAccess = private)
        DeltaT_sim
        DeltaT_cntr
    end
    
    methods
        function self = SimulationOutputs(microgrid,T,deltaT_cntr,deltaT_sim)
            self.Microgrid = microgrid;
            self.Users = microgrid.Users;
            N = length(self.Users);
            self.DeltaT_sim = deltaT_sim;
            self.DeltaT_cntr = deltaT_cntr;
            
            % 'T' is the number of time steps in the simulation
            I_sim = T;
            I_cntr = T*deltaT_sim/deltaT_cntr;
            
            % Time
            self.timeSim = (0:I_sim-1)'*deltaT_sim;
            self.timeCntr = (0:I_cntr-1)'*deltaT_cntr;
            
            % Simulation outputs by user
            self.E = nan(N,I_sim); % Battery state
            self.Pl = nan(N,I_sim); % Power actually used
            self.Pw = nan(N,I_sim); % Power wasted (i.e. solar that could not be used)
            self.P = nan(N,I_sim); % Actual net injection
            self.Pg = nan(N,I_sim); % PV generated
            self.xi = nan(I_sim,1); % Blackout
            %self.chi = nan(N,I_sim); % Number of interruptions
            %self.ci = nan(N,I_sim); % Cost of interruption
            %self.l2 = nan(N,I_sim); % Load limit (potentially modified by system after blackout)
            
            % Simulation outputs by system
            self.xi = nan(I_sim,1);
            self.deltaOmega = nan(I_sim,1);
            
            % Control outputs
            self.l = nan(N,I_cntr); % Load limit (control sent)
            self.Pinj = nan(N,I_cntr); % Injection setpoint
            self.objPred = nan(I_cntr,1);
        end
        
        function saveMG(self,t,yMG)
            % t is the absolute simulation time in seconds
            i = t/self.DeltaT_sim+1;
            
            self.E(:,i) = yMG.E;
            self.P(:,i) = yMG.P;
            self.Pg(:,i) = yMG.Pg;
            self.Pw(:,i) = yMG.Pw;
            self.Pl(:,i) = yMG.Pl;
            self.Pc(:,i) = yMG.Pc;
            
            self.xi(i) = yMG.blackout;
            self.deltaOmega(i) = yMG.deltaOmega;
        end
        
%         function saveUser(self,n,t,yU)
%             i = t/self.DeltaT_sim; % t is the absolute simulation time in seconds
%             self.Pu(n,i) = yU.Pu;
%         end
        
        function saveMGC(self,t,yMGC)
            i = t/self.DeltaT_cntr+1; % t is the absolute simulation time in seconds
            self.l(:,i) = yMGC.loadLimit.energyLimit;
            % self.T(i) = yMGC.loadLimit.T; % Don't save because same for
            % all (equal to DeltaT_cntr) as of current formulation)
            self.Pinj(:,i) = yMGC.Pinj;
            self.objPred(i) = yMGC.objPred;
        end
    end
end