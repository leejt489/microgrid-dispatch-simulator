function [v,f,tSolve] = scenarioFormulationBackRec(gamma,Pg,Pl,ps,Estor0,E_max,P_max,Pl_max,Pc_max)%,Pt_max,Bbus,Mbus)

ts = tic;
try
    cvx_solver('gurobi');
catch
    warning('Could not set solver to gurobi');
end

cvx_solver_settings('TimeLimit',1800) % Set a 30 minute time limit on any run

[N, T, S] = size(Pg);
% Comment in for power flow constraints
% K = size(Mbus,1);

% Check dimensions of inputs
if (any(size(Pg) ~= size(Pl)))
    error('Size of solar and load forecast matrices should be identical');
end
if (~(isvector(ps) && length(ps) == S))
    error('ps should be a vector with length equal to number of scenarios %i',S);
end
if (any(size(Estor0) ~= [N 1]))
    error('x0 should be Nx1');
end
if (any(size(E_max) ~= [N 1]))
    error('x_max should be Nx1');
end
if (any(size(Pl_max) ~= [N 1]))
    error('Pl_max must be scalar');
end

% Comment in for Pinj
% Index of users with batteries
%battInd = E_max > 0;
%Nbatt = sum(battInd);
%beta = Pc_max(battInd); % Approximate participation factor in power sharing

% State space discretization. Consider the state to be one aggregated
% battery
nXbins = 10; % Number of samples per state.

% Set future cost variable, initialize to cost of zero
V1 = zeros(nXbins+1,1);
V0 = V1;

t = T;
while (t > 1)
    
    if (S >= 1)
        M = repmat(max(Pl(:,t,:),[],3),1,S); %Static variable representing max possible for each load over all scenarios at time 1
    end
    
    for xBinInd = 0:nXbins % Iterate over possible values at current state
        
        % Choose the optimal input from current state
        cvx_begin
            cvx_quiet true
            % System variables
            variable EstorTot(1,S) nonnegative % Total energy in storage at next time
            variable PcTot(1,S) % Charge power into batteries
            variable Pu(N,S) nonnegative % Power used in each scenario
            variable Pw(N,S) nonnegative % Power wasted
            
            % Auxiliary variables
            if (S >= 1)
                variable v(N,1) nonnegative %'v' is the load limit; 'l' in the paper, but 'l' looks too much like the number 1
                variable q1(N,S) binary %q1(T,S) = 1 => v(T) <= min(Pl(T,S),x(T+S)+Pg(T,S))
                %variable q2(N,S) binary %Pl(T,S) <= min(v(T),x(T+S)+Pg(T,s)
                
                variable r(S,nXbins+1) nonnegative %for piecewise linear objective. 'r' indicates the weights for where you are on the line
                variable y(S,nXbins+1) binary % for piecewise linear objective. 'y' indicates the bin you are in
            end
            
            minimize(...
                Pu(:)'*kron(diag(sparse(ps)),diag(sparse(gamma./Pl_max/2)))*Pu(:)...
                -kron(ps',gamma')*Pu(:)...
                +ps'*r*V1... % Gives a piecewise linear interpolation of the value function V1(x) from the samples calculated
            )

            subject to
                % Power curtailment must be less than generation
                Pw <= reshape(Pg(:,t,:),[N S]);

                % Flow balance at each node
                PcTot == sum(reshape(Pg(:,t,:),[N S])-Pw-Pu,1);
                % Injection limits
                % Charge limits
                PcTot <= sum(Pc_max)*ones(1,S);
                PcTot >= -sum(Pc_max)*ones(1,S);

                % Conservation of energy in batteries
                PcTot == EstorTot-sum(E_max)*xBinInd/nXbins*ones(1,S); % Second term is Estor0 for a current possible initial state

                % Battery capacity
                EstorTot <= sum(E_max)*ones(1,S);

                Pu <= reshape(Pl(:,t,:),[N S]);

                if (S >= 1)
                    % If more than one scenario, decouple the load limit from used
                    % power
                    % The below constraints (in addition with above on state) enforce: Pu == min(min(v,Pl),x+Pg);
                    Pu <= repmat(v,1,S);
                    %q1+q2 == 1;
                    repmat(v,1,S) <= Pu + M.*(1-q1);
                    repmat(v,1,S) <= reshape(Pl(:,t,:),[N S]) + M.*(1-q1);
                    reshape(Pl(:,t,:),[N S]) <= Pu + M.*q1;
                    reshape(Pl(:,t,:),[N S]) <= repmat(v,1,S) + M.*q1;
%                     reshape(Pl(:,t,:),[N S]) <= Pu + M.*(1-q2);
%                     reshape(Pl(:,t,:),[N S]) <= repmat(v,1,S) + M.*(1-q2);
                end
                
                % SOS constraints
                EstorTot == (0:1/nXbins:1)*(r')*sum(E_max);
                r <= y;
                sum(r,2) == ones(S,1);
                sum(y,2) <= 2*ones(S,1);
                for i = 1:nXbins-1
                    for j = i+2:nXbins+1
                        y(:,i) + y(:,j) <= ones(S,1);
                    end
                end
                
        cvx_end
        
        if ~strcmpi(cvx_status,'solved')
            error('CVX could not solve the optimization problem');
        end
        
        V0(xBinInd+1) = cvx_optval;
        
        % We know that the V(x) is monotonically decreasing in x and that
        % it is strictly monotonic up to a point and then constant. So
        % check if it is constant to avoid unnecessarily computing the
        % rest.
        if (xBinInd > 0 && abs(1-V0(xBinInd)/(V0(xBinInd+1)+0.0001)) < 1e-3)
            V0(xBinInd+1:end) = V0(xBinInd);
            break;
        end
    end
    %Plot for debugging to see how V evolves
    %plot((0:1/nXbins:1)*sum(E_max),V0)%,[0 sum(E_max)],-5*N*(T+1-t)*ones(1,2));
    %xlabel('Total Energy in Storage');
    %ylabel(sprintf('Expected Future Cost To Go From Time %i',t));
    
    % Move back a time step
    V1 = V0;
    t = t-1;
end

% Now, given V1, solve the first-step optimization given battery initial
% state


if (S >= 1)
    M = repmat(max(Pl(:,t,:),[],3),1,S); %Static variable representing max possible for each load over all scenarios at time 1
end

cvx_begin quiet
    % Comment in for Pinj
    % variable Pinj(Nbatt,1) % Power injected
    % System variables
    variable Estor(N,S) nonnegative % Total storage after time 1
    variable P(N,S) % Realized net injection
    variable Pu(N,S) nonnegative % Power used in each scenario
    variable Pc(N,S) % Charge power into batteries
    variable Pw(N,S) nonnegative % Power wasted
    % Auxiliary variables
    if (S >= 1)
        variable v(N,1) nonnegative %'v' is the load limit; 'l' in the paper, but 'l' looks too much like the number 1
        variable q1(N,S) binary %q1(T,S) = 1 => v(T) <= min(Pl(T,S),x(T+S)+Pg(T,S))
        variable q2(N,S) binary %Pl(T,S) <= min(v(T),x(T+S)+Pg(T,s)
                
        variable r(S,nXbins+1) nonnegative %for piecewise linear objective. 'r' indicates the weights for where you are on the line
        variable y(S,nXbins+1) binary % for piecewise linear objective. 'y' indicates the bin you are in
    end
    
    % Quadratic cost summed across scenarios.
    minimize(...
        Pu(:)'*kron(diag(ps),diag(gamma./Pl_max/2))*Pu(:)...
        -kron(ps',gamma')*Pu(:)...
        +ps'*r*V1...
    )
    
    subject to
        % Power curtailment must be less than generation
        Pw <= reshape(Pg(:,t,:),[N S]);
        
        % Comment in for Pinj
        % Relate power injection to setpoint by the beta characteristic
        % (participation factor)
        %repmat(Pinj,1,S)-P(battInd,:) == (beta/sum(beta)*ones(1,Nbatt))*(repmat(Pinj,1,S)-P(battInd,:));
        % Flow balance at each node
        P == reshape(Pg(:,t,:),[N S])-Pw-Pu-Pc;
        % All constraints sum to zero
        sum(P) == zeros(1,S);
        % Injection limits
        P <= repmat(P_max,1,S);
        P >= repmat(-P_max,1,S);
        % Charge limits
        Pc <= repmat(Pc_max,1,S);
        Pc >= repmat(-Pc_max,1,S);
        
        % Conservation of energy in batteries
        Pc == Estor-repmat(Estor0,1,S);
        
        % Battery capacity
        Estor <= repmat(E_max,1,S);        
        
        Pu <= reshape(Pl(:,t,:),[N S]);

        if (S >= 1)
            % If more than one scenario, decouple the load limit from used
            % power
            % The below constraints (in addition with above on state) enforce: Pu == min(min(v,Pl),x+Pg);
            Pu <= repmat(v,1,S);
            q1+q2 == 1;
            repmat(v,1,S) <= Pu + M.*(1-q1);
            repmat(v,1,S) <= reshape(Pl(:,t,:),[N S]) + M.*(1-q1);
            reshape(Pl(:,t,:),[N S]) <= Pu + M.*(1-q2);
            reshape(Pl(:,t,:),[N S]) <= repmat(v,1,S) + M.*(1-q2);
        end

        % SOS constraints
        sum(Estor,1) == (0:1/nXbins:1)*(r')*sum(E_max);
        r <= y;
        sum(r,2) == ones(S,1);
        sum(y,2) <= 2*ones(S,1);
        for i = 1:nXbins-1
            for j = i+2:nXbins+1
                y(:,i) + y(:,j) <= ones(S,1);
            end
        end
            
        
cvx_end

if ~strcmpi(cvx_status,'solved')
    error('CVX could not solve the optimization problem');
end

tSolve = toc(ts);

% With one scenario, we are setting consumption directly
if (S == 1)
    v = Pu(:,1,1);
end

% If load limit set close to max load over all scenarios, clear it
v(max(Pl(:,1,:),[],3) - v < 1e-4) = inf;
f = -cvx_optval/N;

% Comment in for Pinj
% Add nan values for users without batteries
%tmp = nan(N,1);
%tmp(battInd) = Pinj;
%Pinj = tmp;