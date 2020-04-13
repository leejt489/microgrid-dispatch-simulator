function [v,f,tSolve] = scenarioFormulation2Stage1Look(gamma,Pg,Pl,ps,Estor0,E_max,P_max,Pl_max,Pc_max,Pt_max,Bbus,Mbus)

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
% Comment in for power flow constraints
% if (K > 0 && any(size(Pt_max) ~= [K 1]))
%     error('Ptmax must be Kx1 where K is number of rows of Mbus');
% end

% Index of users with batteries
%battInd = E_max > 0;
%Nbatt = sum(battInd);
%beta = Pc_max(battInd); % Approximate participation factor in power sharing

if (S >= 1)
    M = repmat(max(Pl(:,1,:),[],3),1,1,S); %Static variable representing max possible for each load over all scenarios at time 1
end

ts = tic;
cvx_begin
    cvx_quiet true
    % Comment in to include Pinj
    %variable Pinj(Nbatt,1) % Power injected
    
    % System variables
    variable Estor(N,T,S) nonnegative %x0 is fixed so x(N,1,S) is x1; i.e. x is shifted backwards in time by 1 relative to other variables
    variable P(N,T,S) % Realized net injection
    variable Pu(N,T,S) nonnegative % Power used in each scenario
    variable Pc(N,T,S) % Charge power into batteries
    variable Pw(N,T,S) nonnegative % Power wasted
    %variable d(N,T,S)
    % Auxiliary variables
    if (S >= 1)
        variable v(N,1) nonnegative %'v' is the load limit; 'l' in the paper, but 'l' looks too much like the number 1
        variable q1(N,1,S) binary %q1(n,s) = 1 => v(n) <= Pl(n,s)
        %variable q2(N,1,S) binary %Pl(T,S) <= min(v(T),x(T+S)+Pg(T,s)
    end
    
    % Quadratic cost summed across scenarios.
    minimize(...
        Pu(:)'*kron(diag(sparse(ps)),kron(diag(sparse(gamma./Pl_max/2)),sparse(eye(T))))*Pu(:)...
        -kron(ps',kron(gamma',ones(1,T)))*Pu(:)...
        ...%+ d(:)'*kron(diag(ps),kron(eye(T),Bbus'*Bbus)/2*0.1)*d(:)...
        ...%+ Pw(:)'*kron(diag(ps),eye(N*T)/2*0.01)*Pw(:)...
        ...% TODO: add final cost
    )
    
    subject to
        % Power curtailment must be less than generation
        Pw <= Pg;
%         for s = 1:S
%             for t = 1:T
%                 if (t == 1)
%                     x(:,t,s) == x0+Pg(:,t,s)-Pu(:,t,s)-Pw(:,t,s)-Bbus*d(:,t,s);
%                 else
%                     x(:,t,s) == x(:,t-1,s)+Pg(:,t,s)-Pu(:,t,s)-Pw(:,t,s)-Bbus*d(:,t,s);
%                 end
%                 x(:,t,s) <= x_max;
%                 Mbus*d(:,t,s) <= Pt_max;
%                 -Mbus*d(:,t,s) <= Pt_max;
%             end
%         end
 
        %d(1,:,:) == 0
        % Relate power injection to setpoint by the beta characteristic
        % (participation factor)
        %repmat(Pinj,1,S)-reshape(P(battInd,1,:),[Nbatt, S]) == (beta/sum(beta)*ones(1,Nbatt))*(repmat(Pinj,1,S)-reshape(P(battInd,1,:),[Nbatt, S]));
        % Ensure injection setpoint satisfied for all scenarios for users
        % that have a battery to pick up the slack
        
        % Commentin out net injection
        %repmat(Pinj,1,1,S) == P(battInd,1,:);
        
        % Flow balance at each node
        P == Pg-Pw-Pu-Pc;
        % All constraints sum to zero
        %sum(P) == zeros(1,T,S);
        reshape(sum(P,1),[T, S]) == sparse(zeros(T,S));
        % Injection limits
        P <= repmat(P_max,1,T,S);
        P >= repmat(-P_max,1,T,S);
        % Charge limits
        Pc <= repmat(Pc_max,1,T,S);
        Pc >= repmat(-Pc_max,1,T,S);
        
        % Conservation of energy in batteries
        Pc(:,1,:) == Estor(:,1,:)-repmat(Estor0,1,1,S);
        if (T > 1)
            for t = 2:T
                Pc(:,t,:) == Estor(:,t,:)-Estor(:,t-1,:);
            end
        end
        
        % Battery capacity
        Estor <= repmat(E_max,1,T,S);
        
        
        Pu <= Pl;
        
        if (S >= 1)
            % If more than one scenario, decouple the load limit from used
            % power
            
            % The below constraints (in addition with above on state) enforce: Pu == min(min(v,Pl),x+Pg);
            Pu(:,1,:) <= repmat(v,1,1,S);
            %q1+q2 == 1;
            % Ensure Pu = v and v <= min(Pl,x+Pg) when q1 true
            repmat(v,1,1,S) <= Pu(:,1,:) + M.*(1-q1);
            repmat(v,1,1,S) <= Pl(:,1,:) + M.*(1-q1);
            %repmat(v,1,1,S) <= min(Pl(:,1,:),x0+Pg(:,1,:)) + M.*(1-q1);
            % Ensure Pu = min(Pl,x0+Pg) and min(Pl,x0+Pg) <= v when q2 true
            %Pl(:,1,:) <= Pu(:,1,:) + M.*(1-q2);
            %Pl(:,1,:) <= repmat(v,1,1,S) + M.*(1-q2);
            Pl(:,1,:) <= Pu(:,1,:) + M.*q1;
            Pl(:,1,:) <= repmat(v,1,1,S) + M.*q1;
            %min(Pl(:,1,:),x0+Pg(:,1,:)) <= Pu(:,1,:) + M.*(1-q2);
            %min(Pl(:,1,:),x0+Pg(:,1,:)) <= repmat(v,1,1,S) + M.*(1-q2);
        end
            
        
cvx_end

if ~strcmpi(cvx_status,'solved')
    if strcmpi(cvx_status,'suboptimal')
        warning('Suboptimal solution found');
    else
        error('CVX could not solve the optimization problem');
    end
end
tSolve = toc(ts);

% With one scenario, we are setting consumption directly
if (S == 1)
    v = full(Pu(:,1));
end

% If load limit set close to max load over all scenarios, clear it
v(max(Pl(:,1,:),[],3) - v < 1e-4) = inf;
f = -cvx_optval/N; % Returned expected objective value as maximizing benefit per user

% Comment in to set injection setpoints
% Add nan values for users without batteries
%tmp = nan(N,1);
%tmp(battInd) = Pinj;
%Pinj = tmp;