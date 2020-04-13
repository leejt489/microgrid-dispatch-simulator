function [v,f,tSolve] = scenarioFormulationBackRecGurobi(gamma,Pg,Pl,ps,Estor0,E_max,P_max,Pl_max,Pc_max)%,Pt_max,Bbus,Mbus)

params.OutputFlag = 0; % Suppress output
params.TimeLimit = 300; % Set a 5 minute time limit on any run

% State space discretization. Consider the state to be one aggregated
% battery
nXBins = 10; % Number of samples per state.

ts = tic;

[N, T, S] = size(Pg);
NS = N*S;
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

% Set future cost variable, initialize to cost of zero
V1 = zeros(nXBins+1,1);
V0 = V1;

t = T;
    
% Pu = x(1:NS);
% PcTot = x(NS+1:NS+S);
% PwTot = x(NS+S+1:NS+2*S);
% EstorTot = x(NS+2*S+1:NS+3*S); % Next period total energy
% v = x(NS+3*S+1:NS+3*S+N); % Load subject to load limit
% q = x(NS+3*S+N+1:2*NS+3*S+N); % Binary variable whether or not load
% % limit is binding
% r = x(2*NS+3*S+N+1:2*NS+3*S+N+(nXBins+1)*S) % Weight for piecewise linear
indPu = 1:NS;
indPcTot = indPu(end)+(1:S);
indEstorTot = indPcTot(end)+(1:S);
indv = indEstorTot(end)+(1:N);
indq = indv(end)+(1:NS);
indr = indq(end)+(1:(nXBins+1)*S);

if (S > 1)
    numX = indr(end);
    numCon = 8*S + 5*NS;
else
    error('need to check this');
end

% Quadratic part of cost
Q = sparse(numX,numX);
Qu = kron(diag(sparse(ps)),diag(sparse(gamma./Pl_max/2)));
Q(indPu,indPu) = Qu;

% Linear part of cost
c = zeros(numX,1);
c(indPu) = -kron(ps',gamma'); % These elements will be the same every time step, but also need to add linear interpolation of future

% Build constraints
A = sparse(numCon,numX);
b = nan(numCon,1);

% Power balance on the system: PwTot+sum_n(Pu)+PcTot == sum_n(Pg). Rewrite
% dropping PwTot as 0 <= sum_n(Pg) - sum_n(Pu) - PcTot <= sum_n(Pg)
% which goes to sum_n(Pu) + PcTot <= sum_n(Pg) and -sum_n(Pu) - PcTot <= 0
Ai = sparse(S,numX);
Ai(:,indPu) = full(kron(speye(S),ones(1,N)));
Ai(:,indPcTot) = speye(S);

A(1:2*S,1:numX) = [Ai;-Ai];
b(S+1:2*S) = 0;
% upper chunk of b corresponding to Pg will be set for each time

% Charge balance for battery EstorTot - PcTot == E0
Ai = sparse(S,numX);
Ai(:,indPcTot) = -speye(S); % -PcTot
Ai(:,indEstorTot) = speye(S); % EstorTot

A(2*S+1:4*S,1:numX) = [Ai;-Ai];
% b will be set by iterating from E0

if (S > 1)
    % Set Pu(n,1,s) == min(v(n),Pl(n,1,s) for all n,s using inequalities
    % and M(n) and binary q(n,s).
    % Constraints are:
    %   v(n) - Pu(n,s) + q(n,s)*M(n) <= M(n);   row indices 1:NS
    %   v(n) + q(n,s)*M(n) <= Pl(n,s) + M(n);   row indices NS+1:2*NS
    %   -Pu(n,s) - q(n,s)*M(n) <= -Pl(n,s);     row indices 2*NS+1:3*NS
    %   -v(n) - q(n,s)*M(n) <= -Pl(n,s);        row indices 3*NS+1:4*NS
    %   -v(n) + Pu(n,s) <= 0;                   row indices 4*NS+1:5*NS
    Ai = sparse(5*NS,numX);
    % Assign constraints rows depending on v
    Ai(1:2*NS,indv) = repmat(speye(N),2*S,1);
    Ai(3*NS+1:5*NS,indv) = repmat(-speye(N),2*S,1);
    % Assign constraints rows depending on Pu
    Ai(1:NS,indPu) = -speye(NS);
    Ai(2*NS+1:3*NS,indPu) = -speye(NS);
    Ai(4*NS+1:5*NS,indPu) = speye(NS);
    % Assign constraint rows depending on q (this has to be updated for
    % each time as M changes)
    % Ai(1:2*NS,NS+3*S+N+1:2*NS+3*S+N) = repmat(diag(repmat(M,S,1)),2,1);
    % Ai(2*NS+1:4*NS,NS+3*S+N+1:2*NS+3*S+N) = repmat(diag(repmat(-M,S,1)),2,1);
   
    A(4*S+1:4*S+5*NS,1:numX) = Ai;
    
    b(4*S+4*NS+1:4*S+5*NS) = zeros(NS,1); % Other rhs will vary
end

% Piecewise linear constraint ((0:1/nXbins:1)*ETotMax*r)' - EstorTot == 0 where r
% is nXbins+1 x S
Ai = sparse(S,numX);
Ai(:,indEstorTot) = -speye(S);
Ai(:,indr) = kron(speye(S),(0:1/nXBins:1)*sum(E_max));

A(4*S+5*NS+1:6*S+5*NS,1:numX) = [Ai;-Ai];
b(4*S+5*NS+1:6*S+5*NS) = zeros(2*S,1);

% Sum of piecewise linear weights r = 1 for each s.
Ai = sparse(S,numX);
Ai(1:S,indr) = kron(speye(S),ones(1,nXBins+1));

A(6*S+5*NS+1:8*S+5*NS,1:numX) = [Ai;-Ai];
b(6*S+5*NS+1:8*S+5*NS) = [ones(S,1);-ones(S,1)];

% Add SOS2 constraint for linear interpolation
for s = 1:S
    model.sos(s).type = 2;
    model.sos(s).index = indr((s-1)*(nXBins+1)+1:s*(nXBins+1));
end

% Set lower and upper bounds
lb = zeros(numX,1);
lb(indPcTot) = repmat(-sum(Pc_max),S,1); % PcTot >= -sum(Pc_max)

ub = inf(numX,1);
ub(indPcTot) = repmat(sum(Pc_max),S,1); % PcTot <= sum(Pc_max)
ub(indEstorTot) = repmat(sum(E_max),S,1); % EstorTot <= sum(E_max)


model.Q = Q;
model.obj = c;

model.A = A;
model.rhs = b;
model.sense = '<';

model.lb = lb;
model.ub = ub;

% Define variable types
if (S > 1)
    vtype = repmat('C',numX,1);
    vtype(indq) = 'B';
else
    vtype = 'C';
end
model.vtype = vtype;

while (t > 0)
    
    % Piecewise linear interpolation of future value given points in V1,
    % for each scenario
    model.obj(indr) = kron(ps',V1');
    
    if (S > 1)
        M = max(Pl(:,t,:),[],3); %Static variable representing max possible for each load over all scenarios at the time t
    end
    
    % Update constraint rows depending on q for new values of M
    model.A(4*S+1:4*S+2*NS,indq) = repmat(diag(repmat(M,S,1)),2,1);
    model.A(4*S+2*NS+1:4*S+4*NS,indq) = repmat(diag(repmat(-M,S,1)),2,1);
    
    Plt = Pl(:,t,:);
    model.rhs(4*S+1:4*S+4*NS) = [repmat(M,S,1);repmat(M,S,1)+Plt(:);-Plt(:);-Plt(:)];
    
    % Update additional time dependent RHS
    Pgt = sum(Pg(:,t,:));
    model.rhs(1:S) = Pgt(:);
    
    
    model.ub(indPu) = Pl(:,t,:); % Pu <= Pl
    
    if (t > 1)
        % Sample the cost at discrete values
        for xBinInd = 0:nXBins % Iterate over possible values at current state

            % Constraint for current inititial state of charge
            EstorTot0 = sum(E_max)*xBinInd/nXBins;
            model.rhs(2*S+1:4*S) = EstorTot0*[ones(S,1);-ones(S,1)];

            result = gurobi(model,params);

            % See Gurobi solution status codes: https://www.gurobi.com/documentation/9.0/refman/optimization_status_codes.html
            if (~strcmp(result.status,'OPTIMAL'))
                if (strcmp(result.status,'SUBOPTIMAL'))
                    warning('Suboptimal solution found');
                elseif (strcmp(result.status,'INF_OR_UNBD') || strcmp(result.status,'INFEASIBLE'))
                    % Try removing presolve and re-run
                    paramsSafe = params;
                    paramsSafe.DualReductions = 0;
                    paramsSafe.Presolve = 0;
                    result2 = gurobi(model,paramsSafe);
                    if (strcmp(result2.status,'OPTIMAL'))
                        result = result2;
                    else
                        save('gurobidebug');
                        error('Gurobi was unable to solve the problem, result code: %s',result.status);
                    end
                else
                    save('gurobidebug');
                    error('Gurobi was unable to solve the problem, result code: %s',result.status);
                end
            end

            % Store the optimal value for this state
            V0(xBinInd+1) = result.objval;

            % We know that the V(x) is monotonically decreasing in x and that
            % it is strictly monotonic up to a point and then constant. So
            % check if it is constant to avoid unnecessarily computing the
            % rest.
            if (xBinInd > 0 && abs(1-V0(xBinInd)/(V0(xBinInd+1)+0.0001)) < 1e-3)
                V0(xBinInd+1:end) = V0(xBinInd);
                break;
            end
        end
        
        % Set the next step cost samples as the current step cost samples
        V1 = V0;
    else
        % We are down to the first time step, so we only care about the
        % initial state
        EstorTot0 = sum(Estor0);
        
        model.rhs(2*S+1:4*S) = EstorTot0*[ones(S,1);-ones(S,1)];

        result = gurobi(model,params);

        % See Gurobi solution status codes: https://www.gurobi.com/documentation/9.0/refman/optimization_status_codes.html
        if (~strcmp(result.status,'OPTIMAL'))
            if (strcmp(result.status,'SUBOPTIMAL'))
                warning('Suboptimal solution found');
            elseif (strcmp(result.status,'INF_OR_UNBD') || strcmp(result.status,'INFEASIBLE'))
                % Try removing presolve and re-run
                paramsSafe = params;
                paramsSafe.DualReductions = 0;
                paramsSafe.Presolve = 0;
                result2 = gurobi(model,paramsSafe);
                if (strcmp(result2.status,'OPTIMAL'))
                    result = result2;
                else
                    save('gurobidebug');
                    error('Gurobi was unable to solve the problem, result code: %s',result.status);
                end
            else
                save('gurobidebug');
                error('Gurobi was unable to solve the problem, result code: %s',result.status);
            end
        end
    end
    % Move back a time step
    t = t-1;
end

tSolve = toc(ts);

x = result.x;

% With one scenario, we are setting consumption directly
if (S == 1)
    %v = Pu(:,1,1);
else
    v = x(indv);
end

% If load limit set close to max load over all scenarios, clear it
v(max(Pl(:,1,:),[],3) - v < 1e-4) = inf;
f = -result.objval/N;
