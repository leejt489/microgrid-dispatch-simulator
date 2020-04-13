function [v,f,tSolve] = scenarioFormulation2Stage1LookGurobi(gamma,Pg,Pl,ps,Estor0,E_max,P_max,Pl_max,Pc_max,Pt_max,Bbus,Mbus)

params.OutputFlag = 0; % Suppress output
params.TimeLimit = 1800; % Set a 30 minute time limit on any run

ts = tic;

[N, T, S] = size(Pg);
NTS = N*T*S;
TS = T*S;
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
% Comment in for power flow constraints
% if (K > 0 && any(size(Pt_max) ~= [K 1]))
%     error('Ptmax must be Kx1 where K is number of rows of Mbus');
% end

% Index of users with batteries
%battInd = E_max > 0;
%Nbatt = sum(battInd);
%beta = Pc_max(battInd); % Approximate participation factor in power sharing

if (S > 1)
    M = max(Pl(:,1,:),[],3); %Static variable representing max possible for each load over all scenarios at time 1
end

% Pu = zeros(NTS,1);
% Pw = Pg(:);
% Estor = repmat(Estor0,TS,1);
% Pc = zeros(NTS,1);
% P = zeros(NTS,1);
% x = [Pu;P;Pc;Pw;Estor];

% Pu = x(1:NTS);
% P = x(NTS+1:2*NTS);
% Pc = x(2*NTS+1:3*NTS);
% Pw = x(3*NTS+1:4*NTS);
% Estor = x(4*NTS+1:5*NTS);
% v = x(5*NTS+1:5*NTS+N);
% q = x(5*NTS+N+1:5*NTS+N+NS);
indPu = 1:NTS;
indP = NTS+1:2*NTS;
indPc = 2*NTS+1:3*NTS;
indEstor = 3*NTS+1:4*NTS;
indv = 4*NTS+1:4*NTS+N;
indq = 4*NTS+N+1:4*NTS+N+NS;

% For each of the sub vectors y over NTS, y(n,t,s) = ...

% Length of x is 5*NTS+N+NS
if (S > 1)
    numX = 4*NTS+N+NS;
else
    numX = 4*NTS;
end


% Build objective: x'*Q*x+c'*x
Q = sparse(numX,numX);
Qu = kron(diag(sparse(ps)),kron(diag(sparse(gamma./Pl_max/2)),speye(T)));
Q(indPu,indPu) = Qu;

c = zeros(numX,1);
c(indPu) = -kron(ps',kron(gamma',ones(1,T)));

model.Q = Q;
model.obj = c;

% Build constraints
i = 1;
% Power balance at each node: P+Pw+Pu+Pc == Pg. Can drop Pw and write the
% constraint as 0 <= Pg-P-Pu-Pc <= Pg --> P+Pu+Pc <= Pg and -P-Pu-Pc <= 0
Ai = sparse(NTS,numX);
Ai(:,[indPu,indP,indPc]) = repmat(speye(NTS),1,3);
A(i:i+2*NTS-1,1:numX) = [Ai;-Ai];
b(i:i+2*NTS-1,1) = [Pg(:);zeros(NTS,1)];

i = i+2*NTS;
if (i-1 ~= size(A,1))
    error('Mistake in dimensions');
end

% Net injection sums to zero sum_n P == 0
Ai = sparse(TS,numX);
Ai(1:TS,indP) = kron(speye(TS),ones(1,N));
bi = zeros(TS,1);
A(i:i+2*TS-1,1:numX) = [Ai;-Ai];
b(i:i+2*TS-1,1) = [bi;-bi];

i = i+2*TS;
if (i-1 ~= size(A,1))
    error('Mistake in dimensions');
end

% Energy balance in battery -Pc(n,t,s) - Estor(n,t-1,s) + Estor(n,t,s) == 0;
% -Pc(n,t,s) + Estor(n,1,s) == Estor0(n);
Ai = sparse(NTS,numX);
Ai(1:NTS,indEstor) = kron(speye(S),kron(speye(T),speye(N))-[sparse(N,N*T);kron(speye(T-1),speye(N)) sparse(N*(T-1),N)]); % Creates entries for Estor
Ai(1:NTS,indPc) = -speye(N*T*S);
bi = kron(ones(S,1),[Estor0;zeros(N*(T-1),1)]);
A(i:i+2*NTS-1,1:numX) = [Ai;-Ai];
b(i:i+2*NTS-1,1) = [bi;-bi];

i = i+2*NTS;
if (i-1 ~= size(A,1))
    error('Mistake in dimensions');
end

if (S > 1)
    % Set Pu(n,1,s) == min(v(n),Pl(n,1,s) for all n,s using inequalities
    % and M(n) and binary q(n,s).
    % Constraints are:
    %   v(n) - Pu(n,1,s) + q(n,s)*M(n) <= M(n)
    %   v(n) + q(n,s)*M(n) <= Pl(n,1,s) + M(n)
    %   -Pu(n,1,s) - q(n,s)*M(n) <= -Pl(n,1,s)
    %   -v(n) - q(n,s)*M(n) <= -Pl(n,1,s)
    %   -v(n) + Pu(n,1,s) <= 0
    Ai = sparse(5*NS,numX);
    % Assign constraints rows depending on v
    Ai(1:2*NS,indv) = repmat(speye(N),2*S,1);
    Ai(3*NS+1:5*NS,indv) = repmat(-speye(N),2*S,1);
    % Assign constraints rows depending on Pu
    Ai(1:NS,indPu) = kron(speye(S),kron([-1 sparse(1, T-1)],speye(N)));
    Ai(2*NS+1:3*NS,indPu) = kron(speye(S),kron([-1 sparse(1, T-1)],speye(N)));
    Ai(4*NS+1:5*NS,indPu) = kron(speye(S),kron([1 sparse(1, T-1)],speye(N)));
    % Assign constraint rows depending on q
    Ai(1:2*NS,indq) = repmat(diag(repmat(M,S,1)),2,1);
    Ai(2*NS+1:4*NS,indq) = repmat(diag(repmat(-M,S,1)),2,1);
    % Assign right hand side
    t = Pl(:,1,:);
    bi = [repmat(M,S,1);repmat(M,S,1)+t(:);-t(:);-t(:);zeros(N*S,1)];
    
    A(i:i+5*NS-1,1:numX) = Ai;
    b(i:i+5*NS-1,1) = bi;
end

model.A = A;
model.rhs = b;
model.sense = '<';

% Set lower and upper bounds
lb = zeros(numX,1);
lb(indP) = repmat(-P_max,1,T,S); % -P_max <= P
lb(indPc) = repmat(-Pc_max,1,T,S); % -Pc_max <= Pc

ub = inf(numX,1);
ub(indPu) = Pl(:); % Pu <= Pl
ub(indP) = repmat(P_max,1,T,S); % P <= P_max
ub(indPc) = repmat(Pc_max,1,T,S); % Pc <= Pc_max
ub(indEstor) = repmat(E_max,1,T,S); % Estor <= E_max
if (S > 1)
    ub(indv) = M;
end

model.lb = lb;
model.ub = ub;

% Define variable types
if (S > 1)
    vtype = [repmat('C',5*NTS+N,1);repmat('B',NS,1)];
else
    vtype = 'C';
end
model.vtype = vtype;

result = gurobi(model,params);

tSolve = toc(ts);

% See Gurobi solution status codes: https://www.gurobi.com/documentation/9.0/refman/optimization_status_codes.html
if (~strcmp(result.status,'OPTIMAL'))
    if (strcmp(result.status,'SUBOPTIMAL'))
        warning('Suboptimal solution found');
    else
        error('Gurobi was unable to solve the problem');
    end
end

x = result.x;
% With one scenario, we are setting consumption directly
if (S == 1)
    v = full(x(1:N));
else
%     Pu = x(1:NTS);
%     P = x(NTS+1:2*NTS);
%     Pc = x(2*NTS+1:3*NTS);
%     Pw = x(3*NTS+1:4*NTS);
%     Estor = x(4*NTS+1:5*NTS);
%     v = x(5*NTS+1:5*NTS+N);
%     q = x(5*NTS+N+1:5*NTS+N+NS);
    
    v = x(5*NTS+1:5*NTS+N);
end

% If load limit set close to max load over all scenarios, clear it
v(max(Pl(:,1,:),[],3) - v < 1e-4) = inf;

% Record objective value as a benefit (positive) per user
f = -result.objval/N;