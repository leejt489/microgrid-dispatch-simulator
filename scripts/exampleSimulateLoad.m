import MicrogridDispatchSimulator.DataParsing.loadUserTypes
import MicrogridDispatchSimulator.DataParsing.createUsers
import MicrogridDispatchSimulator.Simulation.simLoad
import MicrogridDispatchSimulator.Models.MeterRelay

N = 3; % Number of users
deltaT = 60; % 1 minutes (in seconds)
T = 7*24*3600; % 1 week (in seconds)


% Params struct for users
p = struct;
p.N = N;
p.userTypes = loadUserTypes('data/user_types/');

% Create users
users = createUsers(p);
% Create a dummy bus to connect loads to (they must have a source)
m = MeterRelay();

results = struct; % Struct array, will be length N
for n = 1:N
    user = users(n);
    user.Meter = m; % Connect to the bus
    user.createActivities(T,deltaT); % Need to explicitly create a set of activities over the time window at a certain resolution
    user.initializeState(); % Need to initialize the state
    [results(n).Pl,results(n).Pa,results(n).loadNames] = simLoad(0,struct,T,deltaT,user);
end

% Plot load profiles
f = figure;
time = (0:deltaT:T-deltaT)/3600;
for n = 1:N
    figure(f);
    r = results(n);
    plot(time,r.Pl);
    hold on
    figure;
    plot(time,r.Pa);
    legend(r.loadNames{:});
    xlabel('Time (h)');ylabel('W');
    title(sprintf('User %i Power Consuption By Individual Load',n))
end

figure(f);
xlabel('Time (h)');ylabel('W');
title('Total Power Consumption by User')
legend(arrayfun(@(n) num2str(n),1:N,'UniformOutput',false))
    