function plotTimeSeries(simOut)
% Make a series of plots of time series outputs
%   simOut: A SimulationOutputs object storing time series data
%   microgrid: A Microgrid object with the microgrid parameters

% Variable shorthand
Estor = simOut.E;
Pu = simOut.Pl;
P = simOut.P;
Pg = simOut.Pg;
l = simOut.l;
Pinj = simOut.Pinj;
xi = simOut.xi;
E_max = simOut.Microgrid.Emax;
SOC = Estor./repmat(E_max,1,size(Estor,2));
time = simOut.timeSim/3600; % in hours
timeCntr = simOut.timeCntr/3600; % in hours

% System net load - solar
figure;
plot(time,sum(Pu-Pg)');
title('Total load minus total solar');
xlabel('Hour');
ylabel('W');

% Total load
figure;
plot(time,sum(Pu));
title('Total load');
xlabel('Hour');
ylabel('W');

% Individual load
figure;
plot(time,Pu');
title('Individual load');
xlabel('Hour');
ylabel('W');

% Total battery charge
figure;
yyaxis left
plot(time,sum(Estor));
ylabel('Wh');
yyaxis right
plot(time,xi);
title('Total battery energy')
legend('Estor','blackout');
xlabel('Hour');

% Total SoC
figure;
yyaxis left
plot(time,sum(Estor)/sum(E_max)*100);
ylabel('Percent');
yyaxis right
plot(time,xi);
title('Total battery SoC')
legend('SoC','blackout');
xlabel('Hour');

% Individual SoC
figure;
plot(time,SOC*100)
title('Individual SOC');
xlabel('Hour');
ylabel('Percent');

% Injection into grid
figure;
plot(timeCntr,Pinj)
title('Net injection setpoint Pinj')
xlabel('Hour');
ylabel('W');

figure;
plot(time,P)
title('Net injection actual P')
xlabel('Hour');
ylabel('W');

% Load limits
figure;
plot(timeCntr,l');
title('Individual load limits (control)');
xlabel('Hour');
ylabel('Wh per interval');