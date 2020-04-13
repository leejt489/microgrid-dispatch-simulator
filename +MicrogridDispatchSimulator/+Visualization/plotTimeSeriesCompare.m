function plotTimeSeriesCompare(outputs1,outputs2,label1,label2)
% Compare time series from two results. Currently only looks at energy
% stored, used power, and blackouts, but it can be expanded. The plots show
% a series of subplots for each result and then the difference between
% them.

% Arguments:
%   outputs1: SimulationOutputs object from one simulation
%   outputs2: SimulationOutputs object from second simulation. Should have
%       same time resolution as the first.
%   label1 and label2: optional strings to label 

if (nargin < 3)
    label1 = 'Results 1';
end

if (nargin < 4)
    label2 = 'Results 2';
end

diff_str = sprintf('%s - %s',label1,label2);

Estor = outputs1.E;
Pu = outputs1.Pl;
xi = outputs1.xi;

Estor_j = outputs2.E;
Pu_j = outputs2.Pl;
xj = outputs2.xi;

time = outputs1.timeSim/3600; % in hours; time should be the same for each result

%% Battery charge
figure;
subplot(311);
yyaxis left
plot(time,sum(Estor));
yyaxis right
plot(time,xi);
title(['Total battery Wh: ', label1] )
legend('Estor','blackout');
subplot(312);
yyaxis left
plot(time,sum(Estor_j));
yyaxis right
plot(time,xj);
title(['Total battery Wh: ', label2] )
legend('Estor','blackout');
subplot(313);
yyaxis left
plot(time,sum(Estor)-sum(Estor_j));
yyaxis right
plot(time,xi-xj);
title(['Total battery kWh: ', diff_str] )
legend('Estor','blackout');
xlabel('Hour');

%% Individual load
figure;
subplot(311);
plot(time,Pu');
title(['Load: ', label1]);
ylabel('W');
subplot(312);
plot(time,Pu_j');
title(['Load: ', label2]);
ylabel('W');
subplot(313);
plot(time,Pu'-Pu_j');
title(['Load: ', diff_str]);
ylabel('W');
xlabel('Hour');