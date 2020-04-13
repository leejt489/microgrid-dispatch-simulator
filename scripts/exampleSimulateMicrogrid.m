%% Simulate Microgrid
% This script creates and simulates a microgrid with distributed solar and
% battery storage, and load dispatch
%% Import dependencies
import MicrogridDispatchSimulator.DataParsing.*
import MicrogridDispatchSimulator.Models.Microgrid
import MicrogridDispatchSimulator.Simulation.simRHC
import MicrogridDispatchSimulator.Visualization.plotTimeSeries
%% Define parameters

% Timing parameters
deltaT_sim = 60; % 1 minutes (in seconds) simulation interval
deltaT_cntr = 3600; % 1 hour RHC interval;
T = 7*24*3600; % 1 week (in seconds)
dayOfYear = 90; % Start day of year

% User parameters
N = 5; % Number of users
nBatteries = zeros(N,1);
nBatteries([1 2 5]) = [2 2 4];
nPV([2 5]) = [3 4];
userMaxLoad = 10000; % Max load (W) for each user;

% Set controller
controllerName = 'Reactive'; % Reactive control; TODO: implement stand-alone forecast model for predictive controllers

%% Load and map remaining parameters
% The rest of the parameters are either loaded from files, take default
% values, or are computed from mappings.

% Read DER params
derParams = readKeyValue('data/der.csv');

% Params struct for users
userParams = struct;
userParams.N = N;
userParams.userTypes = loadUserTypes('data/user_types/');
userParams.NBatteries = nBatteries;
userParams.NPV = nPV;
userParams.BatteryUnitCapacity = derParams.batteryEnergyCap*1000;
userParams.PVUnitCapacity = derParams.solarInverterCap*1000;

% Set parameters for the microgrid
microgridParams = struct;
microgridParams.UserMaxLoad = userMaxLoad; % 10 kW max load per user
microgridParams.ERestart = 0.1; % Restart from blackout at 10% SoC
microgridParams.BatteryChargeRate = derParams.batteryInverterCap/derParams.batteryEnergyCap;
microgridParams.Beta = derParams.beta;

% Set up time param struct
timeParams = struct;
timeParams.T = T/deltaT_sim;
timeParams.deltaT_sim = deltaT_sim;
timeParams.deltaT_cntr = deltaT_cntr;

% Set up control param struct
controlParams = struct;
controlParams.loadControllerName = controllerName;
controlParams.K = 0.5;

% Load solar data
disturbances = struct;
disturbances.irradiance = loadSolarData('data/solar/solar_ghi_data', dayOfYear, T, deltaT_sim, 1);

%% Create users and microgrid

% Create users with DERs and initialize activities
users = createUsers(userParams);
for i = 1:length(users)
    users(i).createActivities(T,deltaT_sim);
end
microgridParams.BusMaxInjection = inf;

% Create and initialize the microgrid
microgrid = Microgrid(microgridParams);
microgrid.Users = users; % Connect the users
microgrid.initialize(); % Initialize

%% Run the simulation

results = simRHC(microgrid,timeParams,controlParams,disturbances,[],true);

%% Plot outputs
plotTimeSeries(results);
