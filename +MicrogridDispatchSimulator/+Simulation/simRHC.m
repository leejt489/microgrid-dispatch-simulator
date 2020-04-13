function outputs = simRHC(microgrid,timeParams,controlParams,disturbances,forecasts,verbose)

import MicrogridDispatchSimulator.Models.MicrogridController
import MicrogridDispatchSimulator.Models.SimulationOutputs

users = microgrid.Users; % Shorthand reference to array of users
N = length(users); % Shorthand number of users

T = timeParams.T; % Number of time steps in simulation
deltaT_sim = timeParams.deltaT_sim; % Length of simulation time step in seconds
deltaT_cntr = timeParams.deltaT_cntr; % Length of time between control updates in seconds


% Set up microgrid controller
controlParams.deltaT_cntr = deltaT_cntr;
controlParams.gamma = ones(N,1);
if (isfield(timeParams,'deltaT_hor'))
    controlParams.deltaT_hor = timeParams.deltaT_hor;
end

% Initialize microgrid controller
microgridController = MicrogridController(microgrid,controlParams.loadControllerName,controlParams);


% Initialize data store of outputs
outputs = SimulationOutputs(microgrid,T,deltaT_cntr,deltaT_sim);

% Setup input struct from microgrid controller
uMGC = struct;
if (microgridController.HasMemory)
    % Initialize field to attach the data store of outputs
    uMGC.y = outputs;
end
if (microgridController.RequiresSolarForecast || microgridController.RequiresLoadForecast)
    % Initialize field to attach forecast to microgrid controller inputs
    uMGC.W = struct;
    uMGC.W.ps = forecasts.ps;
    T_hor = timeParams.T_hor;
    deltaT_hor = timeParams.deltaT_hor;
end

printProgress = nargin >= 6 && verbose; % Set to true to print output, good for long simulations
if (printProgress)
    printInterval = 2; % Seconds between print out
    tComp = tic; % Timer for print out
    s = ''; % String for printing
end

% Start simulation loop
i_cntr = 0;
for i_sim = 0:T-1
    t = i_sim*deltaT_sim; % Absolute time in seconds
    
    if (printProgress && (toc(tComp) >= printInterval || isempty(s)))
        % Print progress
        if (~isempty(s))
            fprintf(repmat('\b',1,length(s))); % Delete old line
        end
        s = sprintf('Time: %i of %i; Percent: %2.1f\r',i_sim,T,i_sim/T*100);
        fprintf(s);
        tComp = tic;
    end

    if (mod(t,deltaT_cntr) - deltaT_sim < 0)
        % Take control action
        i_cntr = i_cntr+1;
        
        % TODO: move this to separate forecast object
        if (microgridController.RequiresSolarForecast)
            % Construct forecast object
            % Resample time series to appropriate resolution
            uMGC.W.Pg = forecasts.Pg(:,i_sim/(deltaT_hor/deltaT_sim)+1:min(i_sim/(deltaT_hor/deltaT_sim)+T_hor,end),:);
        end
        if (microgridController.RequiresLoadForecast)
            uMGC.W.Pl = forecasts.Pl(:,i_sim/(deltaT_hor/deltaT_sim)+1:min(i_sim/(deltaT_hor/deltaT_sim)+T_hor,end),:);
        end
        
        % Microgrid controller computes load limits and injection points,
        % sends them to users, microgrid, and meters.
        yMGC = microgridController.update(t,uMGC);
        outputs.saveMGC(t,yMGC);
    end
    
    % Update user and save output
    for n = 1:N
        users(n).update(t,deltaT_sim);
    end
    % We update the user first to present the microgrid with an updated
    % load. This means that if the microgrid registers a service
    % interruption it should be interpreted as happening within the
    % interval. This means that user activities should not move to
    % completed unless they've confirmed there was no interruption, and
    % that the load should reflect what was running during this period.

    % Update microgrid and save output
    uMG.irradiance = disturbances.irradiance(i_sim+1);
    yMG = microgrid.update(deltaT_sim,uMG); % Update microgrid one state
    outputs.saveMG(t,yMG);
    
    % Update meters
    for n = 1:N
        users(n).Meter.update(t,deltaT_sim);
    end

    
end

end

