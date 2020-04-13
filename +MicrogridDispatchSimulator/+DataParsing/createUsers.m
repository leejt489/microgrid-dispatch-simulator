function users = createUsers(params)
% Create users either by reading from table (if params.userTable is
% specified) or by number and type (if params.N and params.userTypes is
% defined). params.userTypes is assumed to be a dict (container.Map)
% Process for now will be to define different user types sequentially (see
% loop below).

% If params.userTable is not set, can specify params.NBatteries and
% params.NPV for number of battery and solar units, respectively, for each
% user as vectors of length N; along with params.BatteryUnitCapacity and
% params.PVUnitCapacity, which can either be scalars (as the same for all
% users) or vectors of length N.
% All of these must be set, otherwise they are all set to NaN

% If params.userTable is set, it can specify NBatteries and NPV as columns.
% The capacities BatteryUnitCapacity and PVUnitCapacity can either be
% specified as columns in the table, or in params as scalars or vectors of
% length N.

import MicrogridDispatchSimulator.Models.User

if (isfield(params,'N') && isfield(params,'userTypes'))
    N = params.N;
    userTypes = params.userTypes;
    if (isfield(params,'NBatteries') && isfield(params,'NPV') && ...
            isfield(params,'BatteryUnitCapacity') && ...
            isfield(params,'PVUnitCapacity'))
        NBatteries = params.NBatteries;
        NPV = params.NPV;
        BatteryUnitCapacity = params.BatteryUnitCapacity;
        PVUnitCapacity = params.PVUnitCapacity;
        
        % Check inputs for right format
        if (length(NBatteries) ~= N)
            error('params.NBatteries must be length N if it is defined');
        end
        if (length(NPV) ~= N)
            error('params.NPV must be length N if it is defined');
        end
        if (isscalar(BatteryUnitCapacity))
            BatteryUnitCapacity = BatteryUnitCapacity*ones(N,1);
        elseif (length(BatteryUnitCapacity) ~= N)
            error('params.BatteryUnitCapacity must be scalar or length N');
        end
        if (isscalar(PVUnitCapacity))
            PVUnitCapacity = PVUnitCapacity*ones(N,1);
        elseif (length(PVUnitCapacity) ~= N)
            error('params.PVUnitCapacity must be scalar or length N');
        end
    else
        NBatteries = nan(N,1);
        NPV = nan(N,1);
        BatteryUnitCapacity = nan(N,1);
        PVUnitCapacity = nan(N,1);
    end
    
    userTypeKeys = keys(userTypes);
    for n = 1:N
        type = userTypes(userTypeKeys{mod(n-1,length(userTypes))+1});
        users(n,1) = User(type,NBatteries(n),BatteryUnitCapacity(n),NPV(n),PVUnitCapacity(n));
    end
elseif (isfield(params,'userTable'))
    % Build from user table
    userTable = params.userTable;
    N = size(userTable,1);

    % Read NBatteries and NPV values from table
    if (ismember('NBatteries', userTable.Properties.VariableNames))
        NBatteries = userTable{:,'NBatteries'};
    else
        NBatteries = nan(N,1);
    end
    if (ismember('NPV', userTable.Properties.VariableNames))
        NPV = userTable{:,'PV'}; 
    else
        NPV = nan(N,1);
    end
    
    % Read capacity values either from table or params
    if (ismember('BatteryUnitCapacity', userTable.Properties.VariableNames))
        BatteryUnitCapacity = userTable{:,'BatteryUnitCapacity'}; 
    elseif (isfield(params,'BatteryUnitCapacity'))
        BatteryUnitCapacity = params.BatteryUnitCapacity;
        if (isscalar(BatteryUnitCapacity))
            BatteryUnitCapacity = BatteryUnitCapacity*ones(N,1);
        elseif (length(BatteryUnitCapacity) ~= N)
            error('params.BatteryUnitCapacity must be scalar or length N');
        end
    else
        BatteryUnitCapacity = nan(N,1);
    end
    if (ismember('PVUnitCapacity', userTable.Properties.VariableNames))
        PVUnitCapacity = userTable{:,'PVUnitCapacity'}; 
    elseif (isfield(params,'PVUnitCapacity'))
        PVUnitCapacity = params.PVUnitCapacity;
        if (isscalar(PVUnitCapacity))
            PVUnitCapacity = PVUnitCapacity*ones(N,1);
        elseif (length(PVUnitCapacity) ~= N)
            error('params.PVUnitCapacity must be scalar or length N');
        end
    else
        PVUnitCapacity = nan(N,1);
    end

    for n=1:N
        type = userTypes(char(userTable{n,'type'})); % Table converts char vector to string, so convert back
        users(n,1) = User(type,NBatteries(n),BatteryUnitCapacity(n),NPV(n),PVUnitCapacity(n));
    end
else
    error('Unrecognized arguments');
end

end

