function userTypes = loadUserTypes(dataFolder)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

import MicrogridDispatchSimulator.Models.ActivityType
import MicrogridDispatchSimulator.Models.LoadType
import MicrogridDispatchSimulator.Models.UserType

% Infer the type from the files
tmp = dir(dataFolder);
userTypeIds = {tmp(3:end).name}; % Returns 1 x numUserTypes cell array of type names
% Enforce unique

dayTypes = 1:3;

userTypes = containers.Map;

for k=1:length(userTypeIds)
    
    act_table = readtable(fullfile(dataFolder,[userTypeIds{k} '/activities.csv']));
    dayTables = cell(length(dayTypes),1);
    for d = dayTypes
        dayTables{d} = readtable(fullfile(dataFolder,sprintf('%s/day_%i.csv',userTypeIds{k},d)));
    end
    
    for a = 1:size(act_table,1)
        name = act_table{a,'Activity'};
        name = name{1}; % Unpack the cell array contents
        loadType = LoadType(act_table{a,'on_W'},name);
        startHourProbability = nan(24,length(dayTypes));
        for d = dayTypes
            dayTable = dayTables{d};
            startHourProbability(:,d) = dayTable{1:24,name};
        end
        activityTypes(a,1) = ActivityType(...
            loadType,act_table{a,'value_USD'},...
            act_table{a,'interrupt_USD'},act_table{a,'min_s'},...
            act_table{a,'max_s'},startHourProbability/100);
    end
    
    userTypes(userTypeIds{k}) = UserType(userTypeIds{k},activityTypes);
    
end
end

