function [Pl,Pa,loadNames] = simLoad(t,u,T,deltaT,user)
% This function simulates the load for a user from start time 't' with time
% step 'deltaT'. 'u' is a struct of inputs, which at this time must have a
% field 'ambientTemp' which is a vector of ambient temperatures assumed to
% be sampled at 'deltaT'. 
% If Pa is requested, it will be a T/deltaT x Number of Loads array of
% power by loads

%tAmbient = u.ambientTemp;

absTime=t:deltaT:t+T-deltaT; % in seconds

if (nargout > 1)
    % Compute power by load
    Pa = nan(length(absTime),length(user.Loads));
    for i = 1:length(absTime)
        Pa(i,:) = user.Pa';
        user.update(absTime(i),deltaT);
    end
    Pl = sum(Pa,2); % Total power is sum across loads
    if (nargout > 2)
        loadNames = user.LoadNames; % Return cell array of strings for load name
    end
else
    % Compute total power only
    Pl = nan(length(absTime),1);
    for i = 1:length(absTime)
        Pl(i) = user.Pl;
        user.update(absTime(i),deltaT);
    end
end

end