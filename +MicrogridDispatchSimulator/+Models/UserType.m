classdef UserType
    properties (SetAccess = immutable)
        ActivityTypes
        Name
    end
    
    methods
        function self = UserType(name,activityTypes)
            self.ActivityTypes = activityTypes;
            self.Name = name;
        end
    end
end