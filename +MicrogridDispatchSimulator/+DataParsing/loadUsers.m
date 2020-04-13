function userTable = loadUsers(caseFolder)
%CREATE_UG returns a structure of microgrid users
%   Each user has a unique sequence of activities and thermal loads 
%   determined probabilistically from their user type description

% 2019, Claudio Vergara, Zola Electric.  

%% Read data   
userTable=readtable(fullfile(caseFolder,'users.csv'));

end

