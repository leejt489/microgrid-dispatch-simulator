function params = readKeyValue(path)
%READKEYVALUE Reads a key value table into a struct
params = struct;

tbl = readtable(path);
% Read key/value pairs into a struct
for i = 1:size(tbl,1)
    params.(string(tbl{i,1})) = tbl{i,2};
end