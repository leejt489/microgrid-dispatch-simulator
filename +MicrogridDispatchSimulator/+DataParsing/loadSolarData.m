function [irradiance] = loadSolarData(dataFilepath, startDayOfYear, horizon, desiredResolution, numYears, nativeResolution)
%SOLAR_FORECAST Generate a set of solar realizations from data
%   Example solar forecast generation from CAMS Radiation Service:
%   Data is manipulated in python with minimal dependencies and output in a
%   .mat file at 1 minute (native) resolution for the current data
%
%   dataFilepath: str, the location of the .mat file to load
%   startDayOfYear: int, day between 1 and 365 (minus time horizon)
%   horizon: int, time window of data requested in seconds
%   desiredResolution: int, seconds per sample in output data
%   numYears: int, number of samples of different years
%   nativeResolution: int, seconds per sample in stored data; set to 60 for 1 minute data, 1 for 1 second data, etc.

import MicrogridDispatchSimulator.Utilities.resampleBasic

solarGHI = load(dataFilepath); % is in cumulative Wh/m^2 per minute (W/m^2) JTL: I think it
solarGHI = solarGHI.solar_ghi*60/1000; % extract the matrix from the struct; rows are minutes of the year and columns are different years

if (nargin < 5) || isnan(numYears) || isempty(numYears)
    numYears = size(solarGHI,2);
end
if (nargin < 6)
    nativeResolution = 60;
end
if (nargin < 4)
    desiredResolution = nativeResolution;
end
startMinute = (startDayOfYear-1)*60*24+1; % get the starting time in minutes. Minus 1 is for start of that day
irradiance = solarGHI(startMinute:startMinute+horizon/nativeResolution-1,end-numYears+1:end); % Get most recent years
% Takes the 1 minute data and transforms it to the simulation frequency
if (nativeResolution ~= desiredResolution)
    irradiance = resampleBasic(irradiance, nativeResolution, desiredResolution); % Average irradiance over each time period in kW/m^2
end