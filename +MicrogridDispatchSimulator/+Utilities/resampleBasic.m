function data = resampleBasic(data, input_ts, output_ts,interpOrder)

if nargin < 4
    interpOrder = 1;
end

resampleRate = output_ts / input_ts;
if input_ts < output_ts
    % Downsample / average
    data2 = nan(size(data)*[1/resampleRate 0;0 1]);
    for i = 1:size(data2,1)
        data2(i,:) = mean(data((i-1)*resampleRate+1:i*resampleRate,:));
    end
    data = data2;
elseif input_ts > output_ts
    % Upsample / interpolate
    if (round(1/resampleRate) ~= 1/resampleRate)
        error('Solar resampling must be done by integer multiple');
    end
    data2 = nan(size(data)*[1/resampleRate 0;0 1]);
    filler_points = linspace(0,1,1/resampleRate);
    ones_vector = ones(length(filler_points),1);
    % Linear interpoloation
    for i = 1:size(data,1)-1
        if (interpOrder == 1)
            data2((i-1)/resampleRate+1:i/resampleRate,:) = ones_vector*data(i,:)+filler_points'*(data(i+1,:)-data(i,:));
        elseif (interpOrder == 0)
            data2((i-1)/resampleRate+1:i/resampleRate,:) = ones_vector*data(i,:);
        else
            error('Interpolation order not recognized');
        end
    end
    % For last block, no endpoint given, so assume constant
    data2(i/resampleRate+1:(i+1)/resampleRate,:) = ones_vector*data(i+1,:);
    data = data2;
   
end