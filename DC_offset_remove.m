function [output] = DC_offset_remove(data)
% Removes DC offsets

% Input is a dataset of column vector (n by m)
DC = mean(data);
output = data - DC;


end

