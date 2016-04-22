function [output] = ave_convol (data,leaflet)
% Calculate the average of the data points before and after to acquire
% trend in the data

% For this particular project, a leaflet of 75 will be taking roughly 10ms
% of the data

% Input can be a dataset of column vector (n by m)
datalength = size(data,2);
numset = size(data,1);
output = zeros(numset,datalength);

for i = 1:numset
    for j = 1:datalength
        % Check for the end-cases
        if j < leaflet+1
            output(i,j) = mean(data(i,1:j+leaflet));
        elseif j > datalength - leaflet
            output(i,j) = mean(data(i,j-leaflet:end));
        % Find an average of points
        else
            output(i,j) = mean(data(i,j-leaflet:j+leaflet));
        end
    end
end




