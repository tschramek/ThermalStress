function [tout,output] = weeklyMean(tin,input)
% [tout,output] = weeklyMean(tin,input)
%--------------------------------------------------------------------------
% Discription: computes weekly means using time and data vector inputs.
% Outputs both the time of the weekly means and the means themselves.
%
% Author:   T.Schramek
% Date:     2017.11.16
%--------------------------------------------------------------------------


% form weekly average
[yr,mo,dy,hr,mi,se] = datevec(tin);
iyr = unique(yr);
w = weeknum(tin);

% loop through each year
jj = 1;
for j = 1:length(iyr)
    % loops through each week 
    for j2 = 1:weeknum(datenum(iyr(j),12,31))
        k = find(tin>=datenum(iyr(j),1,1,0,0,0) & tin < datenum(iyr(j),12,31,11,59,59) & w == j2);
        output(jj) = nanmean(input(k));
        tout(jj) = nanmean(tin(k));
        if isnan(tout(jj))
            tout(jj) = datenum(iyr(j),1,1)+(j2*7)-3.5;
        end
        jj = jj+1;
    end
end
