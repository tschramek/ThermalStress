function [DHW, tDHW, MMM, MMC, degT] = computeDHW(t,temp,varargin)
% [DHW, MMM, MM, degT] = computeDHW(t,temp,MMMin,sstStd)
%--------------------------------------------------------------------------
% Description: Compute DHW for a time series where the time vector and 
% temperature record are provided by the user. Input should be provided as
% weekly averages.
%
% Inputs:   t       - time vector
%           temp    - temperature time series
%           MMMin   - max of the mean monthly climatology
%                     MMMin = 0 - use internally computed climatology
%                     MMMin = 1 - use CRW-NOAA Palau value
%                     MMMin > 1 - use inputted MMMin as the value
%           sstStd  - standard deviation of SST to be used to normalize var
%
% Outputs:  DHW     - Degree Heating Weeks
%           MMM     - max of the mean monthly climatology
%           MMC     - mean monthly climatology
%           degT    - degree threshold, either 1deg (no sstStd inputted) 
%                     or std based (sstStd included in input)
%
% Author: T.Schramek & M. Merrifield
% Date: 2017.11.30
% Edits:
% 2017.12.15 - TS - integrating MMC edits, interp to daily values 
% 2018.04.04 - TS - constrain MMC computation to years 1985 - 2012, per CRW
%--------------------------------------------------------------------------

%interpolate to daily values
tday = (t(1):t(end))';
k = find(diff(t)~=0);     %flag to avoid duplicate times in the record
Tday = interp1(t(k),temp(k),tday,'linear');
  
% make a monthly climatology
MMC = monthlyClimatology(tday,Tday);

% check for MMM input (0 = use internal calc, 1 = use CRW, >1 = use input)
if nargin > 2 % check if it was inputted
    if varargin{1} > 1
        MMM = varargin{1};
    else
        if varargin{1} == 1
            % MMM from CRW https://coralreefwatch.noaa.gov/vs/data/palau.txt
            % OVERWRITES THE CLIMATOLOGY
            MMM = 29.2309;
        end
        if varargin{1} == 0
            MMM = max(MMC); % max monthly climalology value
        end
    end
end


% if the std(stt) is inputted then use that to determine characteristic
% temperature scale to compute the threshold off of. else, use 1degC.
if nargin > 3
    % degThreshold = temperature threshold for bleaching threat
    % also can normalize it by a sst standard deviation
    sstStd = varargin{2};
    degT = (1/sstStd)*nanstd(Tday);
else
    % % OVERWRITES THE std normalized values
    degT = 1;
end

%DHW = values that >= MMM+degT for previous 12 weeks = 84 days
DHW = NaN(length(tday)-83,1);
tDHW = tday(84:end);
for j = 1:length(tday)-83 
    jj = j:j+83; 
    hotspot = Tday(jj) - MMM;
    k = find(hotspot >= degT); % only selecting hotspot >= 1
    if ~isempty(k)
         DHW(j) = nansum(hotspot(k))/7;
    else
        DHW(j) = 0;
    end
end

function MMC = monthlyClimatology(t,temp)
% T_climatology = monthlyClimatology(t,temp);
%--------------------------------------------------------------------------
% Description: load temp data. compute a monthly climatology
%
% Author: M.Merrifield
% Date: 2017.12.13
% Edits:
% 2018.04.08 - TS - limit climatology to 1985-2012, per CRW base period
%--------------------------------------------------------------------------

[yy,mo,~] = datevec(t);
% make a monthly climatology
for imonth = 1:12
    MMC(imonth) = nanmean(temp(mo==imonth & yy < 2013 & yy > 1984));
end
