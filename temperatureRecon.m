function [Tout, Ttout, R2all] = temperatureRecon(SST_in, SSTt_in, SL_in, SLt_in, T_in, Tt_in)
% temperatureRecon.m
%--------------------------------------------------------------------------
% Description: This is example code for computing a temperature
% reconstruction from observations. All data are turned into weekly
% averages. SST, SLA and SLA^2 are used as inputs for the multiple linear
% regression 
%
% Inputs:   SST  - sea surface temperature record
%           SSTt - ^associated time vector
%           SL   - local sea level record
%           SLt  - ^associated time vector
%           Tin  - insitu temperature record 
%           Ttin - ^associated time vector
%
% Outputs:  Tout  - reconstruction of insitu temperature record
%           Ttout - ^associated time vector
%           R2all - R^2 of model fit and for each component separately
%                      where R2all = [full model; SST; SLA; SLA^2]
%
%--------------------------------------------------------------------------
% Author: T.Schramek
% Date: 2018.05.08
%--------------------------------------------------------------------------

% demean SL and SST
SLA_in = SL_in - nanmean(SL_in);
SST_in = SST_in - nanmean(SST_in);
T_m = nanmean(T_in);
T = T_in - T_m;

% determine weekly averages - if not already computed 
[SSTt,SST] = weeklyMean(SSTt_in,SST_in);
[SLAt,SLA] = weeklyMean(SLt_in,SLA_in);
[Tt,T] = weeklyMean(Tt_in,T);

% verify unique time vectors
[~,ia,~] = unique(SSTt);SSTt = SSTt(ia);SST = SST(ia);
[~,ia,~] = unique(SLAt);SLAt = SLAt(ia);SLA = SLA(ia);
[~,ia,~] = unique(Tt);Tt = Tt(ia);T = T(ia);

% tCommon - set common time step for SSTA and SLA
Ttout = round(min(SLAt):1:max(SLAt));
SSTi = interp1(SSTt,SST,Ttout,'linear');
SLAi = interp1(SLAt,SLA,Ttout,'linear');
Ti = interp1(Tt,T,Ttout,'linear');

% add in SLA^2 term to account for non-linear response for SLA
SLA2i = SLAi.^2;
SLA2i = SLA2i-nanmean(SLA2i);

% create linear regression model
regin = [SSTi(:) SLAi(:) SLA2i(:) ones(size(SLAi(:)))];
[b,~,~,~,stats] = regress(Ti',regin);

% reconstruct the time series
Tfulldm = regin*b;

% add mean back in
Tout = Tfulldm + T_m;

% grab stats
R2 = stats(1);

% determine fraction of variance (R^2) for each input
R2_inputs = zeros(size(size(regin,2)-1));
for ii = 1:size(regin,2)-1
    [~,~,~,~,statsall] = regress(Ti',regin(:,[ii,end]));
    R2_inputs(ii) = statsall(1);
end

% combine the R^2 values 
R2all = [R2,R2_inputs];

