% Example_temperatureRecon_DHW.m
%--------------------------------------------------------------------------
% Description: This is example code for computing a temperature
% reconstruction from observations throughout the water column
%
% Inputs:   SST  - sea surface temperature record
%           SSTt - ^associated time vector
%           SL   - local sea level record
%           SLt  - ^associated time vector
%           T    - insitu temperature record 
%           Tt   - ^associated time vector
%
%--------------------------------------------------------------------------
% Author: T.Schramek
% Date: 2018.05.08
%--------------------------------------------------------------------------
%% organize data

close all;clear;clc

examplepath = '.\ThermalStress';
addpath(genpath(examplepath));


%% SST

% load text file
hdrlns = 16; dl = ','; fname = 'short_drop_off.txt';
[sst,delimiterOut,headerlinesOut] = importdata(fname,dl,hdrlns);
sst.time = datenum(sst.data(:,1),sst.data(:,2),sst.data(:,3));
sst.sstd  = sst.data(:,4);
% compute weekly means 
[SSTt,SST] = weeklyMean(sst.time,sst.sstd);

%% SLA 

% SL - time
sl.time = ncread('OS_UH-FDD007_20170628_R.nc','time');
sl.time = sl.time+datenum(1700,1,1,0,0,0);
% SL - obs
sl.sl = ncread('OS_UH-FDD007_20170628_R.nc','sea_surface_height_above_reference_level');
sl.sl = squeeze(sl.sl)/1000;
% weekly averaged 
[sl.timew,sl.slw] = weeklyMean(sl.time,sl.sl);
% SLA - final form
SLAt = sl.timew;
SLA  = sl.slw-nanmean(sl.slw);

%% T - observations 

load Tobs_02mPalau

% weekly averaged 
[Ttobs,Tobs] = weeklyMean(Ttobs,Tobs);

%% Compute temperature reconstruction

[Trec, Ttrec, R2all] = temperatureRecon(SST, SSTt, SLA, SLAt, Tobs, Ttobs);

%% Compute Degree Heat Weeks for both T-recon and T-obs

MMMin = 1; % use hard coded MMM - value is for CRW Palau - virtual station
[DHWobs, tDHWobs, MMM, MMC, degT] = computeDHW(Ttobs,Tobs,MMMin);
[DHWrec, tDHWrec, MMM, MMC, degT] = computeDHW(Ttrec,Trec,MMMin);

%% plot results

trng = datenum(1997:4:2020,1,1);

figure,
% plot SST
subplot(411)
plot(SSTt,SST,'k-');
xlabel('Date (Jan-yyyy)');
ylabel('SST (\circC)');
xlim([trng(1) trng(end)]);set(gca,'XTick',trng);
datetick('x','yyyy','keepticks','keeplimits');
grid on;
ylim([25 32]);

% plot SLA
subplot(412)
plot(SLAt,SLA,'k-');
xlabel('Date (Jan-yyyy)');
ylabel('SLA (m)');
xlim([trng(1) trng(end)]);set(gca,'XTick',trng);
datetick('x','yyyy','keepticks','keeplimits');
grid on;

% plot temperature reconstructions vs observations 
subplot(413),
plot(Ttobs,Tobs,Ttrec,Trec);
xlim([trng(1) trng(end)]);set(gca,'XTick',trng);
ylim([25 32]);
xlabel('Date (Jan-yyyy)');
ylabel('T_{2m}(\circC)');
legend('Obs','Recon','location','best','orientation','horizontal');
datetick('x','yyyy','keepticks','keeplimits');
grid on;

% plot DHW time series - use both obs and recons
subplot(414),
plot(tDHWobs,DHWobs,tDHWrec,DHWrec);
xlim([trng(1) trng(end)]);set(gca,'XTick',trng);
xlabel('Date (Jan-yyyy)');
ylabel('DHW');
legend('Obs','Recon','location','best','orientation','horizontal');
datetick('x','yyyy','keepticks','keeplimits');
grid on;

