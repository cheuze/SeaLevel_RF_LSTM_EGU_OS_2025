%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last step pre-machine learning, for the LSTM
% Copied a lot from dataprep_RF.m
% For each city:
% Cut all datasets to period common to all
% Min-max normalise
% Select the months with the most high sea level values
%
% I am tempted to ignore the rivers - RF showed they don't matter, and they
% are the ones with the shortest coverage in general (and dubious
% interpolation daily to hourly) - currently they're in
%
% Save the sea level values as "target" csv
% Save all the other variables as "predictors" csv
%
% Second section verifies the correlations between the predictors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all

pathdata='';
pathsave='';
load('levels_homogenised.mat','city_all')

for icity=1:length(city_all)

eval(sprintf('load([pathdata ''levels_homogenised.mat''],''%s'');',city_all{icity}));
eval(sprintf('Sealvl=%s;',city_all{icity}));eval(sprintf('clear %s;',city_all{icity}));

%removing outliers (threshold identified by eye)
Sealvl.Level_notide_notrend(find(abs(Sealvl.Level_notide_notrend)>2.5))=NaN;


%then set NaNs to zero - alt would be to set them to climatology
Sealvl.Level_notide_notrend(find(isnan(Sealvl.Level_notide_notrend)))=0;
% common years - not all are complete years, must do as for rivers
[jY,jM,jD]=datevec(Sealvl.Time);
yrZ(1,1)=jY(find(jM==1 & jD==1,1,'first'));
yrZ(1,2)=jY(find(jM==12 & jD==31,1,'last')); clear jY jM jD
%and there's ERA5: 1940-2023
yrZ(2,1)=1940; yrZ(2,2)=2023;

%rivers 
fileZ=dir([pathdata city_all{icity} '_rivers_*']);
timeriv=ncread([pathdata fileZ(1).name],'Datenum');
[jY,jM,jD]=datevec(timeriv);
yrZ(3,1)=jY(find(jM==1 & jD==1,1,'first'));
yrZ(3,2)=jY(find(jM==12 & jD==31,1,'last'));
clear jY jM jD fileZ timeriv
yr1=nanmax(yrZ(:,1),[],1);yr2=nanmin(yrZ(:,2),[],1); clear yrZ

% during the common years, months with most high sea levels
% do a rolling mean
pos=find(year(Sealvl.Time)>=yr1 & year(Sealvl.Time)<=yr2);
junk=Sealvl.Level_notide_notrend(pos); 
%do min-max normalisation here, on full series
junk=(junk-nanmin(junk))./(nanmax(junk)-nanmin(junk));

junktime=Sealvl.Datenum(pos); clear pos
for iwt=1:length(junk)-2*30*24 %cutting one month on either side
junk_mth(iwt,1)=nanmean(junk(30*24+iwt:60*24-1+iwt));
end
high3std=find(junk_mth>=nanmean(junk_mth)+3*nanstd(junk_mth))+30*24;
clear iwt junk_mth

%that's 2-3 thousands of points per city if I take them all 
high3std=high3std(~isnan(junk(high3std)));
Sealvl=junk(high3std); 

T=table(Sealvl);
writetable(T,[pathsave city_all{icity} '_target.csv']); 
clear T Sealvl junk

% Now for the other predictors
%load only from yr1 to yr2
yrera1=1979; yrera2=2019;
time_era5=datenum([yrera1,1,1]):1/24:datenum([yrera2,12,31]);
posera=find(time_era5>=datenum(yr1,1,1) & time_era5<=datenum(yr2,12,31));

e=ncread([pathdata city_all{icity} '_e_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'e_detrend',[posera(1)],[numel(posera)]);
msl=ncread([pathdata city_all{icity} '_msl_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'msl_detrend',[posera(1)],[numel(posera)]);
tp=ncread([pathdata city_all{icity} '_tp_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'tp_detrend',[posera(1)],[numel(posera)]);
sst=ncread([pathdata city_all{icity} '_sst_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'sst_detrend',[posera(1)],[numel(posera)]);
upos=ncread([pathdata city_all{icity} '_winds_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'upos',[posera(1)],[numel(posera)]);
uneg=ncread([pathdata city_all{icity} '_winds_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'uneg',[posera(1)],[numel(posera)]);
vpos=ncread([pathdata city_all{icity} '_winds_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'vpos',[posera(1)],[numel(posera)]);
vneg=ncread([pathdata city_all{icity} '_winds_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'vneg',[posera(1)],[numel(posera)]);
wdir=ncread([pathdata city_all{icity} '_winds_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'wdir',[posera(1)],[numel(posera)]);
wspeed=ncread([pathdata city_all{icity} '_winds_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'wspeed',[posera(1)],[numel(posera)]);
mslazo=ncread([pathdata 'remote_mslNAO_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'msl_azo_detrend',[posera(1)],[numel(posera)]);
mslice=ncread([pathdata 'remote_mslNAO_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'msl_ice_detrend',[posera(1)],[numel(posera)]);
sstena=ncread([pathdata 'remote_sstENA_' num2str(yrera1) '_' num2str(yrera2) '_hourly_ERA5.nc'],'sst_ena_detrend',[posera(1)],[numel(posera)]);
clear posera 

fileR=dir([pathdata city_all{icity} '_rivers_*.nc']);
time_rivers=ncread([pathdata fileR(1).name],'Datenum'); 
posriv=find(time_rivers>=datenum(yr1,1,1) & time_rivers<=datenum(yr2,12,31));
rivers=ncread([pathdata fileR(1).name],'runoff',posriv(1),numel(posriv));
%set NaN in rivers to 0 - again, could be more elegant, with 0 if always NaN and clim if only occasional missing value
rivers(isnan(rivers))=0;
time_era5=time_rivers(posriv);
clear time_rivers posriv fileR yr1 yr2

%I need to do the min-max normalising (and abs) before selection
parZ={'e';'msl';'sst';'tp';'upos';'uneg';'vpos';'vneg';'wdir';'wspeed';'mslazo';'mslice';'sstena';'rivers'};
for ipar=1:length(parZ)
eval(sprintf('%s=abs(%s);',parZ{ipar},parZ{ipar}));
if strcmp(parZ{ipar},'rivers')==0 || nanmax(rivers)~=nanmin(rivers) %if not river set to 0 everywhere
eval(sprintf('%s=(%s-nanmin(%s))./(nanmax(%s)-nanmin(%s));',parZ{ipar},parZ{ipar},parZ{ipar},parZ{ipar},parZ{ipar}));
end
end

%for all of the above, select simultaneous 
%mismatch for a reason I can't figure out - dirty loop then
for iday=1:length(high3std)
postime(iday)=find(time_era5==junktime(high3std(iday))); 
end
clear iday

for ipar=1:length(parZ)
    eval(sprintf('%s=%s(postime);',parZ{ipar},parZ{ipar}))
end

clear delayZ eventZ idel iev ipar time_era5 posEv parZ junktime high3std postime yrera*

C=who; pos=find(strcmp(C,'e')==0 & strcmp(C,'city_all')==0 & strcmp(C,'icity')==0 & strcmp(C,'pathsave')==0 & strcmp(C,'pathdata')==0);
C=C(pos); 
T=table(e); %initialise the table with one of the variables
for iC=1:length(C)
eval(sprintf('T=addvars(T,%s);',C{iC}))
eval(sprintf('clear %s',C{iC}))
end
writetable(T,[pathsave city_all{icity} '_predictors.csv'])
clear pos C T DD iC e


end

%% Verify the correlations

pathLSTMseries=''; %should be same as pathsave of previous section
city_all={'den_helder'; 'esbjerg'; 'gedser';'Got'; 'helsinki'; 'lowestoft'; 'KlagshamnMalmoe'; 'oslo'; 'RatanUmeaa'};

%load the list of all potential predictors
Tref=readtable([pathLSTEMseries 'lowestoft.csv']);
AllFeatures=Tref.Properties.VariableNames; clear Tref

mat_corr=NaN(14,14,length(city_all));

for icity=1:length(city_all)
    T=table2array(readtable([pathLSTEMseries city_all{icity} '_predictors.csv']));
    for iv1=1:13
        for iv2=iv1+1:14
            [RR,pp]=corrcoef(T(:,iv1),T(:,iv2));
            if pp(1,2)<0.05
                mat_corr(iv1,iv2,icity)=RR(1,2);
            end
            clear RR pp
        end
    end
    clear T
end