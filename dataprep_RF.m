%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last step pre-machine learning, for the RF
% For each city:
% Cut all datasets to period common to all
% Obs: no need to normalise the data for RF 
% Select the X highest sea level events - have a time interval constraint
% Get the variables at that time, -1 h, -3h, -6h, -12h, -24h, -48h, -72h, -5d and -7d 
% Save the sea level values as "target" csv
% Save all the other variables as "predictors" csv
%
% Second section verifies the correlations between the predictors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all

pathdata='';
pathsave='';
city_all={'den_helder'; 'esbjerg'; 'gedser';'Got'; 'helsinki'; 'lowestoft'; 'KlagshamnMalmoe'; 'oslo'; 'RatanUmeaa'};

for icity=1:length(city_all)
    yrera1=1940; yrera2=2023;

eval(sprintf('load([pathdata ''levels_homogenised.mat''],''%s'');',city_all{icity}));
eval(sprintf('Sealvl=%s;',city_all{icity}));eval(sprintf('clear %s;',city_all{icity}));

% common years
yrZ(1,1)=nanmin(year(Sealvl.Time));yrZ(1,2)=nanmax(year(Sealvl.Time));
fileZ=dir([pathdata city_all{icity} '_rivers_*']);
timeriv=ncread([pathdata fileZ(1).name],'Datenum');
[jY,jM,jD]=datevec(timeriv);
yrZ(2,1)=jY(find(jM==1 & jD==1,1,'first'));
yrZ(2,2)=jY(find(jM==12 & jD==31,1,'last'));
clear jY jM jD timeriv fileZ
%and there's ERA5: 1940-2023
yrZ(3,1)=yrera1; yrZ(3,2)=yrera2;
yr1=nanmax(yrZ(:,1),[],1);yr2=nanmin(yrZ(:,2),[],1); clear yrZ

% during the common years, highest  sea level values
pos=find(year(Sealvl.Time)>=yr1 & year(Sealvl.Time)<=yr2);
junk=Sealvl.Level_notide_notrend(pos); junk(abs(junk)>2)=NaN;
high3std=find(junk>=nanmean(junk)+3*nanstd(junk)); %yields ca 4000 points
high3std=high3std(high3std>7*24); %issue if flood in very first week of series
% identify 7 day periods with more than 1 value; keep maximum sea level
junktime=Sealvl.Datenum(pos); clear pos
eventZ=[];
compt=1;
while compt<length(high3std)
    j1=junk(high3std(compt):min(high3std(compt)+7*24,length(junk)));
    j2=junktime(high3std(compt):min(high3std(compt)+7*24,length(junk)));
    [A,B]=max(j1);
    eventZ=[eventZ;A j2(B)]; clear A B j1 j2
    compt=find(high3std>high3std(compt)+7*24,1,'first');
    if isempty(compt)
        compt=length(high3std);
    end
end

clear compt junk junktime high3std Sealvl
Sealvl=eventZ(:,1); T=table(Sealvl);
writetable(T,[pathsave city_all{icity} '_target.csv']); clear T Sealvl

% Now for the predictors
%load only from yr1 to yr2
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
time_era5=time_rivers(posriv);
clear time_rivers posriv fileR yr1 yr2

%for all of the above, select simultaneous and up to 1 week before
for iev=1:size(eventZ,1)
posEv(iev,1)=find(time_era5==eventZ(iev,2));
end

delayZ=[0 1 3 6 12 24 48 72 5*24 7*24];
parZ={'e';'msl';'sst';'tp';'upos';'uneg';'vpos';'vneg';'wdir';'wspeed';'mslazo';'mslice';'sstena';'rivers'};
for ipar=1:length(parZ)
    for idel=1:length(delayZ)
        eval(sprintf('%s_h%.0f=%s(posEv-delayZ(idel));',parZ{ipar},delayZ(idel),parZ{ipar}))
    end
    eval(sprintf('clear %s',parZ{ipar}))
end

%time components
[YYYY,MM,DD]=datevec(eventZ(:,2));

clear delayZ eventZ idel iev ipar time_era5 posEv parZ yrera1 yrera2

C=who; pos=find(strcmp(C,'DD')==0 & strcmp(C,'city_all')==0 & strcmp(C,'icity')==0 & strcmp(C,'pathsave')==0 & strcmp(C,'pathdata')==0);
C=C(pos); 
T=table(DD);
for iC=1:length(C)
eval(sprintf('T=addvars(T,%s);',C{iC}))
eval(sprintf('clear %s',C{iC}))
end

writetable(T,[pathsave city_all{icity} '_predictors.csv'])
clear pos C T DD iC


end


%% Verify the correlation between the predictors

clear all; close all;

pathRFseries=''; %should be same as pathsave of previous section
load('levels_homogenised.mat','city_all')

matmat=NaN(14,14,length(city_all));

for icity=1:length(city_all)
    T=readtable([pathRFseries city_all{icity} '_predictors.csv']);
    if icity==1
        AllFeatures=T.Properties.VariableNames;
        for ifeat=1:14
        AllFeatures_leg{ifeat}=AllFeatures{4+(ifeat-1)*10};
        end
    end
    subZ=T(:,4:10:end);
    for ivar1=1:size(subZ,2)-1
        for ivar2=ivar1+1:size(subZ,2)
            [RR,pp]=corrcoef(table2array(subZ(:,ivar1)),table2array(subZ(:,ivar2)));
            if pp(1,2)<=0.05
                matmat(ivar1,ivar2,icity)=RR(1,2);
            end
            clear RR pp
        end
    end
    clear T subZ
end

