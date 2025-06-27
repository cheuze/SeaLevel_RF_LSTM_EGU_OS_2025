%% Removing the tidal signal, using the UTide package
% Load the sea level data
% De-tides and de-trends them, using the Matlab package UTide (see paper)
% Converts them both to m (SWE in cm; GESLA/EU in m)
% Saves all nine locations into common .mat
% Keeps both "full" sea level series, and the de-tided de-trended one

%fill these two paths
pathrawdata='';
addpath([path 'UTideCurrentVersion'])

%the Swedish locations
city_Swe={'Got';'KlagshamnMalmoe';'RatanUmeaa'};
lat_Swe=[57.6847;55.5222;63.9861];
lon_Swe=[11.7906;12.8936;20.895];

%the other six
city_EU={'lowestoft';'den_helder';'esbjerg';'gedser';'helsinki';'oslo'};
lat_EU=[52.473083;52.964357;55.4667;54.567;60.153599;59.908559];
lon_EU=[1.75025;4.74499;8.4333;11.9330;24.95620;10.734510];

for iSwe=1:3
    fileZ=dir([pathrawdata 'smhi-*' city_Swe{iSwe} '.csv']);
AA=readtable([pathrawdata fileZ(1).name],'NumHeaderLines',6);
BB=timetable(AA.DatumTid_UTC_,AA.Havsvattenst_nd);
BB=retime(BB,'hourly'); %add missing time steps
AA=timetable2table(BB); clear BB 
AA.Datenum=datenum(AA.Time);
coef = ut_solv (AA.Datenum,AA.Var1,[],lat_Swe(iSwe),'auto','RunTimeDisp','nnn');
[ sl_fit, ~ ] = ut_reconstr(AA.Datenum,coef); 
AA=renamevars(AA,'Var1','SeaLevel');
AA.Level_notide_notrend=AA.SeaLevel-sl_fit;
%convert both to m
AA.SeaLevel=AA.SeaLevel./100;AA.Level_notide_notrend=AA.Level_notide_notrend./100;
eval(sprintf('%s=AA;',city_Swe{iSwe}))
if iSwe==1
eval(sprintf('save(''levels_homogenised.mat'',''%s'',''-v7.3'')',city_Swe{iSwe}))    
else
eval(sprintf('save(''levels_homogenised.mat'',''%s'',''-append'')',city_Swe{iSwe}))
end
clear AA fileZ coef sl_fit
eval(sprintf('clear %s',city_Swe{iSwe}))
end


for iEU=1:6
    fileZ=dir([pathrawdata city_EU{iEU} '*_GESLA']);
    AA=readtable([pathrawdata fileZ(1).name],'NumHeaderLines',41);
    %date and time in two separate columns
    [YY,MM,DD,~,~,~]=datevec(AA.Var1); 
    [~,~,~,hh,mmm,~]=datevec(AA.Var2);
    AA.Mergedate=datenum(YY,MM,DD,hh,mmm,0); clear YY MM DD hh mmm
    %values not in the right columns; there should only be 5 column btw
    SL5=AA.Var5; SL3=AA.Var3;
    SL3(isnan(SL3))=SL5(isnan(SL3)); clear SL5
    %missing values are at -99.9999 (bit stupid since it can be a SL value)
    SL3(SL3==-99.9999)=NaN;
    AA.SL=SL3; clear SL3
    BB=timetable(datetime(AA.Mergedate,'ConvertFrom','datenum'),AA.SL);

    %now becomes the same as the Swedish script
    BB=retime(BB,'hourly'); 
    AA=timetable2table(BB); clear BB
    AA.Datenum=datenum(AA.Time);
    coef = ut_solv (AA.Datenum,AA.Var1,[],lat_EU(iEU),'auto','RunTimeDisp','nnn');
    [ sl_fit, ~ ] = ut_reconstr(AA.Datenum,coef); 
    AA=renamevars(AA,'Var1','SeaLevel');
    AA.Level_notide_notrend=AA.SeaLevel-sl_fit;
    eval(sprintf('%s=AA;',city_EU{iEU}))
    eval(sprintf('save(''levels_homogenised.mat'',''%s'',''-append'')',city_EU{iEU}))
    clear AA fileZ coef sl_fit
    eval(sprintf('clear %s',city_EU{iEU}))
end

