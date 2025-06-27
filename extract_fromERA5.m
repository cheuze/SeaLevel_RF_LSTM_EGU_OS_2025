%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce the timeseries of the ERA5 variables
% Uses ERA5 grid cell nearest to lat/lon of the city, as per GESLA files
% Saves timeseries as netcdf
%
% Section 2: detrends the data, adds to the netcdf
%
% For MSLP and SST, also produces "remote" timeseries:
% min SLP over Iceland (lon -30:-10, lat 56:66)
% max SLP over Azores (lon -37:-22, lat 32:42)
% mean SST over eastern North Atlantic (lon -40:-10; lat 40:60)
% lat/lon limits somewhat arbitrary. For SST, dictated by what I downloaded
% - stopped too far east to include the Gulf Stream region
%
% Section 3: rivers, detrend 
% Section 4: and turn daily into hourly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all

%Set these two paths
pathERA=''; %Where your ERA5 files are
pathsave=''; %Where you want to save the extracted series

load('levels_homogenised.mat','lon_*','lat_*','city_*')

varZ={'e';'msl';'sst';'tp';'u10';'v10'};
coord_ice=[-30 -10; 56 66];
coord_azo=[-37 -22; 32 42];
coord_ena=[-40 -10; 40 60];

for ivar=1:length(varZ)
    fileZ=dir([pathERA varZ{ivar} '_*.nc']);
    for ifile=1:length(fileZ)
        lastyr(ifile)=str2double(fileZ(ifile).name(end-6:end-3));
    end
    [~,loadme]=sort(lastyr); clear lastyr

    lat=ncread([pathERA fileZ(1).name],'latitude');
    lon=ncread([pathERA fileZ(1).name],'longitude');

    for iSwe=1:length(lon_Swe)
        loadlat(iSwe,:)=find(abs(lat-lat_Swe(iSwe))==nanmin(abs(lat-lat_Swe(iSwe))),1,'first');
        loadlon(iSwe,:)=find(abs(lon-lon_Swe(iSwe))==nanmin(abs(lon-lon_Swe(iSwe))),1,'first');
    end
    for iEU=1:length(lon_EU)
        loadlat(iEU+length(lon_Swe),:)=find(abs(lat-lat_EU(iEU))==nanmin(abs(lat-lat_EU(iEU))),1,'first');
        loadlon(iEU+length(lon_Swe),:)=find(abs(lon-lon_EU(iEU))==nanmin(abs(lon-lon_EU(iEU))),1,'first');
    end

    if strcmp(varZ{ivar},'sst') %special case for Oslo because of the fjord
    loadlat(9,:)=28; loadlon(9,:)=203; 
    end

    for icity=1:length(loadlat)
        for ifile=loadme
            if exist('Par','var')==0
                Par=squeeze(ncread([pathERA fileZ(ifile).name],varZ{ivar},[loadlon(icity,1) loadlat(icity,1) 1],[1 1 Inf]));
            else
                Par=cat(1,Par,squeeze(ncread([pathERA fileZ(ifile).name],varZ{ivar},[loadlon(icity,1) loadlat(icity,1) 1],[1 1 Inf])));
            end
        end

        itry=1;
        while isempty(find(~isnan(Par),1)) && itry<=3
            clear Par
            for ifile=loadme
                if exist('Par','var')==0
                    Par=squeeze(ncread([pathERA fileZ(ifile).name],varZ{ivar},[loadlon(icity,1)+itry loadlat(icity,1)+itry 1],[1 1 Inf]));
                else
                    Par=cat(1,Par,squeeze(ncread([pathERA fileZ(ifile).name],varZ{ivar},[loadlon(icity,1)+itry loadlat(icity,1)+itry 1],[1 1 Inf])));
                end
            end
            if isempty(find(~isnan(Par),1))
                clear Par
                for ifile=loadme
                    if exist('Par','var')==0
                        Par=squeeze(ncread([pathERA fileZ(ifile).name],varZ{ivar},[loadlon(icity,1)-itry loadlat(icity,1)-itry 1],[1 1 Inf]));
                    else
                        Par=cat(1,Par,squeeze(ncread([pathERA fileZ(ifile).name],varZ{ivar},[loadlon(icity,1)-itry loadlat(icity,1)-itry 1],[1 1 Inf])));
                    end
                end
            end
            itry=itry+1;
        end
        clear itry
        

        if icity<=3
            nccreate([pathsave city_Swe{icity} '_' varZ{ivar} '_1940_2023_hourly_ERA5.nc'],varZ{ivar},'Dimensions',{'Time',length(Par)});
            ncwrite([pathsave city_Swe{icity} '_' varZ{ivar} '_1940_2023_hourly_ERA5.nc'],varZ{ivar},Par);
        else
            nccreate([pathsave city_EU{icity-4} '_' varZ{ivar} '_1940_2023_hourly_ERA5.nc'],varZ{ivar},'Dimensions',{'Time',length(Par)});
            ncwrite([pathsave city_EU{icity-4} '_' varZ{ivar} '_1940_2023_hourly_ERA5.nc'],varZ{ivar},Par);
        end

        clear Par
    end

    if strcmp(varZ{ivar},'msl')
        poslon=find(lon>=coord_ice(1,1) & lon<=coord_ice(1,2));
        poslat=find(lat>=coord_ice(2,1) & lat<=coord_ice(2,2));
        for ifile=loadme
        junk=ncread([pathERA fileZ(ifile).name],varZ{ivar},[poslon(1) poslat(1) 1],[numel(poslon) numel(poslat) Inf]);
        junk=reshape(junk,[size(junk,1)*size(junk,2) size(junk,3)]);
        if exist('msl_ice','var')
            msl_ice=cat(2,msl_ice,squeeze(nanmin(junk,[],1)));
        else
            msl_ice=squeeze(nanmin(junk,[],1));
        end
        clear junk
        end
        msl_ice=msl_ice(:);
        clear poslon poslat

        poslon=find(lon>=coord_azo(1,1) & lon<=coord_azo(1,2));
        poslat=find(lat>=coord_azo(2,1) & lat<=coord_azo(2,2));
        for ifile=loadme
        junk=ncread([pathERA fileZ(ifile).name],varZ{ivar},[poslon(1) poslat(1) 1],[numel(poslon) numel(poslat) Inf]);
        junk=reshape(junk,[size(junk,1)*size(junk,2) size(junk,3)]);
        if exist('msl_azo','var')
            msl_azo=cat(2,msl_azo,squeeze(nanmax(junk,[],1)));
        else
            msl_azo=squeeze(nanmax(junk,[],1));
        end
        clear junk 
        end
        msl_azo=msl_azo(:);
        clear poslon poslat

        nccreate([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_ice','Dimensions',{'Time',length(msl_ice)});
        nccreate([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_azo','Dimensions',{'Time',length(msl_ice)});

        ncwrite([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_ice',msl_ice);
        ncwrite([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_azo',msl_azo);

        clear msl_*

    elseif strcmp(varZ{ivar},'sst')
        poslon=find(lon>=coord_ena(1,1) & lon<=coord_ena(1,2));
        poslat=find(lat>=coord_ena(2,1) & lat<=coord_ena(2,2));
        for ifile=loadme
        junk=ncread([pathERA fileZ(ifile).name],varZ{ivar},[poslon(1) poslat(1) 1],[numel(poslon) numel(poslat) Inf]);
        junk=reshape(junk,[size(junk,1)*size(junk,2) size(junk,3)]);
        if exist('sst_ena','var')
            sst_ena=cat(2,sst_ena,squeeze(nanmean(junk,1)));
        else
            sst_ena=squeeze(nanmean(junk,1));
        end
        clear junk 
        end
        sst_ena=sst_ena(:);

        nccreate([pathsave 'remote_sstENA_1940_2023_hourly_ERA5.nc'],'sst_ena','Dimensions',{'Time',length(sst_ena)});
        ncwrite([pathsave 'remote_sstENA_1940_2023_hourly_ERA5.nc'],'sst_ena',sst_ena);

        clear sst_ena poslon poslat

    end

    clear loadlat loadlon lat lon loadme fileZ

end




%% Detrend
clear all; close all

varZ={'e';'msl';'sst';'tp';'siconc'};
pathsave='';

for ivar=1:length(varZ)
    fileZ=dir([pathsave '*_' varZ{ivar} '*.nc']);
    for ifile=1:length(fileZ)
        try
        Par=ncread([pathsave fileZ(ifile).name],varZ{ivar});
        Time=1:length(Par); X=[ones(length(Time),1) Time'];
        [b,~,~,~,stats]=regress(Par(:),X);
        if stats(3)<=0.05
            Par_trend=b(1)+b(2).*Time;
            Par_detrend=Par(:)-Par_trend(:)+Par(1);
        else
            Par_detrend=Par(:);
        end

        eval(sprintf('nccreate([pathsave fileZ(ifile).name],''%s_detrend'',''Dimensions'',{''Time'',length(Par_detrend)});',varZ{ivar}))
        eval(sprintf('ncwrite([pathsave fileZ(ifile).name],''%s_detrend'',Par_detrend);',varZ{ivar}))

        catch %remote files
            if strcmp(varZ{ivar},'msl')
                Par1=ncread([pathsave fileZ(ifile).name],'msl_ice');
                Par2=ncread([pathsave fileZ(ifile).name],'msl_azo');
                Time=1:length(Par1); X=[ones(length(Time),1) Time'];
                [b,~,~,~,stats]=regress(Par1(:),X);
                if stats(3)<=0.05
                    Par_trend=b(1)+b(2).*Time;
                    Par1_detrend=Par1(:)-Par_trend(:)+Par1(1);
                else
                    Par1_detrend=Par1(:);
                end
                clear b stats Par_trend
                [b,~,~,~,stats]=regress(Par2(:),X);
                if stats(3)<=0.05
                    Par_trend=b(1)+b(2).*Time;
                    Par2_detrend=Par2(:)-Par_trend(:)+Par2(1);
                else
                    Par2_detrend=Par2(:);
                end
                nccreate([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_ice_detrend','Dimensions',{'Time',length(Par1_detrend)});
                nccreate([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_azo_detrend','Dimensions',{'Time',length(Par1_detrend)});

                ncwrite([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_ice_detrend',Par1_detrend);
                ncwrite([pathsave 'remote_mslNAO_1940_2023_hourly_ERA5.nc'],'msl_azo_detrend',Par2_detrend);


            else %sst
                Par=ncread([pathsave fileZ(ifile).name],'sst_ena');
                Time=1:length(Par); X=[ones(length(Time),1) Time'];
                [b,~,~,~,stats]=regress(Par(:),X);
                if stats(3)<=0.05
                    Par_trend=b(1)+b(2).*Time;
                    Par_detrend=Par(:)-Par_trend(:)+Par(1);
                else
                    Par_detrend=Par(:);
                end
                nccreate([pathsave 'remote_sstENA_1940_2023_hourly_ERA5.nc'],'sst_ena_detrend','Dimensions',{'Time',length(Par_detrend)});
                ncwrite([pathsave 'remote_sstENA_1940_2023_hourly_ERA5.nc'],'sst_ena_detrend',Par_detrend);
 
            end
        end
        
    end
    clear fileZ
end


%% Rivers

clear all
close all

pathdata=''; %Where your river files are stored

fileZ=dir([pathdata '*_GRDC.nc']);
for ifile=1:length(fileZ)
    runoff=ncread([pathdata fileZ(ifile).name],'runoff_mean');
    runoff(runoff<0)=NaN; 

    time=ncread([pathdata fileZ(ifile).name],'time');
    [YYYY,MM,DD]=datevec(double(time)+datenum(1700,1,1));
    time=double(time)+datenum(1700,1,1); %as Matlab datenum

    runoff_seas=NaN(size(runoff));
    for iriv=1:size(runoff,1)
        % calculate the seasonal cycle
        % fill gap with seasonal cycle value
        % remove trend, if any
        % save as is
        % save sum of all rivers  
        
        for imth=1:12
            for iday=1:31
                pos=find(DD==iday & MM==imth);
                if ~isempty(pos)
                    runoff_seas(iriv,pos)=repmat(nanmedian(runoff(iriv,pos),'all'),[1 length(pos)]);
                end
                clear pos
            end
        end 
        pos=find(isnan(runoff(iriv,:)));
        if ~isempty(pos)
            runoff(iriv,pos)=runoff_seas(iriv,pos);
        end
        clear pos imth iday runoff_seas

        Time=1:length(time); X=[ones(length(Time),1) Time'];
        [b,~,~,~,stats]=regress(runoff(iriv,:)',X);
        if stats(3)<=0.05
            runoff_trend=b(1)+b(2).*Time;
            runoff_detrend(iriv,:)=runoff(iriv,:)-runoff_trend+runoff(iriv,1);
        else
            runoff_detrend(iriv,:)=runoff(iriv,:);
        end

        clear runoff_trend stats b Time
    end

    nccreate([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'Datenum','Dimensions',{'Time',length(time)});
    nccreate([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'runoff','Dimensions',{'Rivers',size(runoff,1),'Time',length(time)});
    nccreate([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'runoff_detrended','Dimensions',{'Rivers',size(runoff,1),'Time',length(time)});
    nccreate([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'runoff_detrended_sum','Dimensions',{'Time',length(time)});
    
    ncwrite([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'Datenum',time);
    ncwrite([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'runoff',runoff);
    ncwrite([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'runoff_detrended',runoff_detrend);
    junk=nansum(runoff_detrend,1);
    ncwrite([pathdata fileZ(ifile).name(1:end-7) 'rivers_' num2str(nanmin(YYYY)) '_' num2str(nanmax(YYYY)) '.nc'],'runoff_detrended_sum',junk(:));
     
    clear time runoff YYYY MM DD runoff_detrend junk
end

%% Rivers: turn daily sum into hourly
% Make files of 0 for cities with no river
clear all; close all;

pathdata='';
load('levels_homogenised.mat','city_all')

for ic=1:length(city_all)
    fileZ=dir([pathdata city_all{ic} '_rivers_*.nc']);
    if isempty(fileZ)
        time=datenum(1940,1,1):1/24:datenum(2023,12,31);
        runoff=zeros(length(time),1);
    nccreate([pathdata city_all{ic} '_rivers_1940_2023_hourly.nc'],'Datenum','Dimensions',{'Time',length(time)});
    nccreate([pathdata city_all{ic} '_rivers_1940_2023_hourly.nc'],'runoff','Dimensions',{'Time',length(time)});

    ncwrite([pathdata city_all{ic} '_rivers_1940_2023_hourly.nc'],'Datenum',time(:));
    ncwrite([pathdata city_all{ic} '_rivers_1940_2023_hourly.nc'],'runoff',runoff(:));

    else
        runoff_ds=ncread(fileZ(1).name,'runoff_detrended_sum');
        time_ds=ncread(fileZ(1).name,'Datenum');

        time=nanmin(time_ds):1/24:nanmax(time_ds);
        runoff=interp1(time_ds,runoff_ds,time);

    nccreate([pathdata fileZ(1).name(1:end-3) '_hourly.nc'],'Datenum','Dimensions',{'Time',length(time)});
    nccreate([pathdata fileZ(1).name(1:end-3) '_hourly.nc'],'runoff','Dimensions',{'Time',length(time)});

    ncwrite([pathdata fileZ(1).name(1:end-3) '_hourly.nc'],'Datenum',time(:));
    ncwrite([pathdata fileZ(1).name(1:end-3) '_hourly.nc'],'runoff',runoff(:));
    end

    clear time time_ds runoff runoff_ds fileZ
end

%% u10 and v10
% make and detrend wind speed, u10 pos, u10 neg etc

pathdata='';
load('levels_homogenised.mat','city_all')

yr1=1940; yr2=2023; %time period to consider

for icity=1:length(city_all)
    u=ncread([pathdata city_all{icity} '_u10_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'u10');
    v=ncread([pathdata city_all{icity} '_v10_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'v10');

    wspeed=sqrt(u.^2+v.^2);
    wdir=mod(180+atan2d(u,v),360); %ref: "from the north" = 0; eq by ecmwf

    upos=u; upos(upos<0)=0;uneg=u; uneg(uneg>0)=0;
    vpos=v; vpos(vpos<0)=0;vneg=v; vneg(vneg>0)=0;

    Time=1:length(u); X=[ones(length(Time),1) Time'];
    [b,~,~,~,stats]=regress(wspeed,X);
    if stats(3)<=0.05
        wspeed_trend=b(1)+b(2).*Time;
        wspeed_detrend=wspeed(:)-wspeed_trend(:)+wspeed(1);
    else
        wspeed_detrend=wspeed(:);
    end

    [b,~,~,~,stats]=regress(upos,X);
    if stats(3)<=0.05
        upos_trend=b(1)+b(2).*Time;
        upos_detrend=upos(:)-upos_trend(:)+upos(1);
    else
        upos_detrend=upos(:);
    end
    [b,~,~,~,stats]=regress(uneg,X);
    if stats(3)<=0.05
        uneg_trend=b(1)+b(2).*Time;
        uneg_detrend=uneg(:)-uneg_trend(:)+uneg(1);
    else
        uneg_detrend=uneg(:);
    end

    [b,~,~,~,stats]=regress(vpos,X);
    if stats(3)<=0.05
        vpos_trend=b(1)+b(2).*Time;
        vpos_detrend=vpos(:)-vpos_trend(:)+vpos(1);
    else
        vpos_detrend=vpos(:);
    end
    [b,~,~,~,stats]=regress(vneg,X);
    if stats(3)<=0.05
        vneg_trend=b(1)+b(2).*Time;
        vneg_detrend=vneg(:)-vneg_trend(:)+vneg(1);
    else
        vneg_detrend=vneg(:);
    end

    nccreate([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'vpos','Dimensions',{'Time',length(u)});
    nccreate([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'vneg','Dimensions',{'Time',length(u)});
    nccreate([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'upos','Dimensions',{'Time',length(u)});
    nccreate([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'uneg','Dimensions',{'Time',length(u)});
    nccreate([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'wspeed','Dimensions',{'Time',length(u)});
    nccreate([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'wdir','Dimensions',{'Time',length(u)});

    ncwrite([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'vpos',vpos_detrend);
    ncwrite([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'upos',upos_detrend);
    ncwrite([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'vneg',vneg_detrend);
    ncwrite([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'uneg',uneg_detrend);
    ncwrite([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'wspeed',wspeed_detrend);
    ncwrite([pathdata city_all{icity} '_winds_' num2str(yr1) '_' num2str(yr2) '_hourly_ERA5.nc'],'wdir',wdir);
    
clear upos uneg vpos vneg wspeed wdir *_detrend u v

end