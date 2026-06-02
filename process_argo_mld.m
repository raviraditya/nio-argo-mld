%% =======================================================================
%  process_argo_mld.m
%  North Indian Ocean Argo MLD threshold-sensitivity study.
%
%  STAGE 1  Walk every *_prof.nc, compute per-profile MLD/ILD/BLT,
%           build ONE flat table  ->  .mat + .xlsx
%  STAGE 2  Apply physical caps, then write summary files:
%             (A) overall  min/mean/max/std  by ALL / AS / BoB   -> xlsx + txt
%             (B) region x season(4)  (AS-N/S, BoB-N/S)          -> xlsx + txt
%             (C) monsoon  PreMon / PostMon / Other  x group     -> xlsx + txt
%             (D) temp-criterion bias vs density + BLT corr       -> txt
%
%  Reference depth 10 dbar; criteria & BLT per de Boyer Montegut (2004),
%  Toyoda (2015), Jaffres (2013), Kara (2000),
%  Sprintall & Tomczak (1992); density via TEOS-10 (GSW).
%
%  HOW TO RUN:
%    1. Set ARGO_ROOT to the folder containing the Argo YYYY/MM/*_prof.nc files.
%    2. Set OUT_DIR to where you want outputs written (default: ./output).
%    3. Set GSW_PATH to the installed Gibbs SeaWater (TEOS-10) toolbox folder.
%    4. Press the green Run button (F5), or type  process_argo_mld
%    Do NOT use Run Section -- the local function needs the whole file.
%
%  To re-do only summaries without the slow loop: set RUN_STAGE1=false.
% =======================================================================

clear; clc;

RUN_STAGE1 = true;   % set false to skip processing and only rebuild summaries

%% ======================= PATHS / SETTINGS =============================
% --- EDIT THESE THREE PATHS FOR YOUR SYSTEM ---------------------------
ARGO_ROOT = fullfile('data','argo_indian_ocean');  % parent of YYYY/MM/*_prof.nc
OUT_DIR   = fullfile('output');                     % all outputs go here
GSW_PATH  = '';   % e.g. fullfile('toolbox','GSW-Matlab'); leave '' if already on path
% ----------------------------------------------------------------------

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

OUT_MAT  = fullfile(OUT_DIR,'nio_mld_profiles.mat');
OUT_XLSX = fullfile(OUT_DIR,'nio_mld_profiles.xlsx');

SUM_OVERALL_XLSX = fullfile(OUT_DIR,'nio_summary_overall.xlsx');
SUM_OVERALL_TXT  = fullfile(OUT_DIR,'nio_summary_overall.txt');
SUM_REGION_XLSX  = fullfile(OUT_DIR,'nio_summary_region_season.xlsx');
SUM_REGION_TXT   = fullfile(OUT_DIR,'nio_summary_region_season.txt');
SUM_MONSOON_XLSX = fullfile(OUT_DIR,'nio_summary_monsoon.xlsx');
SUM_MONSOON_TXT  = fullfile(OUT_DIR,'nio_summary_monsoon.txt');
SUM_BIAS_TXT     = fullfile(OUT_DIR,'nio_summary_bias.txt');

if ~isempty(GSW_PATH) && exist(GSW_PATH,'dir'), addpath(genpath(GSW_PATH)); end
assert(exist('gsw_SA_from_SP','file')==2, ...
    'GSW toolbox not on path. Set GSW_PATH or add it to the MATLAB path.');

% --- region clip: North Indian Ocean -----------------------------------
LAT_MIN=0; LAT_MAX=30; LON_MIN=30; LON_MAX=100;
LON_ASBOB_SPLIT = 78;     % <78E -> Arabian Sea, >=78E -> Bay of Bengal
LAT_NS_SPLIT    = 12;     % >=12N -> North, else South (independent of MLD)

% --- thresholds (literature-defensible) --------------------------------
REF_DEPTH = 10;                                  % dbar (de Boyer Montegut 2004)
DSIG_LIST = [0.01 0.02 0.03 0.04 0.05 0.125];    % kg m-3
DT_LIST   = [0.2 0.5 0.8];                        % deg C
ILD_DT    = 0.5;                                  % deg C (Sprintall & Tomczak)

% --- QC / screening ----------------------------------------------------
GOOD_QC=['1' '2']; MAX_REF_GAP=25; MIN_LEVELS=5; MAX_TOP_SPACING=25; DEEP_ENOUGH=200;

% --- physical sanity caps applied before summaries ---------------------
CAP_MLD = 300;   % dbar
CAP_ILD = 300;   % dbar
CAP_BLT = 100;   % dbar (abs)

dsig_names = arrayfun(@(x) sprintf('MLDdens_%03d',round(x*1000)),DSIG_LIST,'uni',0);
dt_names   = arrayfun(@(x) sprintf('MLDtemp_%02d', round(x*10) ),DT_LIST,  'uni',0);

%% ========================= STAGE 1: PROCESS ===========================
if RUN_STAGE1
    files = dir(fullfile(ARGO_ROOT,'**','*_prof.nc'));
    nF = numel(files);
    fprintf('Found %d profile files.\n', nF);

    REC = {}; t0=tic; nSeen=0; nKept=0;

    for k=1:nF
        fp = fullfile(files(k).folder, files(k).name);
        try
            P_all=ncread(fp,'PRES');  T_all=ncread(fp,'TEMP');  S_all=ncread(fp,'PSAL');
            Pa=ncread(fp,'PRES_ADJUSTED'); Ta=ncread(fp,'TEMP_ADJUSTED'); Sa=ncread(fp,'PSAL_ADJUSTED');
            Pq=ncread(fp,'PRES_ADJUSTED_QC'); Tq=ncread(fp,'TEMP_ADJUSTED_QC'); Sq=ncread(fp,'PSAL_ADJUSTED_QC');
            Pqr=ncread(fp,'PRES_QC'); Tqr=ncread(fp,'TEMP_QC'); Sqr=ncread(fp,'PSAL_QC');
            LAT=ncread(fp,'LATITUDE'); LON=ncread(fp,'LONGITUDE'); JULD=ncread(fp,'JULD');
            DM=ncread(fp,'DATA_MODE'); WMOc=ncread(fp,'PLATFORM_NUMBER'); CYC=ncread(fp,'CYCLE_NUMBER');
            POSq=ncread(fp,'POSITION_QC'); JULq=ncread(fp,'JULD_QC');
            nprof=numel(LAT);
        catch ME
            warning('skip %s (%s)',files(k).name,ME.message); continue
        end

        for ip=1:nprof
            nSeen=nSeen+1;
            lat=LAT(ip); lon=LON(ip);
            if isnan(lat)||isnan(lon), continue; end
            if lon>180, lon=lon-360; end
            if lat<LAT_MIN||lat>LAT_MAX||lon<LON_MIN||lon>LON_MAX, continue; end
            if ~any(POSq(ip)==GOOD_QC)||~any(JULq(ip)==GOOD_QC), continue; end

            dm=DM(ip); useAdj=(dm=='D'||dm=='A');
            if useAdj
                P=Pa(:,ip);T=Ta(:,ip);S=Sa(:,ip); pq=Pq(:,ip);tq=Tq(:,ip);sq=Sq(:,ip);
                if all(isnan(P))||all(P==99999)
                    P=P_all(:,ip);T=T_all(:,ip);S=S_all(:,ip); pq=Pqr(:,ip);tq=Tqr(:,ip);sq=Sqr(:,ip); useAdj=false;
                end
            else
                P=P_all(:,ip);T=T_all(:,ip);S=S_all(:,ip); pq=Pqr(:,ip);tq=Tqr(:,ip);sq=Sqr(:,ip);
            end

            good=ismember(pq,GOOD_QC)&ismember(tq,GOOD_QC)&ismember(sq,GOOD_QC);
            P(P==99999)=NaN;T(T==99999)=NaN;S(S==99999)=NaN;
            good=good&~isnan(P)&~isnan(T)&~isnan(S)&P>=0&P<=12000;
            P=P(good);T=T(good);S=S(good);
            if numel(P)<MIN_LEVELS, continue; end
            [P,si]=sort(P);T=T(si);S=S(si);
            [P,ui]=unique(P,'stable');T=T(ui);S=S(ui);
            if numel(P)<MIN_LEVELS, continue; end
            if P(1)>MAX_REF_GAP||P(end)<DEEP_ENOUGH, continue; end
            top=P(P<=50); if numel(top)>=2&&max(diff(top))>MAX_TOP_SPACING, continue; end

            SA=gsw_SA_from_SP(S,P,lon,lat); CT=gsw_CT_from_t(SA,T,P); sig0=gsw_sigma0(SA,CT);

            Tref=interp1(P,T,REF_DEPTH,'linear'); Sref=interp1(P,sig0,REF_DEPTH,'linear');
            if isnan(Tref)||isnan(Sref), Tref=T(1); Sref=sig0(1); end

            mld_dens=nan(1,numel(DSIG_LIST));
            for j=1:numel(DSIG_LIST), mld_dens(j)=first_crossing(P,sig0,Sref+DSIG_LIST(j),REF_DEPTH); end
            mld_temp=nan(1,numel(DT_LIST));
            for j=1:numel(DT_LIST),   mld_temp(j)=first_crossing(P,abs(T-Tref),DT_LIST(j),REF_DEPTH); end
            ild=first_crossing(P,abs(T-Tref),ILD_DT,REF_DEPTH);
            mld03=mld_dens(DSIG_LIST==0.03); blt=ild-mld03;

            if lon<LON_ASBOB_SPLIT, basin='AS'; else, basin='BoB'; end
            if lat>=LAT_NS_SPLIT, ns='N'; else, ns='S'; end

            r=struct();
            r.file=files(k).name; r.wmo=strtrim(WMOc(:,ip)'); r.cycle=double(CYC(ip));
            r.juld=JULD(ip); r.datenum=datenum(1950,1,1)+JULD(ip);
            r.datestr=datestr(r.datenum,'yyyy-mm-dd');
            r.lat=lat; r.lon=lon; r.basin=basin; r.northsouth=ns;
            r.region=[basin '-' ns]; r.data_mode=dm;
            r.n_levels_used=numel(P); r.p_shallow=P(1); r.p_deepest=P(end);
            r.sst=T(1); r.sss=S(1); r.ref_depth=REF_DEPTH;
            for j=1:numel(DSIG_LIST), r.(dsig_names{j})=mld_dens(j); end
            for j=1:numel(DT_LIST),   r.(dt_names{j})  =mld_temp(j); end
            r.ild=ild; r.blt=blt;
            REC{end+1}=r; nKept=nKept+1; %#ok<SAGROW>
        end
        if mod(k,250)==0, fprintf('  %d/%d files | %d kept | %.1f min\n',k,nF,nKept,toc(t0)/60); end
    end

    fprintf('\n%d seen, %d kept (%.1f%%).\n',nSeen,nKept,100*nKept/max(nSeen,1));
    T_tab = struct2table([REC{:}]);
    save(OUT_MAT,'T_tab','-v7.3');
    writetable(T_tab,OUT_XLSX,'FileType','spreadsheet');
    fprintf('Saved full table (%d rows): .mat + .xlsx\n',height(T_tab));
end

%% ====================== STAGE 2: SUMMARIES ============================
load(OUT_MAT,'T_tab'); T=T_tab;

dsig_cols=T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'MLDdens_'));
dt_cols  =T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'MLDtemp_'));
vars_all =[dsig_cols, dt_cols, {'ild','blt'}];        % everything for overall
vars_key ={'MLDdens_030','MLDtemp_02','ild','blt'};   % key set for region/season
REFCOL='MLDdens_030';

% --- season labels: 4-bin and 3-bin -----------------------------------
mo=month(datetime(string(T.datestr),'InputFormat','yyyy-MM-dd'));
s4=strings(height(T),1);
s4(ismember(mo,[12 1 2]))="Winter"; s4(ismember(mo,[3 4 5]))="PreMon";
s4(ismember(mo,[6 7 8 9]))="Monsoon"; s4(ismember(mo,[10 11]))="PostMon";
T.season4=s4;
s3=strings(height(T),1);
s3(ismember(mo,[3 4 5]))="PreMon"; s3(ismember(mo,[10 11]))="PostMon";
s3(s3=="")="Other";                     % Winter+Monsoon
T.season3=s3;

% --- physical sanity caps (remove junk profiles from stats) -----------
mldcols = [dsig_cols, dt_cols, {'ild'}];
bad = abs(T.blt) > CAP_BLT;
for c = 1:numel(mldcols), bad = bad | T.(mldcols{c}) > CAP_MLD; end
T(bad,:) = [];
fprintf('Removed %d non-physical profiles; %d remain for summaries.\n',sum(bad),height(T));

regions={'AS-N','AS-S','BoB-N','BoB-S'};

%% (A) OVERALL min/mean/max/std  by ALL / AS / BoB ---------------------
grp={'ALL','AS','BoB'};
Aout={};
for g=1:numel(grp)
    if g==1, m=true(height(T),1); else, m=strcmp(T.basin,grp{g}); end
    row={grp{g},sum(m)};
    for v=1:numel(vars_all)
        x=T.(vars_all{v})(m); x=x(~isnan(x));
        if isempty(x),mn=NaN;me=NaN;mx=NaN;sd=NaN;
        else,mn=min(x);me=mean(x);mx=max(x);sd=std(x);end
        row=[row,{round(mn,1),round(me,1),round(mx,1),round(sd,1)}]; %#ok<AGROW>
    end
    Aout(end+1,:)=row; %#ok<SAGROW>
end
hdrA={'Group','N'};
for v=1:numel(vars_all), hdrA=[hdrA,{[vars_all{v} '_min'],[vars_all{v} '_mean'],[vars_all{v} '_max'],[vars_all{v} '_std']}]; end %#ok<AGROW>
writetable(cell2table(Aout,'VariableNames',hdrA),SUM_OVERALL_XLSX,'Sheet','overall');
fid=fopen(SUM_OVERALL_TXT,'w');
fprintf(fid,'OVERALL summary (min/mean/max/std, dbar). N profiles = %d\n\n',height(T));
for g=1:size(Aout,1)
    fprintf(fid,'== %s  (N=%d) ==\n',Aout{g,1},Aout{g,2});
    c=3;
    for v=1:numel(vars_all)
        fprintf(fid,'  %-12s  min %6.1f  mean %6.1f  max %6.1f  std %5.1f\n', ...
            vars_all{v},Aout{g,c},Aout{g,c+1},Aout{g,c+2},Aout{g,c+3}); c=c+4;
    end
    fprintf(fid,'\n');
end
fclose(fid);

%% (B) REGION x SEASON(4)  min/mean/max ---------------------------------
seasons4={'Winter','PreMon','Monsoon','PostMon','ALL'};
Bout={};
fid=fopen(SUM_REGION_TXT,'w');
fprintf(fid,'REGION x SEASON summary (min/mean/max, dbar). N/S split at %d N.\n\n',LAT_NS_SPLIT);
for r=1:numel(regions)
    fprintf(fid,'===== %s =====\n',regions{r});
    for s=1:numel(seasons4)
        if strcmp(seasons4{s},'ALL'), m=strcmp(T.region,regions{r});
        else, m=strcmp(T.region,regions{r})&strcmp(T.season4,seasons4{s}); end
        fprintf(fid,'  %-8s N=%6d |',seasons4{s},sum(m));
        row={regions{r},seasons4{s},sum(m)};
        for v=1:numel(vars_key)
            x=T.(vars_key{v})(m); x=x(~isnan(x));
            if isempty(x),mn=NaN;me=NaN;mx=NaN;else,mn=min(x);me=mean(x);mx=max(x);end
            fprintf(fid,' %s %5.1f/%5.1f/%5.1f',vars_key{v},mn,me,mx);
            row=[row,{round(mn,1),round(me,1),round(mx,1)}]; %#ok<AGROW>
        end
        fprintf(fid,'\n'); Bout(end+1,:)=row; %#ok<SAGROW>
    end
    fprintf(fid,'\n');
end
fclose(fid);
hdrB={'Region','Season','N'};
for v=1:numel(vars_key), hdrB=[hdrB,{[vars_key{v} '_min'],[vars_key{v} '_mean'],[vars_key{v} '_max']}]; end %#ok<AGROW>
writetable(cell2table(Bout,'VariableNames',hdrB),SUM_REGION_XLSX,'Sheet','region_season');

%% (C) MONSOON: PreMon / PostMon / Other  x group ----------------------
mons={'PreMon','PostMon','Other','ALL'};
groupsM=[{'ALL','AS','BoB'},regions];
Cout={};
fid=fopen(SUM_MONSOON_TXT,'w');
fprintf(fid,'MONSOON-SEASON summary (PreMon=MAM, PostMon=ON, Other=Winter+Monsoon).\n');
fprintf(fid,'min/mean/max, dbar.\n\n');
for g=1:numel(groupsM)
    gn=groupsM{g};
    fprintf(fid,'===== %s =====\n',gn);
    for s=1:numel(mons)
        if strcmp(gn,'ALL'), base=true(height(T),1);
        elseif any(strcmp(gn,{'AS','BoB'})), base=strcmp(T.basin,gn);
        else, base=strcmp(T.region,gn); end
        if strcmp(mons{s},'ALL'), m=base; else, m=base&strcmp(T.season3,mons{s}); end
        fprintf(fid,'  %-8s N=%6d |',mons{s},sum(m));
        row={gn,mons{s},sum(m)};
        for v=1:numel(vars_key)
            x=T.(vars_key{v})(m); x=x(~isnan(x));
            if isempty(x),mn=NaN;me=NaN;mx=NaN;else,mn=min(x);me=mean(x);mx=max(x);end
            fprintf(fid,' %s %5.1f/%5.1f/%5.1f',vars_key{v},mn,me,mx);
            row=[row,{round(mn,1),round(me,1),round(mx,1)}]; %#ok<AGROW>
        end
        fprintf(fid,'\n'); Cout(end+1,:)=row; %#ok<SAGROW>
    end
    fprintf(fid,'\n');
end
fclose(fid);
hdrC={'Group','MonsoonSeason','N'};
for v=1:numel(vars_key), hdrC=[hdrC,{[vars_key{v} '_min'],[vars_key{v} '_mean'],[vars_key{v} '_max']}]; end %#ok<AGROW>
writetable(cell2table(Cout,'VariableNames',hdrC),SUM_MONSOON_XLSX,'Sheet','monsoon');

%% (D) BIAS of temp criterion vs density + BLT correlation -------------
ref=T.(REFCOL);
fid=fopen(SUM_BIAS_TXT,'w');
fprintf(fid,'Temperature-criterion (DeltaT=0.2) bias vs %s, and bias~BLT.\n\n',REFCOL);
fprintf(fid,'%-6s %9s %12s %8s %9s\n','Group','meanBias','r(bias,BLT)','slope','N');
for g=1:numel(grp)
    if g==1, m=true(height(T),1); else, m=strcmp(T.basin,grp{g}); end
    m=m&~isnan(T.blt)&~isnan(T.MLDtemp_02)&~isnan(ref);
    bias=T.MLDtemp_02(m)-ref(m); blt=T.blt(m);
    R=corrcoef(blt,bias); rr=R(1,2); p=polyfit(blt,bias,1);
    fprintf(fid,'%-6s %9.2f %12.2f %8.3f %9d\n',grp{g},mean(bias),rr,p(1),sum(m));
end
fclose(fid);

fprintf('\nAll summaries written to %s:\n',OUT_DIR);
fprintf('  overall : nio_summary_overall.(xlsx/txt)\n');
fprintf('  region  : nio_summary_region_season.(xlsx/txt)\n');
fprintf('  monsoon : nio_summary_monsoon.(xlsx/txt)\n');
fprintf('  bias    : nio_summary_bias.txt\n');

%% ========================= LOCAL FUNCTION =============================
function z=first_crossing(P,V,thresh,minz)
    if nargin<4, minz=0; end
    z=NaN; above=V<thresh;
    if all(above), return; end
    idx=find(~above,1,'first');
    if idx==1, z=NaN; return; end          % gradient at/above first level -> invalid
    v1=V(idx-1);v2=V(idx);p1=P(idx-1);p2=P(idx);
    if v2==v1, z=p2; else, z=p1+(thresh-v1)*(p2-p1)/(v2-v1); end
    if z<minz, z=NaN; end                  % shallower than reference depth -> invalid
end