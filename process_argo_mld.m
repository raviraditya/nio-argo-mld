%% =======================================================================
%  process_argo_mld.m
%  North Indian Ocean Argo MLD threshold-sensitivity study.
%
%  Walks every Argo *_prof.nc file, applies quality control, computes
%  per-profile MLD (multiple density and temperature criteria), isothermal
%  layer depth (ILD) and barrier layer thickness (BLT), and writes ONE flat
%  per-profile table  ->  .mat (checkpoint) + .xlsx (the table used by the
%  figure script and shared in the repository).
%
%  Reference depth 10 dbar; criteria & BLT per de Boyer Montegut (2004) and
%  Sprintall & Tomczak (1992); density via TEOS-10 (GSW).
%
%  HOW TO RUN:
%    1. Set ARGO_ROOT to the folder containing the Argo YYYY/MM/*_prof.nc files.
%    2. Set OUT_DIR to where outputs should be written (default: ./output).
%    3. Set GSW_PATH to the installed Gibbs SeaWater (TEOS-10) toolbox folder.
%    4. Press the green Run button (F5), or type  process_argo_mld
%    Do NOT use Run Section -- the local function needs the whole file.
% =======================================================================

clear; clc;

%% ======================= PATHS / SETTINGS =============================
% --- EDIT THESE THREE PATHS FOR YOUR SYSTEM ---------------------------
ARGO_ROOT = fullfile('data','argo_indian_ocean');  % parent of YYYY/MM/*_prof.nc
OUT_DIR   = fullfile('output');                     % all outputs go here
GSW_PATH  = '';   % e.g. fullfile('toolbox','GSW-Matlab'); leave '' if already on path
% ----------------------------------------------------------------------

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

OUT_MAT  = fullfile(OUT_DIR,'nio_mld_profiles.mat');
OUT_XLSX = fullfile(OUT_DIR,'nio_mld_profiles.xlsx');

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

dsig_names = arrayfun(@(x) sprintf('MLDdens_%03d',round(x*1000)),DSIG_LIST,'uni',0);
dt_names   = arrayfun(@(x) sprintf('MLDtemp_%02d', round(x*10) ),DT_LIST,  'uni',0);

%% ========================= PROCESS ====================================
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
save(OUT_MAT,'T_tab','-v7.3');                       % checkpoint (not shared)
writetable(T_tab,OUT_XLSX,'FileType','spreadsheet'); % per-profile table (shared)
fprintf('Saved per-profile table (%d rows): .mat + .xlsx in %s\n',height(T_tab),OUT_DIR);

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
