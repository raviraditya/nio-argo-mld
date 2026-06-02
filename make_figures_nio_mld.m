%% =======================================================================
%  make_figures_nio_mld.m
%  Figures + Table 1 for the NIO Argo MLD threshold study.
%  Reads nio_mld_profiles.xlsx (output of process_argo_mld.m).
%
%  Fig 1  Argo profile-density map + AS/BoB split
%  Fig 2  Example profiles: (a) BoB barrier layer, (b) AS deep mixing
%  Fig 3  Seasonal (4-bin) sub-region means: MLD0.03, BLT, dT=0.2 bias
%  Fig 4  Bias & RMSE of each criterion vs density 0.03, by sub-basin
%  Table 1  error metrics: criterion x sub-basin x season -> CSV
%
%  Reference MLD = density at 0.03 kg m-3 (de Boyer Montegut 2004).
%
%  HOW TO RUN: set IN_TABLE, OUT_DIR, ARGO_ROOT, GSW_PATH below; press Run.
% =======================================================================

clear; clc; close all;

%% ---------------------- SETTINGS (EDIT) ------------------------------
IN_TABLE  = fullfile('output','nio_mld_profiles.xlsx');  % per-profile table
OUT_DIR   = fullfile('output','figures');                % figure output folder
ARGO_ROOT = fullfile('data','argo_indian_ocean');        % needed for Fig 2
GSW_PATH  = '';   % GSW toolbox folder; '' if already on path
% ----------------------------------------------------------------------

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
if ~isempty(GSW_PATH) && exist(GSW_PATH,'dir'), addpath(genpath(GSW_PATH)); end

REFCOL = 'MLDdens_030';

% Okabe-Ito colour-blind-safe palette
C.blue=[0 114 178]/255; C.orange=[230 159 0]/255; C.green=[0 158 115]/255;
C.verm=[213 94 0]/255;  C.purple=[204 121 167]/255; C.sky=[86 180 233]/255;
C.yellow=[240 228 66]/255; C.grey=[120 120 120]/255;
PAL=[C.blue;C.orange;C.green;C.verm;C.purple;C.sky;C.yellow;C.grey];

set(0,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
      'DefaultAxesFontSize',9,'DefaultAxesLineWidth',0.8, ...
      'DefaultAxesBox','on','DefaultAxesLayer','top');

%% ---------------------- LOAD + CAPS ----------------------------------
T = readtable(IN_TABLE);
T(T.MLDdens_030>300 | T.ild>300 | abs(T.blt)>100, :) = [];
fprintf('Profiles after physical cap: %d\n', height(T));

dsig_cols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'MLDdens_'));
dt_cols   = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'MLDtemp_'));
all_crit  = [dsig_cols, dt_cols];
critLabel = cellfun(@local_label, all_crit, 'uni',0);

mo = month(datetime(T.datestr,'InputFormat','yyyy-MM-dd'));
seas = strings(height(T),1);
seas(ismember(mo,[12 1 2]))="Winter(DJF)"; seas(ismember(mo,[3 4 5]))="PreMon(MAM)";
seas(ismember(mo,[6 7 8 9]))="Monsoon(JJAS)"; seas(ismember(mo,[10 11]))="PostMon(ON)";
T.season = seas;

basins  = {'AS','BoB'};
seasons = {'Winter(DJF)','PreMon(MAM)','Monsoon(JJAS)','PostMon(ON)'};
ref = T.(REFCOL);

%% ---------------------- FIG 1: map -----------------------------------
f1 = figure('Units','centimeters','Position',[2 2 11 7],'Color','w');
ax1 = axes(f1); hold(ax1,'on');
loe=30:0.5:100; lae=0:0.5:30;
N=histcounts2(T.lon,T.lat,loe,lae); N(N==0)=NaN;
pcolor(ax1,loe(1:end-1)+0.25,lae(1:end-1)+0.25,log10(N')); shading(ax1,'flat');
colormap(ax1,local_seq(C.blue));
cb=colorbar(ax1); cb.Label.String='log_{10}(profiles per 0.5\circ cell)';
if exist('coastlines','file')
    load coastlines coastlat coastlon
    plot(ax1,coastlon,coastlat,'k-','LineWidth',0.6);
end
xline(ax1,78,'--','Color',C.grey,'LineWidth',1.2);
text(ax1,50,28.5,'Arabian Sea','FontWeight','bold','Color',C.verm);
text(ax1,80,28.5,'Bay of Bengal','FontWeight','bold','Color',C.green);
xlabel(ax1,'Longitude (\circE)'); ylabel(ax1,'Latitude (\circN)');
xlim(ax1,[30 100]); ylim(ax1,[0 30]); daspect(ax1,[1 1 1]);
local_panel(ax1,'a');
local_export(f1,fullfile(OUT_DIR,'Fig1_map'),true);

%% ---------------------- FIG 2: example profiles ----------------------
isB=strcmp(T.basin,'BoB')&~isnan(T.blt); [~,EX_BOB]=max(T.blt.*isB);
isA=strcmp(T.basin,'AS')&~isnan(ref); tmp=ref; tmp(~isA)=NaN; [~,EX_AS]=max(tmp);
f2 = figure('Units','centimeters','Position',[2 2 16 8],'Color','w');
for pp=1:2
    if pp==1, rr=EX_BOB; lab='a'; ttl='Bay of Bengal'; cc=C.green;
    else,     rr=EX_AS;  lab='b'; ttl='Arabian Sea';  cc=C.verm; end
    [P,Tt,~,sig0]=local_reload_profile(T(rr,:),ARGO_ROOT);
    axL=subplot(1,2,pp); hold(axL,'on');
    plot(axL,Tt,P,'-','Color',cc,'LineWidth',1.6);
    set(axL,'YDir','reverse'); ylim(axL,[0 200]);
    xlabel(axL,'Temperature (\circC)'); if pp==1, ylabel(axL,'Pressure (dbar)'); end
    axT=axes('Position',axL.Position,'Color','none','XAxisLocation','top', ...
             'YAxisLocation','right','YDir','reverse'); hold(axT,'on'); ylim(axT,[0 200]);
    plot(axT,sig0,P,'-','Color',C.grey,'LineWidth',1.2);
    axT.YTickLabel=[]; xlabel(axT,'\sigma_\theta (kg m^{-3})','Color',C.grey);
    for j=1:numel(all_crit)
        yline(axL,T.(all_crit{j})(rr),':','Color',PAL(mod(j-1,size(PAL,1))+1,:),'LineWidth',1.0);
    end
    yline(axL,T.ild(rr),'-.','Color',C.orange,'LineWidth',1.3);
    text(axL,min(Tt)+0.3,8,ttl,'FontWeight','bold','Color',cc);
    local_panel(axL,lab);
end
legend(axL,[critLabel,{'ILD (\DeltaT=0.5)'}],'Location','eastoutside','FontSize',7,'Box','off');
local_export(f2,fullfile(OUT_DIR,'Fig2_profiles'));

%% ---------------------- FIG 3: seasonal panels -----------------------
regions5={'AS-N','AS-S','BoB-N','BoB-S'}; regPAL=[C.verm;C.orange;C.green;C.sky];
mldM=nan(4,4); bltM=nan(4,4); biasM=nan(4,4);
for s=1:4, for r=1:4
    m=strcmp(T.region,regions5{r})&strcmp(T.season,seasons{s});
    x=ref(m); x=x(~isnan(x)); if ~isempty(x), mldM(s,r)=mean(x); end
    b=T.blt(m); b=b(~isnan(b)); if ~isempty(b), bltM(s,r)=mean(b); end
    d=T.MLDtemp_02(m)-ref(m); d=d(~isnan(d)); if ~isempty(d), biasM(s,r)=mean(d); end
end, end
seasTick={'Win','Pre','Mon','Post'};
f3 = figure('Units','centimeters','Position',[2 2 17 16],'Color','w');
mats={mldM,bltM,biasM}; ylabs={'Mean MLD_{0.03} (dbar)','Mean BLT (dbar)','\DeltaT=0.2 bias (dbar)'};
for q=1:3
    ax=subplot(3,1,q); h=bar(ax,mats{q},'grouped');
    for r=1:4, h(r).FaceColor=regPAL(r,:); end
    set(ax,'XTick',1:4,'XTickLabel',seasTick); ylabel(ax,ylabs{q});
    if q>1, yline(ax,0,'k-'); end
    if q==1, legend(ax,regions5,'Box','off','Location','northwest','FontSize',7,'NumColumns',2); end
    if q==3, xlabel(ax,'Season'); end
    local_panel(ax,char('a'+q-1));
end
local_export(f3,fullfile(OUT_DIR,'Fig3_seasonal'));

%% ---------------------- FIG 4: bias & RMSE ---------------------------
crit_use=setdiff(all_crit,REFCOL,'stable');
critUseLab=cellfun(@local_label,crit_use,'uni',0);
bias=nan(numel(crit_use),2); rmse=nan(numel(crit_use),2);
for b=1:2
    m=strcmp(T.basin,basins{b});
    for j=1:numel(crit_use)
        d=T.(crit_use{j})(m)-ref(m); d=d(~isnan(d));
        bias(j,b)=mean(d); rmse(j,b)=sqrt(mean(d.^2));
    end
end
f4 = figure('Units','centimeters','Position',[2 2 17 7],'Color','w');
mets={bias,rmse}; ylabs2={'Bias vs \Delta\sigma=0.03 (dbar)','RMSE vs \Delta\sigma=0.03 (dbar)'};
for q=1:2
    ax=subplot(1,2,q); h=bar(ax,mets{q},'grouped');
    h(1).FaceColor=C.verm; h(2).FaceColor=C.green;
    set(ax,'XTick',1:numel(crit_use),'XTickLabel',critUseLab,'XTickLabelRotation',40);
    ylabel(ax,ylabs2{q});
    if q==1, yline(ax,0,'k-'); legend(ax,basins,'Box','off','Location','northwest','FontSize',8); end
    local_panel(ax,char('a'+q-1));
end
local_export(f4,fullfile(OUT_DIR,'Fig4_bias_rmse'));

%% ---------------------- TABLE 1 --------------------------------------
rowsT={};
for j=1:numel(crit_use), for b=1:2, for s=1:4
    m=strcmp(T.basin,basins{b})&strcmp(T.season,seasons{s});
    d=T.(crit_use{j})(m)-ref(m); d=d(~isnan(d));
    if isempty(d), bi=NaN; rm=NaN; si=NaN; n=0;
    else, bi=mean(d); rm=sqrt(mean(d.^2)); si=rm/mean(ref(m&~isnan(ref)),'omitnan'); n=numel(d); end
    rowsT(end+1,:)={crit_use{j},basins{b},seasons{s},n,round(bi,2),round(rm,2),round(si,3)}; %#ok<SAGROW>
end, end, end
Tab1=cell2table(rowsT,'VariableNames', ...
    {'Criterion','Basin','Season','N','Bias_dbar','RMSE_dbar','ScatterIndex'});
writetable(Tab1,fullfile(OUT_DIR,'Table1_error_metrics.csv'));
fprintf('Figures + Table 1 written to %s\n',OUT_DIR);

%% ---------------------- LOCAL FUNCTIONS ------------------------------
function local_panel(ax,letter)
    text(ax,0.02,0.98,['(' letter ')'],'Units','normalized','FontWeight','bold', ...
        'FontSize',11,'VerticalAlignment','top','BackgroundColor','w','Margin',0.5);
end

function local_export(fig,base,raster)
    if nargin<3, raster=false; end
    if raster, exportgraphics(fig,[base '.pdf'],'ContentType','image','Resolution',300);
    else,      exportgraphics(fig,[base '.pdf'],'ContentType','vector'); end
    exportgraphics(fig,[base '.png'],'Resolution',300);
end

function lab=local_label(c)
    if startsWith(c,'MLDdens_'), lab=sprintf('\\Delta\\sigma=%.3g',str2double(c(9:end))/1000);
    else, lab=sprintf('\\DeltaT=%.1f',str2double(c(9:end))/10); end
end

function cm=local_seq(base)
    n=64; cm=[linspace(1,base(1),n)',linspace(1,base(2),n)',linspace(1,base(3),n)'];
end

function [P,Tt,S,sig0]=local_reload_profile(rrow,argoRoot)
    d=dir(fullfile(argoRoot,'**',rrow.file{1})); fp=fullfile(d(1).folder,d(1).name);
    LAT=ncread(fp,'LATITUDE'); LON=ncread(fp,'LONGITUDE'); CYC=ncread(fp,'CYCLE_NUMBER');
    [~,ip]=min(abs(LAT-rrow.lat)+abs(LON-rrow.lon)+abs(double(CYC)-rrow.cycle));
    P=ncread(fp,'PRES_ADJUSTED'); Tt=ncread(fp,'TEMP_ADJUSTED'); S=ncread(fp,'PSAL_ADJUSTED');
    P=P(:,ip); Tt=Tt(:,ip); S=S(:,ip);
    if all(isnan(P))||all(P==99999)
        P=ncread(fp,'PRES'); Tt=ncread(fp,'TEMP'); S=ncread(fp,'PSAL'); P=P(:,ip); Tt=Tt(:,ip); S=S(:,ip);
    end
    P(P==99999)=NaN; Tt(Tt==99999)=NaN; S(S==99999)=NaN;
    ok=~isnan(P)&~isnan(Tt)&~isnan(S); P=P(ok); Tt=Tt(ok); S=S(ok);
    [P,si]=sort(P); Tt=Tt(si); S=S(si);
    if exist('gsw_SA_from_SP','file')==2
        SA=gsw_SA_from_SP(S,P,rrow.lon,rrow.lat); CT=gsw_CT_from_t(SA,Tt,P); sig0=gsw_sigma0(SA,CT);
    else, sig0=nan(size(P)); end
end
