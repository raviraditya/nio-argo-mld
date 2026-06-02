%% =======================================================================
%  make_figures_nio_mld.m
%  Journal-quality figures + Table 1 for the NIO Argo MLD threshold study.
%  Reads nio_mld_profiles.xlsx (output of process_argo_mld.m).
%
%  Figures (4 figures + 1 table):
%    Fig 1  Map of NIO: Argo profile density + AS/BoB sub-basin boxes
%    Fig 2  Example profiles (a) BoB barrier-layer  (b) AS deep-mixing,
%           showing where each criterion places the MLD
%    Fig 3  Seasonal (4-bin) sub-region panels: MLD0.03, BLT, T-bias
%    Fig 4  Bias & RMSE of each criterion by sub-basin (grouped bars)
%    Table 1  error metrics: criterion x sub-basin x season  -> CSV + TXT
%
%  NOTE: figure numbering here follows the MANUSCRIPT order. In the paper,
%  the seasonal panel is Fig. 3 and the bias/RMSE bars are Fig. 4. The
%  filenames below keep the script's build order; relabel in the paper.
%
%  STYLE: panel letters a/b/c, NO in-figure titles, colour-blind-safe
%         palette, Arial, vector/raster PDF + 300-dpi PNG.
%  EVERY plotted value is also written to a sidecar .txt for reproducibility.
%
%  Reference MLD for bias/RMSE = density MLD at 0.03 kg m-3 (de Boyer
%  Montegut 2004).
%
%  HOW TO RUN:
%    1. Run process_argo_mld.m first to create nio_mld_profiles.xlsx.
%    2. Set IN_TABLE to that file, OUT_DIR to a figure folder.
%    3. Set GSW_PATH (needed only for the Fig 2 example profiles).
%    4. Press Run (F5).
% =======================================================================

clear; clc; close all;

%% ---------------------- SETTINGS --------------------------------------
% --- EDIT THESE FOR YOUR SYSTEM ---------------------------------------
IN_TABLE  = fullfile('output','nio_mld_profiles.xlsx');  % per-profile table
OUT_DIR   = fullfile('output','figures');                % figure output folder
ARGO_ROOT = fullfile('data','argo_indian_ocean');        % needed for Fig 2 reload
GSW_PATH  = '';   % GSW toolbox folder; leave '' if already on path
% ----------------------------------------------------------------------

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
if ~isempty(GSW_PATH) && exist(GSW_PATH,'dir'), addpath(genpath(GSW_PATH)); end

% Two example profiles for Fig 2 (leave [] to auto-pick:
% BoB = largest BLT; AS = deepest MLD0.03).
EX_BOB_ROW = [];
EX_AS_ROW  = [];

% Reference criterion column for bias/RMSE
REFCOL = 'MLDdens_030';      % 0.03 kg m-3

% Colour-blind-safe palette (Okabe-Ito) -------------------------------
C.blue   = [0    114 178]/255;
C.orange = [230 159 0  ]/255;
C.green  = [0    158 115]/255;
C.verm   = [213 94  0  ]/255;
C.purple = [204 121 167]/255;
C.sky    = [86  180 233]/255;
C.yellow = [240 228 66 ]/255;
C.grey   = [120 120 120]/255;
PAL = [C.blue;C.orange;C.green;C.verm;C.purple;C.sky;C.yellow;C.grey];

set(0,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial', ...
      'DefaultAxesFontSize',9,'DefaultAxesLineWidth',0.8, ...
      'DefaultAxesBox','on','DefaultAxesLayer','top');

%% ---------------------- LOAD ------------------------------------------
T = readtable(IN_TABLE);
fprintf('Loaded %d profiles.\n', height(T));

bad = T.MLDdens_030>300 | T.ild>300 | abs(T.blt)>100;
T(bad,:) = [];
fprintf('After physical cap: %d profiles.\n', height(T));

% criterion columns
dsig_cols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'MLDdens_'));
dt_cols   = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'MLDtemp_'));
all_crit  = [dsig_cols, dt_cols];

critLabel = cellfun(@(c) local_label(c), all_crit, 'uni',0);

% season from month (NH meteorological + monsoon-aware)
mo = month(datetime(T.datestr,'InputFormat','yyyy-MM-dd'));
seas = strings(height(T),1);
seas(ismember(mo,[12 1 2]))  = "Winter(DJF)";
seas(ismember(mo,[3 4 5]))   = "PreMon(MAM)";
seas(ismember(mo,[6 7 8 9])) = "Monsoon(JJAS)";
seas(ismember(mo,[10 11]))   = "PostMon(ON)";
T.season = seas;

basins = {'AS','BoB'};
seasons = {'Winter(DJF)','PreMon(MAM)','Monsoon(JJAS)','PostMon(ON)'};

ref = T.(REFCOL);

%% =======================================================================
%  FIGURE 1  —  Map: Argo profile density + sub-basin boxes
% =======================================================================
f1 = figure('Units','centimeters','Position',[2 2 11 7],'Color','w');
ax1 = axes(f1); hold(ax1,'on');

loe = 30:0.5:100;  lae = 0:0.5:30;
N = histcounts2(T.lon, T.lat, loe, lae);
N(N==0) = NaN;
pcolor(ax1, loe(1:end-1)+0.25, lae(1:end-1)+0.25, log10(N')); shading(ax1,'flat');
cmap = local_seq(C.blue); colormap(ax1,cmap);
cb = colorbar(ax1); cb.Label.String = 'log_{10}(profiles per 0.5\circ cell)';
cb.LineWidth = 0.8;

if exist('coastlines','file')
    load coastlines coastlat coastlon
    plot(ax1, coastlon, coastlat, 'k-', 'LineWidth',0.6);
end

xline(ax1,78,'--','Color',C.grey,'LineWidth',1.2);
text(ax1, 50, 28.5, 'Arabian Sea',  'FontWeight','bold','FontSize',9,'Color',C.verm);
text(ax1, 80, 28.5, 'Bay of Bengal','FontWeight','bold','FontSize',9,'Color',C.green);

xlabel(ax1,'Longitude (\circE)'); ylabel(ax1,'Latitude (\circN)');
xlim(ax1,[30 100]); ylim(ax1,[0 30]);
daspect(ax1,[1 1 1]);
local_panel(ax1,'a');
local_export(f1, fullfile(OUT_DIR,'Fig1_map'), true);

fid=fopen(fullfile(OUT_DIR,'Fig1_map_values.txt'),'w');
fprintf(fid,'Fig1: Argo profile density (0.5-deg bins).\n');
fprintf(fid,'lon_center\tlat_center\tn_profiles\n');
for i=1:numel(loe)-1, for j=1:numel(lae)-1
    if ~isnan(N(i,j)), fprintf(fid,'%.2f\t%.2f\t%d\n',loe(i)+0.25,lae(j)+0.25,N(i,j)); end
end, end
fclose(fid);

%% =======================================================================
%  FIGURE 2  —  Example profiles (intuition)
% =======================================================================
if isempty(EX_BOB_ROW)
    isB = strcmp(T.basin,'BoB') & ~isnan(T.blt);
    [~,ix] = max(T.blt .* isB); EX_BOB_ROW = ix;
end
if isempty(EX_AS_ROW)
    isA = strcmp(T.basin,'AS') & ~isnan(ref);
    tmp = ref; tmp(~isA)=NaN; [~,ix]=max(tmp); EX_AS_ROW = ix;
end
f2 = figure('Units','centimeters','Position',[2 2 16 8],'Color','w');
for pp = 1:2
    if pp==1, rrow=EX_BOB_ROW; lab='a'; ttl='Bay of Bengal'; cc=C.green;
    else,     rrow=EX_AS_ROW;  lab='b'; ttl='Arabian Sea';  cc=C.verm; end
    [P,Tt,~,sig0] = local_reload_profile(T(rrow,:), ARGO_ROOT);
    axL = subplot(1,2,pp); hold(axL,'on');
    plot(axL, Tt, P, '-','Color',cc,'LineWidth',1.6);
    set(axL,'YDir','reverse'); ylim(axL,[0 200]);
    xlabel(axL,'Temperature (\circC)'); if pp==1, ylabel(axL,'Pressure (dbar)'); end
    axT = axes('Position',axL.Position,'Color','none', ...
               'XAxisLocation','top','YAxisLocation','right','YDir','reverse');
    hold(axT,'on'); ylim(axT,[0 200]);
    plot(axT, sig0, P, '-','Color',C.grey,'LineWidth',1.2);
    axT.YTickLabel = []; xlabel(axT,'\sigma_\theta (kg m^{-3})','Color',C.grey);
    crit_vals = [];
    for j=1:numel(all_crit), crit_vals(j) = T.(all_crit{j})(rrow); end %#ok<SAGROW>
    yl = crit_vals;
    for j=1:numel(all_crit)
        yline(axL, yl(j), ':','Color',PAL(mod(j-1,size(PAL,1))+1,:),'LineWidth',1.0);
    end
    yline(axL, T.ild(rrow),'-.','Color',C.orange,'LineWidth',1.3);
    text(axL, min(Tt)+0.3, 8, ttl,'FontWeight','bold','Color',cc,'FontSize',9);
    local_panel(axL,lab);
end
lh = legend(axL, [critLabel,{'ILD (\DeltaT=0.5)'}],'Location','eastoutside', ...
            'FontSize',7,'Box','off'); %#ok<NASGU>
local_export(f2, fullfile(OUT_DIR,'Fig2_profiles'));

fid=fopen(fullfile(OUT_DIR,'Fig2_profiles_values.txt'),'w');
fprintf(fid,'Fig2 example profiles. Rows used: BoB=%d AS=%d\n',EX_BOB_ROW,EX_AS_ROW);
for rrow=[EX_BOB_ROW EX_AS_ROW]
    fprintf(fid,'--- row %d  basin=%s  date=%s  lat=%.2f lon=%.2f ---\n', ...
        rrow, T.basin{rrow}, string(T.datestr(rrow)), T.lat(rrow), T.lon(rrow));
    for j=1:numel(all_crit), fprintf(fid,'%s = %.1f dbar\n', all_crit{j}, T.(all_crit{j})(rrow)); end
    fprintf(fid,'ILD = %.1f  BLT = %.1f\n', T.ild(rrow), T.blt(rrow));
end
fclose(fid);

%% =======================================================================
%  FIGURE 3 (paper)  —  Seasonal (4-bin) sub-region panels
%  (a) mean MLD0.03  (b) mean BLT  (c) DeltaT=0.2 bias, grouped by sub-region
% =======================================================================
regions5 = {'AS-N','AS-S','BoB-N','BoB-S'};
seas5    = seasons;
regPAL   = [C.verm; C.orange; C.green; C.sky];

mldM  = nan(numel(seas5), numel(regions5));
bltM  = nan(numel(seas5), numel(regions5));
biasM = nan(numel(seas5), numel(regions5));
for s = 1:numel(seas5)
    for r = 1:numel(regions5)
        m = strcmp(T.region,regions5{r}) & strcmp(T.season,seas5{s});
        x = ref(m); x = x(~isnan(x));
        if ~isempty(x), mldM(s,r) = mean(x); end
        b = T.blt(m); b = b(~isnan(b));
        if ~isempty(b), bltM(s,r) = mean(b); end
        d = T.MLDtemp_02(m) - ref(m); d = d(~isnan(d));
        if ~isempty(d), biasM(s,r) = mean(d); end
    end
end
seasTick = {'Win','Pre','Mon','Post'};

f3 = figure('Units','centimeters','Position',[2 2 17 16],'Color','w');
axa = subplot(3,1,1);
ha = bar(axa, mldM, 'grouped');
for r=1:numel(regions5), ha(r).FaceColor = regPAL(r,:); end
set(axa,'XTick',1:numel(seas5),'XTickLabel',seasTick);
ylabel(axa,'Mean MLD_{0.03} (dbar)');
legend(axa, regions5, 'Box','off','Location','northwest','FontSize',7,'NumColumns',2);
local_panel(axa,'a');

axb = subplot(3,1,2);
hb = bar(axb, bltM, 'grouped');
for r=1:numel(regions5), hb(r).FaceColor = regPAL(r,:); end
set(axb,'XTick',1:numel(seas5),'XTickLabel',seasTick);
ylabel(axb,'Mean BLT (dbar)'); yline(axb,0,'k-');
local_panel(axb,'b');

axc = subplot(3,1,3);
hc = bar(axc, biasM, 'grouped');
for r=1:numel(regions5), hc(r).FaceColor = regPAL(r,:); end
set(axc,'XTick',1:numel(seas5),'XTickLabel',seasTick);
ylabel(axc,'\DeltaT=0.2 bias (dbar)'); yline(axc,0,'k-');
xlabel(axc,'Season');
local_panel(axc,'c');
local_export(f3, fullfile(OUT_DIR,'Fig3_seasonal'));

fid=fopen(fullfile(OUT_DIR,'Fig3_seasonal_values.txt'),'w');
fprintf(fid,'Fig3 (paper): seasonal (4-bin) means by sub-basin.\n');
fprintf(fid,'panel\tseason\tregion\tvalue_dbar\n');
labs = {'a_MLD0.03','b_BLT','c_biasT0.2'};
mats = {mldM, bltM, biasM};
for q=1:3
    for s=1:numel(seas5), for r=1:numel(regions5)
        fprintf(fid,'%s\t%s\t%s\t%.2f\n',labs{q},seas5{s},regions5{r},mats{q}(s,r));
    end, end
end
fclose(fid);

%% =======================================================================
%  FIGURE 4 (paper)  —  Bias & RMSE by sub-basin
% =======================================================================
crit_use = setdiff(all_crit, REFCOL, 'stable');
critUseLab = cellfun(@(c) local_label(c), crit_use, 'uni',0);

bias = nan(numel(crit_use), numel(basins));
rmse = nan(numel(crit_use), numel(basins));
for b=1:numel(basins)
    m = strcmp(T.basin,basins{b});
    for j=1:numel(crit_use)
        d = T.(crit_use{j})(m) - ref(m);
        d = d(~isnan(d));
        bias(j,b) = mean(d);
        rmse(j,b) = sqrt(mean(d.^2));
    end
end

f4 = figure('Units','centimeters','Position',[2 2 17 7],'Color','w');
axa = subplot(1,2,1);
hbar = bar(axa, bias,'grouped'); hbar(1).FaceColor=C.verm; hbar(2).FaceColor=C.green;
set(axa,'XTick',1:numel(crit_use),'XTickLabel',critUseLab,'XTickLabelRotation',40);
ylabel(axa,'Bias vs \Delta\sigma=0.03 (dbar)'); yline(axa,0,'k-');
legend(axa,basins,'Box','off','Location','northwest','FontSize',8);
local_panel(axa,'a');
axb = subplot(1,2,2);
hr = bar(axb, rmse,'grouped'); hr(1).FaceColor=C.verm; hr(2).FaceColor=C.green;
set(axb,'XTick',1:numel(crit_use),'XTickLabel',critUseLab,'XTickLabelRotation',40);
ylabel(axb,'RMSE vs \Delta\sigma=0.03 (dbar)');
local_panel(axb,'b');
local_export(f4, fullfile(OUT_DIR,'Fig4_bias_rmse'));

fid=fopen(fullfile(OUT_DIR,'Fig4_bias_rmse_values.txt'),'w');
fprintf(fid,'Fig4 (paper): bias & RMSE vs %s, by sub-basin.\n',REFCOL);
fprintf(fid,'criterion\tbasin\tbias_dbar\trmse_dbar\n');
for j=1:numel(crit_use), for b=1:numel(basins)
    fprintf(fid,'%s\t%s\t%.2f\t%.2f\n',crit_use{j},basins{b},bias(j,b),rmse(j,b));
end, end
fclose(fid);

%% =======================================================================
%  TABLE 1  —  error metrics: criterion x sub-basin x season
% =======================================================================
rowsT = {};
for j=1:numel(crit_use)
  for b=1:numel(basins)
    for s=1:numel(seasons)
      m = strcmp(T.basin,basins{b}) & strcmp(T.season,seasons{s});
      d = T.(crit_use{j})(m) - ref(m); d=d(~isnan(d));
      if isempty(d), bi=NaN; rm=NaN; si=NaN; n=0;
      else
          bi=mean(d); rm=sqrt(mean(d.^2));
          si=rm/mean(ref(m & ~isnan(ref)),'omitnan');
          n=numel(d);
      end
      rowsT(end+1,:) = {crit_use{j}, basins{b}, seasons{s}, n, ...
                        round(bi,2), round(rm,2), round(si,3)}; %#ok<SAGROW>
    end
  end
end
Tab1 = cell2table(rowsT,'VariableNames', ...
    {'Criterion','Basin','Season','N','Bias_dbar','RMSE_dbar','ScatterIndex'});
writetable(Tab1, fullfile(OUT_DIR,'Table1_error_metrics.csv'));
fid=fopen(fullfile(OUT_DIR,'Table1_error_metrics.txt'),'w');
fprintf(fid,'Table 1. MLD-criterion error metrics vs %s.\n',REFCOL);
fprintf(fid,'%-14s %-5s %-14s %6s %10s %10s %12s\n', ...
    'Criterion','Basin','Season','N','Bias','RMSE','ScatterIdx');
for i=1:height(Tab1)
    fprintf(fid,'%-14s %-5s %-14s %6d %10.2f %10.2f %12.3f\n', ...
        Tab1.Criterion{i},Tab1.Basin{i},Tab1.Season{i},Tab1.N(i), ...
        Tab1.Bias_dbar(i),Tab1.RMSE_dbar(i),Tab1.ScatterIndex(i));
end
fclose(fid);

fprintf('\nAll figures + Table 1 written to:\n  %s\n', OUT_DIR);

%% ---------------------- LOCAL FUNCTIONS -------------------------------
function local_panel(ax,letter)
    text(ax,0.02,0.98,['(' letter ')'],'Units','normalized', ...
        'FontWeight','bold','FontSize',11,'VerticalAlignment','top', ...
        'BackgroundColor','w','Margin',0.5);
end

function local_export(fig, base, raster)
    if nargin<3, raster=false; end
    if raster
        exportgraphics(fig,[base '.pdf'],'ContentType','image','Resolution',300);
    else
        exportgraphics(fig,[base '.pdf'],'ContentType','vector');
    end
    exportgraphics(fig,[base '.png'],'Resolution',300);
end

function lab = local_label(c)
    if startsWith(c,'MLDdens_')
        v = str2double(c(9:end))/1000;
        lab = sprintf('\\Delta\\sigma=%.3g', v);
    else
        v = str2double(c(9:end))/10;
        lab = sprintf('\\DeltaT=%.1f', v);
    end
end

function cm = local_seq(base)
    n=64; w=[1 1 1];
    cm = [linspace(w(1),base(1),n)', linspace(w(2),base(2),n)', linspace(w(3),base(3),n)'];
end

function [P,Tt,S,sig0] = local_reload_profile(rrow, argoRoot)
    % Re-open the original NetCDF file for one example profile to get the
    % full vertical T/S/P (not stored in the summary table).
    d = dir(fullfile(argoRoot,'**',rrow.file{1}));
    fp = fullfile(d(1).folder, d(1).name);
    LAT=ncread(fp,'LATITUDE'); LON=ncread(fp,'LONGITUDE'); CYC=ncread(fp,'CYCLE_NUMBER');
    [~,ip] = min(abs(LAT-rrow.lat)+abs(LON-rrow.lon)+abs(double(CYC)-rrow.cycle));
    P=ncread(fp,'PRES_ADJUSTED'); Tt=ncread(fp,'TEMP_ADJUSTED'); S=ncread(fp,'PSAL_ADJUSTED');
    P=P(:,ip);Tt=Tt(:,ip);S=S(:,ip);
    if all(isnan(P))||all(P==99999)
        P=ncread(fp,'PRES');Tt=ncread(fp,'TEMP');S=ncread(fp,'PSAL');P=P(:,ip);Tt=Tt(:,ip);S=S(:,ip);
    end
    P(P==99999)=NaN;Tt(Tt==99999)=NaN;S(S==99999)=NaN;
    ok=~isnan(P)&~isnan(Tt)&~isnan(S); P=P(ok);Tt=Tt(ok);S=S(ok);
    [P,si]=sort(P);Tt=Tt(si);S=S(si);
    if exist('gsw_SA_from_SP','file')==2
        SA=gsw_SA_from_SP(S,P,rrow.lon,rrow.lat); CT=gsw_CT_from_t(SA,Tt,P);
        sig0=gsw_sigma0(SA,CT);
    else
        sig0 = nan(size(P));
    end
end