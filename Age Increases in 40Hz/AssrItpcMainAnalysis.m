
%% Main ASSR_ITPC_Analysis
%  Full statistical analysis pipeline for ASSR 40Hz ITPC manuscript
%  Shlomit Beker shlomitbeker@gmail.com
%  Models: 0A (KW: TD vs ASD vs PMS, collapsed),
%          0B (OLS: all groups all ages),
%          1 (TD vs ASD w/o ID), 2 (TD + ASD w/ID + PMS under 18),
%          3S (IQ covariate), Sex effects post-hoc,
%          Trial count analysis
%
%  Data file: DATA_SB.xlsx
%  Columns: DX, Age, Gender, ITPC, Segments, ID, IQ
%
%  EXCLUSIONS:
%  - PMS adults (n=2) excluded
%  - ASD w/ID adults (n=1) excluded
%  - Total analytic N=128

clear; clc;

%% ============================================================
%  PLOT SETTINGS
%% ============================================================

col_TD  = [0.13 0.47 0.71];   % blue
col_ASD = [0.84 0.15 0.16];   % red
col_PMS = [0.17 0.63 0.17];   % green
col_M   = [0.12 0.47 0.71];
col_F   = [0.84 0.15 0.16];
jitter  = 0.12;

age_range_fine = linspace(2, 38, 300);


%% ============================================================
%  LOAD AND PREPARE DATA
%% ============================================================

T = readtable('DATA_SB.xlsx');

% Convert IQ to numeric
T.IQ = str2double(string(T.IQ));

% Exclude PMS adults and ASD w/ID adults
T = T(~(strcmp(T.DX,'PMS') & T.Age >= 18), :);
T = T(~(strcmp(T.DX,'ASD') & (strcmp(T.ID,'ID') | strcmp(T.ID,'PMSID')) & T.Age >= 18), :);

% Log-transform age
T.logAge = log(T.Age);

% Gender: M=0, F=1
T.Sex = double(strcmp(T.Gender, 'F'));

% Group flags
is_TD      = strcmp(T.DX, 'TD');
is_ASD     = strcmp(T.DX, 'ASD');
is_PMS     = strcmp(T.DX, 'PMS');
is_noID    = strcmp(T.ID, 'noID');
is_ID      = strcmp(T.ID, 'ID') | strcmp(T.ID, 'PMSID');
is_adult   = T.Age >= 18;
is_under18 = T.Age < 18;

fprintf('=== Data loaded: N=%d ===\n', height(T));
fprintf('TD=%d | ASD=%d | PMS=%d\n', sum(is_TD), sum(is_ASD), sum(is_PMS));


%% ============================================================
%  DESCRIPTIVE STATISTICS
%% ============================================================

fprintf('\n=== Descriptive Statistics ===\n');
grp_list = {
    'TD (all)',              T(is_TD, :);
    'TD Under 18',           T(is_TD & is_under18, :);
    'TD Adults',             T(is_TD & is_adult, :);
    'ASD w/o ID Under 18',  T(is_ASD & is_noID & is_under18, :);
    'ASD w/o ID Adults',    T(is_ASD & is_noID & is_adult, :);
    'ASD w/ID Under 18',    T(is_ASD & is_ID & is_under18, :);
    'PMS Under 18',         T(is_PMS & is_under18, :);
};

for g = 1:size(grp_list,1)
    name = grp_list{g,1};
    grp  = grp_list{g,2};
    fprintf('%s: n=%d, Age M=%.1f (SD=%.1f), ITPC M=%.3f (SD=%.3f)\n', ...
        name, height(grp), mean(grp.Age), std(grp.Age), ...
        mean(grp.ITPC), std(grp.ITPC));
end


%% ============================================================
%  NORMALITY CHECK
%% ============================================================

fprintf('\n=== Normality (Lilliefors) ===\n');
[h, p] = lillietest(T.ITPC);
fprintf('h=%d, p=%.4f\n', h, p);


%% ============================================================
%  MODEL 0A: KW — TD vs ASD (all) vs PMS (ages collapsed)
%% ============================================================

fprintf('\n=== MODEL 0A: KW — TD vs ASD vs PMS (ages collapsed) ===\n');

td_itpc  = T.ITPC(is_TD);
asd_itpc = T.ITPC(is_ASD);
pms_itpc = T.ITPC(is_PMS);

[p_kw0a, tbl_kw0a] = kruskalwallis([td_itpc; asd_itpc; pms_itpc], ...
    [repmat({'TD'},sum(is_TD),1); repmat({'ASD'},sum(is_ASD),1); repmat({'PMS'},sum(is_PMS),1)], 'off');
fprintf('KW H=%.3f, p=%.4f\n', tbl_kw0a{2,5}, p_kw0a);
fprintf('TD:  n=%d, M=%.3f (SD=%.3f)\n', sum(is_TD),  mean(td_itpc),  std(td_itpc));
fprintf('ASD: n=%d, M=%.3f (SD=%.3f)\n', sum(is_ASD), mean(asd_itpc), std(asd_itpc));
fprintf('PMS: n=%d, M=%.3f (SD=%.3f)\n', sum(is_PMS), mean(pms_itpc), std(pms_itpc));

% Pairwise (Bonferroni k=3)
fprintf('\nPairwise Mann-Whitney U (Bonferroni k=3):\n');
pairs = {td_itpc, asd_itpc, 'TD','ASD';
         td_itpc, pms_itpc, 'TD','PMS';
         asd_itpc,pms_itpc, 'ASD','PMS'};
for k = 1:3
    [p_mw,~] = ranksum(pairs{k,1}, pairs{k,2});
    fprintf('  %s vs %s: p=%.4f, p_bonf=%.4f\n', pairs{k,3}, pairs{k,4}, p_mw, min(p_mw*3,1));
end


%% ---- MODEL 0A PLOT ----

figure('Name','Model 0A: All Groups Collapsed','Color','w','Position',[50 50 500 480]);
hold on;

data_0a   = [td_itpc; asd_itpc; pms_itpc];
labels_0a = [repmat({'TD'},length(td_itpc),1); ...
             repmat({'ASD'},length(asd_itpc),1); ...
             repmat({'PMS'},length(pms_itpc),1)];
bp = boxplot(data_0a, labels_0a, 'Colors',[col_TD;col_ASD;col_PMS], ...
    'Symbol','o','Widths',0.5);
set(bp,'LineWidth',1.5);

scatter(1+(rand(length(td_itpc),1)-0.5)*jitter,  td_itpc,  30,col_TD,  'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(length(asd_itpc),1)-0.5)*jitter, asd_itpc, 30,col_ASD, 'filled','MarkerFaceAlpha',0.5);
scatter(3+(rand(length(pms_itpc),1)-0.5)*jitter, pms_itpc, 30,col_PMS, 'filled','MarkerFaceAlpha',0.5);

ylabel('40Hz ITPC','FontSize',12);
title(sprintf('All Groups — Ages Collapsed\n(KW H=%.3f, p=%.3f)', tbl_kw0a{2,5}, p_kw0a), ...
    'FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',10);
ylim([0 0.7]); hold off;


%% ============================================================
%  MODEL 0B: OLS — ITPC ~ logAge * DX + Sex (all groups, all ages)
%% ============================================================

fprintf('\n=== MODEL 0B: OLS — ITPC ~ logAge * DX + Sex (all groups) ===\n');

% DX: TD=reference
T.DX_3 = categorical(T.DX, {'TD','ASD','PMS'});
mdl0b = fitlm(T, 'ITPC ~ logAge * DX_3 + Sex');
disp(mdl0b);
fprintf('N=%d, R2=%.3f\n', height(T), mdl0b.Rsquared.Ordinary);

% Within-group slopes
fprintf('\nWithin-group slopes (with Sex):\n');
T0b_TD  = T(is_TD, :);
T0b_ASD = T(is_ASD, :);
T0b_PMS = T(is_PMS, :);

mdl0b_TD  = fitlm(T0b_TD,  'ITPC ~ logAge + Sex');
mdl0b_ASD = fitlm(T0b_ASD, 'ITPC ~ logAge + Sex');
mdl0b_PMS = fitlm(T0b_PMS, 'ITPC ~ logAge + Sex');

fprintf('TD:  beta=%.3f, SE=%.3f, p=%.4f\n', ...
    mdl0b_TD.Coefficients.Estimate(2),  mdl0b_TD.Coefficients.SE(2),  mdl0b_TD.Coefficients.pValue(2));
fprintf('ASD: beta=%.3f, SE=%.3f, p=%.4f\n', ...
    mdl0b_ASD.Coefficients.Estimate(2), mdl0b_ASD.Coefficients.SE(2), mdl0b_ASD.Coefficients.pValue(2));
fprintf('PMS: beta=%.3f, SE=%.3f, p=%.4f\n', ...
    mdl0b_PMS.Coefficients.Estimate(2), mdl0b_PMS.Coefficients.SE(2), mdl0b_PMS.Coefficients.pValue(2));

% KW by age group
T0b_adults   = T(is_adult, :);
T0b_children = T(is_under18, :);

fprintf('\nKW adults (TD vs ASD only — no PMS adults):\n');
[p_kw0b_ad, tbl_kw0b_ad] = kruskalwallis( ...
    [T0b_adults.ITPC(strcmp(T0b_adults.DX,'TD')); T0b_adults.ITPC(strcmp(T0b_adults.DX,'ASD'))], ...
    [repmat({'TD'},sum(strcmp(T0b_adults.DX,'TD')),1); ...
     repmat({'ASD'},sum(strcmp(T0b_adults.DX,'ASD')),1)], 'off');
fprintf('  H=%.3f, p=%.4f\n', tbl_kw0b_ad{2,5}, p_kw0b_ad);

fprintf('KW children (TD vs ASD vs PMS):\n');
[p_kw0b_ch, tbl_kw0b_ch] = kruskalwallis( ...
    [T0b_children.ITPC(strcmp(T0b_children.DX,'TD')); ...
     T0b_children.ITPC(strcmp(T0b_children.DX,'ASD')); ...
     T0b_children.ITPC(strcmp(T0b_children.DX,'PMS'))], ...
    [repmat({'TD'},sum(strcmp(T0b_children.DX,'TD')),1); ...
     repmat({'ASD'},sum(strcmp(T0b_children.DX,'ASD')),1); ...
     repmat({'PMS'},sum(strcmp(T0b_children.DX,'PMS')),1)], 'off');
fprintf('  H=%.3f, p=%.4f\n', tbl_kw0b_ch{2,5}, p_kw0b_ch);


%% ---- MODEL 0B PLOTS ----

figure('Name','Model 0B: All Groups All Ages','Color','w','Position',[50 50 1100 480]);

% Panel 1: Scatter + regression lines
subplot(1,2,1); hold on;

scatter(T0b_TD.Age,  T0b_TD.ITPC,  40,col_TD,  'o','filled','MarkerFaceAlpha',0.5);
scatter(T0b_ASD.Age, T0b_ASD.ITPC, 40,col_ASD, 's','filled','MarkerFaceAlpha',0.5);
scatter(T0b_PMS.Age, T0b_PMS.ITPC, 40,col_PMS, '^','filled','MarkerFaceAlpha',0.5);

mean_sex0 = mean(T.Sex);
b_TD0  = mdl0b_TD.Coefficients.Estimate;
b_ASD0 = mdl0b_ASD.Coefficients.Estimate;
b_PMS0 = mdl0b_PMS.Coefficients.Estimate;

plot(age_range_fine, b_TD0(1)  + b_TD0(2)*log(age_range_fine)  + b_TD0(3)*mean_sex0,  '-',  'Color',col_TD,  'LineWidth',2.5);
plot(age_range_fine, b_ASD0(1) + b_ASD0(2)*log(age_range_fine) + b_ASD0(3)*mean_sex0, '--', 'Color',col_ASD, 'LineWidth',2.5);
plot(age_range_fine, b_PMS0(1) + b_PMS0(2)*log(age_range_fine) + b_PMS0(3)*mean_sex0, ':',  'Color',col_PMS, 'LineWidth',2.5);

xlabel('Age (years)','FontSize',11); ylabel('40Hz ITPC','FontSize',11);
title('ITPC vs Age — All Groups','FontSize',11,'FontWeight','bold');
legend({'TD','ASD','PMS'},'Location','northwest','FontSize',9,'Box','off');
set(gca,'Box','off','TickDir','out','FontSize',10);
ylim([0 0.65]); xlim([1 40]); hold off;

% Panel 2: Boxplot split by adults vs children
subplot(1,2,2); hold on;

% Adults: TD vs ASD only
TD_ad  = T0b_adults(strcmp(T0b_adults.DX,'TD'), :);
ASD_ad = T0b_adults(strcmp(T0b_adults.DX,'ASD'), :);
% Children: all 3
TD_ch  = T0b_children(strcmp(T0b_children.DX,'TD'), :);
ASD_ch = T0b_children(strcmp(T0b_children.DX,'ASD'), :);
PMS_ch = T0b_children(strcmp(T0b_children.DX,'PMS'), :);

data_split   = [TD_ch.ITPC; ASD_ch.ITPC; PMS_ch.ITPC; TD_ad.ITPC; ASD_ad.ITPC];
labels_split = [repmat({'TD U18'},   height(TD_ch), 1); ...
                repmat({'ASD U18'},  height(ASD_ch), 1); ...
                repmat({'PMS U18'},  height(PMS_ch), 1); ...
                repmat({'TD Adult'}, height(TD_ad),  1); ...
                repmat({'ASD Adult'},height(ASD_ad), 1)];
col_split = [col_TD*0.7+[0.3 0.3 0.3]*0.3; col_ASD*0.7+[0.3 0.3 0.3]*0.3; col_PMS*0.7+[0.3 0.3 0.3]*0.3; col_TD; col_ASD];

bp2 = boxplot(data_split, labels_split, 'Symbol','o','Widths',0.5);
set(bp2,'LineWidth',1.5);

scatter(1+(rand(height(TD_ch),1)-0.5)*jitter,  TD_ch.ITPC,  20,col_TD*0.7+[0.3 0.3 0.3]*0.3,'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(height(ASD_ch),1)-0.5)*jitter, ASD_ch.ITPC, 20,col_ASD*0.7+[0.3 0.3 0.3]*0.3,'filled','MarkerFaceAlpha',0.5);
scatter(3+(rand(height(PMS_ch),1)-0.5)*jitter, PMS_ch.ITPC, 20,col_PMS*0.7+[0.3 0.3 0.3]*0.3,'filled','MarkerFaceAlpha',0.5);
scatter(4+(rand(height(TD_ad),1)-0.5)*jitter,  TD_ad.ITPC,  20,col_TD, 'filled','MarkerFaceAlpha',0.5);
scatter(5+(rand(height(ASD_ad),1)-0.5)*jitter, ASD_ad.ITPC, 20,col_ASD,'filled','MarkerFaceAlpha',0.5);

% Significance bar for adults
if p_kw0b_ad < 0.05
    y_sig = max([TD_ad.ITPC; ASD_ad.ITPC])+0.03;
    plot([4 5],[y_sig y_sig],'k-','LineWidth',1.2);
    text(4.5,y_sig+0.01,'*','HorizontalAlignment','center','FontSize',14);
end

ylabel('40Hz ITPC','FontSize',11);
title('ITPC by Group and Age','FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',9,'XTickLabelRotation',15);
ylim([0 0.7]); hold off;

sgtitle('Model 0B: All Groups All Ages','FontSize',13,'FontWeight','bold');


%% ============================================================
%  MODEL 1: TD vs ASD w/o ID (all ages)
%  ITPC ~ log(Age) * DX + Sex
%% ============================================================

fprintf('\n=== MODEL 1: TD vs ASD w/o ID (all ages) ===\n');

idx1      = (is_TD | (is_ASD & is_noID));
T1        = T(idx1, :);
T1.DX_bin = categorical(strcmp(T1.DX,'ASD'), [0 1], {'TD','ASD'});

mdl1 = fitlm(T1, 'ITPC ~ logAge * DX_bin + Sex');
disp(mdl1);
fprintf('R2=%.3f\n', mdl1.Rsquared.Ordinary);

T1_TD  = T1(strcmp(T1.DX,'TD'), :);
T1_ASD = T1(strcmp(T1.DX,'ASD'), :);
mdl1_TD  = fitlm(T1_TD,  'ITPC ~ logAge + Sex');
mdl1_ASD = fitlm(T1_ASD, 'ITPC ~ logAge + Sex');

fprintf('\nWithin-group TD:  beta=%.3f, SE=%.3f, p=%.4f\n', ...
    mdl1_TD.Coefficients.Estimate(2),  mdl1_TD.Coefficients.SE(2),  mdl1_TD.Coefficients.pValue(2));
fprintf('Within-group ASD: beta=%.3f, SE=%.3f, p=%.4f\n', ...
    mdl1_ASD.Coefficients.Estimate(2), mdl1_ASD.Coefficients.SE(2), mdl1_ASD.Coefficients.pValue(2));

T1_adults   = T1(T1.Age >= 18, :);
T1_children = T1(T1.Age <  18, :);

[p_kw1, tbl_kw1] = kruskalwallis(T1.ITPC, T1.DX_bin, 'off');
fprintf('\nKW all ages: H=%.3f, p=%.4f\n', tbl_kw1{2,5}, p_kw1);

[p_kw1a, tbl_kw1a] = kruskalwallis(T1_adults.ITPC,   T1_adults.DX_bin,   'off');
fprintf('KW adults:   H=%.3f, p=%.4f\n', tbl_kw1a{2,5}, p_kw1a);

[p_kw1c, tbl_kw1c] = kruskalwallis(T1_children.ITPC, T1_children.DX_bin, 'off');
fprintf('KW children: H=%.3f, p=%.4f\n', tbl_kw1c{2,5}, p_kw1c);


%% ---- MODEL 1 PLOTS ----

figure('Name','Model 1','Color','w','Position',[50 50 1200 480]);

subplot(1,3,1); hold on;
scatter(T1_TD.Age,  T1_TD.ITPC,  40,col_TD,  'o','filled','MarkerFaceAlpha',0.5);
scatter(T1_ASD.Age, T1_ASD.ITPC, 40,col_ASD, 's','filled','MarkerFaceAlpha',0.5);

mean_sex1 = mean(T1.Sex);
b_TD1  = mdl1_TD.Coefficients.Estimate;
b_ASD1 = mdl1_ASD.Coefficients.Estimate;
plot(age_range_fine, b_TD1(1)+b_TD1(2)*log(age_range_fine)+b_TD1(3)*mean_sex1,   '-',  'Color',col_TD,  'LineWidth',2.5);
plot(age_range_fine, b_ASD1(1)+b_ASD1(2)*log(age_range_fine)+b_ASD1(3)*mean_sex1,'--', 'Color',col_ASD, 'LineWidth',2.5);

xlabel('Age (years)','FontSize',11); ylabel('40Hz ITPC','FontSize',11);
title('ITPC vs Age (all ages)','FontSize',11,'FontWeight','bold');
legend({'TD','ASD w/o ID'},'Location','northwest','FontSize',9,'Box','off');
set(gca,'Box','off','TickDir','out','FontSize',10);
ylim([0 0.65]); xlim([1 40]); hold off;

subplot(1,3,2); hold on;
TD_ad1  = T1_adults(strcmp(T1_adults.DX,'TD'), :);
ASD_ad1 = T1_adults(strcmp(T1_adults.DX,'ASD'), :);
data_ad1   = [TD_ad1.ITPC; ASD_ad1.ITPC];
labels_ad1 = [repmat({'TD'},height(TD_ad1),1); repmat({'ASD'},height(ASD_ad1),1)];
bp3 = boxplot(data_ad1,labels_ad1,'Colors',[col_TD;col_ASD],'Symbol','o','Widths',0.5);
set(bp3,'LineWidth',1.5);
scatter(1+(rand(height(TD_ad1),1)-0.5)*jitter,  TD_ad1.ITPC,  25,col_TD,  'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(height(ASD_ad1),1)-0.5)*jitter, ASD_ad1.ITPC, 25,col_ASD, 'filled','MarkerFaceAlpha',0.5);
if p_kw1a < 0.05
    y_sig = max(data_ad1)+0.03;
    plot([1 2],[y_sig y_sig],'k-','LineWidth',1.2);
    text(1.5,y_sig+0.01,'*','HorizontalAlignment','center','FontSize',14);
end
ylabel('40Hz ITPC','FontSize',11);
title(sprintf('Adults only\n(KW p=%.3f)',p_kw1a),'FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',10); ylim([0 0.7]); hold off;

subplot(1,3,3); hold on;
TD_ch1  = T1_children(strcmp(T1_children.DX,'TD'), :);
ASD_ch1 = T1_children(strcmp(T1_children.DX,'ASD'), :);
data_ch1   = [TD_ch1.ITPC; ASD_ch1.ITPC];
labels_ch1 = [repmat({'TD'},height(TD_ch1),1); repmat({'ASD'},height(ASD_ch1),1)];
bp4 = boxplot(data_ch1,labels_ch1,'Colors',[col_TD;col_ASD],'Symbol','o','Widths',0.5);
set(bp4,'LineWidth',1.5);
scatter(1+(rand(height(TD_ch1),1)-0.5)*jitter,  TD_ch1.ITPC,  25,col_TD,  'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(height(ASD_ch1),1)-0.5)*jitter, ASD_ch1.ITPC, 25,col_ASD, 'filled','MarkerFaceAlpha',0.5);
ylabel('40Hz ITPC','FontSize',11);
title(sprintf('Children only\n(KW p=%.3f)',p_kw1c),'FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',10); ylim([0 0.7]); hold off;

sgtitle('Model 1: TD vs ASD w/o ID','FontSize',13,'FontWeight','bold');


%% ============================================================
%  MODEL 2: TD + ASD w/ID + PMS (under 18 only)
%% ============================================================

fprintf('\n=== MODEL 2: TD + ASD w/ID + PMS (under 18) ===\n');

idx2 = is_under18 & (is_TD | (is_ASD & is_ID) | is_PMS);
T2   = T(idx2, :);
T2.DX3 = categorical(T2.DX, {'TD','ASD','PMS'});

mdl2 = fitlm(T2, 'ITPC ~ logAge * DX3 + Sex');
disp(mdl2);
fprintf('R2=%.3f\n', mdl2.Rsquared.Ordinary);

for grpname = {'TD','ASD','PMS'}
    sub = T2(strcmp(string(T2.DX), grpname{1}), :);
    if height(sub) > 3
        m = fitlm(sub, 'ITPC ~ logAge + Sex');
        fprintf('Within-group %s: beta=%.3f, SE=%.3f, p=%.4f, n=%d\n', ...
            grpname{1}, m.Coefficients.Estimate(2), m.Coefficients.SE(2), ...
            m.Coefficients.pValue(2), height(sub));
    end
end

[p_kw2, tbl_kw2] = kruskalwallis(T2.ITPC, T2.DX3, 'off');
fprintf('\nKW Model 2: H=%.3f, p=%.4f\n', tbl_kw2{2,5}, p_kw2);


%% ---- MODEL 2 PLOTS ----

figure('Name','Model 2','Color','w','Position',[50 50 1100 480]);

subplot(1,2,1); hold on;
T2_TD  = T2(strcmp(T2.DX,'TD'), :);
T2_ASD = T2(strcmp(T2.DX,'ASD'), :);
T2_PMS = T2(strcmp(T2.DX,'PMS'), :);
scatter(T2_TD.Age,  T2_TD.ITPC,  40,col_TD,  'o','filled','MarkerFaceAlpha',0.5);
scatter(T2_ASD.Age, T2_ASD.ITPC, 40,col_ASD, 's','filled','MarkerFaceAlpha',0.5);
scatter(T2_PMS.Age, T2_PMS.ITPC, 40,col_PMS, '^','filled','MarkerFaceAlpha',0.5);

age_u18   = linspace(1.5,18,200);
mean_sex2 = mean(T2.Sex);
lsmap = {'TD',col_TD,'-'; 'ASD',col_ASD,'--'; 'PMS',col_PMS,':'};
for k = 1:3
    sub = T2(strcmp(string(T2.DX), lsmap{k,1}), :);
    if height(sub) > 3
        m = fitlm(sub, 'ITPC ~ logAge + Sex');
        b = m.Coefficients.Estimate;
        plot(age_u18, b(1)+b(2)*log(age_u18)+b(3)*mean_sex2, ...
            lsmap{k,3},'Color',lsmap{k,2},'LineWidth',2.5);
    end
end
xlabel('Age (years)','FontSize',11); ylabel('40Hz ITPC','FontSize',11);
title('ITPC vs Age (under 18)','FontSize',11,'FontWeight','bold');
legend({'TD','ASD w/ID','PMS'},'Location','northwest','FontSize',9,'Box','off');
set(gca,'Box','off','TickDir','out','FontSize',10);
ylim([0 0.65]); xlim([0 19]); hold off;

subplot(1,2,2); hold on;
data_m2   = [T2_TD.ITPC; T2_ASD.ITPC; T2_PMS.ITPC];
labels_m2 = [repmat({'TD'},height(T2_TD),1); repmat({'ASD'},height(T2_ASD),1); repmat({'PMS'},height(T2_PMS),1)];
bp5 = boxplot(data_m2,labels_m2,'Colors',[col_TD;col_ASD;col_PMS],'Symbol','o','Widths',0.5);
set(bp5,'LineWidth',1.5);
scatter(1+(rand(height(T2_TD),1)-0.5)*jitter,  T2_TD.ITPC,  25,col_TD,  'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(height(T2_ASD),1)-0.5)*jitter, T2_ASD.ITPC, 25,col_ASD, 'filled','MarkerFaceAlpha',0.5);
scatter(3+(rand(height(T2_PMS),1)-0.5)*jitter, T2_PMS.ITPC, 25,col_PMS, 'filled','MarkerFaceAlpha',0.5);
ylabel('40Hz ITPC','FontSize',11);
title(sprintf('ITPC by Group (under 18)\n(KW p=%.3f)',p_kw2),'FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',10); ylim([0 0.65]); hold off;

sgtitle('Model 2: TD + ASD w/ID + PMS (Under 18)','FontSize',13,'FontWeight','bold');


%% ============================================================
%  MODEL 3S: IQ sensitivity (ASD + PMS, IQ available)
%% ============================================================

fprintf('\n=== MODEL 3S: IQ sensitivity (ASD + PMS) ===\n');

idx3 = (is_ASD | (is_PMS & is_under18)) & ~isnan(T.IQ);
T3   = T(idx3, :);
T3.DX_ap = categorical(strcmp(T3.DX,'PMS'), [0 1], {'ASD','PMS'});

mdl3 = fitlm(T3, 'ITPC ~ logAge * DX_ap + Sex + IQ');
disp(mdl3);
fprintf('R2=%.3f, N=%d (ASD=%d, PMS=%d)\n', mdl3.Rsquared.Ordinary, height(T3), ...
    sum(strcmp(T3.DX,'ASD')), sum(strcmp(T3.DX,'PMS')));

for grpname = {'ASD','PMS'}
    sub = T3(strcmp(T3.DX, grpname{1}), :);
    m = fitlm(sub, 'ITPC ~ logAge + Sex + IQ');
    fprintf('Within-group %s: beta=%.3f, SE=%.3f, p=%.4f\n', ...
        grpname{1}, m.Coefficients.Estimate(2), m.Coefficients.SE(2), m.Coefficients.pValue(2));
end

[p_kw3, tbl_kw3] = kruskalwallis(T3.ITPC, T3.DX_ap, 'off');
fprintf('KW: H=%.3f, p=%.4f\n', tbl_kw3{2,5}, p_kw3);


%% ---- MODEL 3S PLOTS ----

figure('Name','Model 3S','Color','w','Position',[50 50 1100 480]);

subplot(1,2,1); hold on;
T3_ASD = T3(strcmp(T3.DX,'ASD'), :);
T3_PMS = T3(strcmp(T3.DX,'PMS'), :);
scatter(T3_ASD.Age, T3_ASD.ITPC, 40,col_ASD,'s','filled','MarkerFaceAlpha',0.5);
scatter(T3_PMS.Age, T3_PMS.ITPC, 40,col_PMS,'^','filled','MarkerFaceAlpha',0.5);

b3        = mdl3.Coefficients.Estimate;
mean_iq3  = mean(T3.IQ,'omitnan');
mean_sex3 = mean(T3.Sex);
y_ASD3 = b3(1) + b3(2)*log(age_range_fine) + b3(4)*mean_sex3 + b3(5)*mean_iq3;
y_PMS3 = b3(1) + (b3(2)+b3(6))*log(age_range_fine) + b3(3) + b3(4)*mean_sex3 + b3(5)*mean_iq3;
plot(age_range_fine, y_ASD3,'--','Color',col_ASD,'LineWidth',2.5);
plot(age_range_fine, y_PMS3,':','Color',col_PMS,'LineWidth',2.5);
xlabel('Age (years)','FontSize',11); ylabel('40Hz ITPC','FontSize',11);
title('ITPC vs Age (IQ covariate)','FontSize',11,'FontWeight','bold');
legend({'ASD','PMS'},'Location','northwest','FontSize',9,'Box','off');
set(gca,'Box','off','TickDir','out','FontSize',10);
ylim([0 0.65]); xlim([1 40]); hold off;

subplot(1,2,2); hold on;
data_m3   = [T3_ASD.ITPC; T3_PMS.ITPC];
labels_m3 = [repmat({'ASD'},height(T3_ASD),1); repmat({'PMS'},height(T3_PMS),1)];
bp6 = boxplot(data_m3,labels_m3,'Colors',[col_ASD;col_PMS],'Symbol','o','Widths',0.5);
set(bp6,'LineWidth',1.5);
scatter(1+(rand(height(T3_ASD),1)-0.5)*jitter, T3_ASD.ITPC, 25,col_ASD,'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(height(T3_PMS),1)-0.5)*jitter, T3_PMS.ITPC, 25,col_PMS,'filled','MarkerFaceAlpha',0.5);
ylabel('40Hz ITPC','FontSize',11);
title(sprintf('ASD vs PMS (IQ controlled)\n(KW p=%.3f)',p_kw3),'FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',10); ylim([0 0.65]); hold off;

sgtitle('Model 3S: IQ Sensitivity','FontSize',13,'FontWeight','bold');


%% ============================================================
%  SEX EFFECTS POST-HOC (N=131, all participants)
%% ============================================================

fprintf('\n=== SEX EFFECTS POST-HOC ===\n');

% Reload full data for sex model (includes PMS and ASD w/ID adults)
T_sex = readtable('DATA_SB.xlsx');
T_sex.IQ     = str2double(string(T_sex.IQ));
T_sex.logAge = log(T_sex.Age);
T_sex.Sex    = double(strcmp(T_sex.Gender,'F'));
T_sex = T_sex(~isnan(T_sex.Age), :);   % N=131

T_sex.Sex_cat = categorical(T_sex.Sex, [0 1], {'M','F'});
mdl_sex = fitlm(T_sex, 'ITPC ~ logAge * Sex_cat');
disp(mdl_sex);
fprintf('N=%d (M=%d, F=%d)\n', height(T_sex), sum(T_sex.Sex==0), sum(T_sex.Sex==1));

is_TD_sex  = strcmp(T_sex.DX,'TD');
is_ASD_sex = strcmp(T_sex.DX,'ASD');
is_PMS_sex = strcmp(T_sex.DX,'PMS');

TD_M = T_sex(is_TD_sex  & T_sex.Sex==0, :);
TD_F = T_sex(is_TD_sex  & T_sex.Sex==1, :);
ASD_M= T_sex(is_ASD_sex & T_sex.Sex==0, :);
ASD_F= T_sex(is_ASD_sex & T_sex.Sex==1, :);
PMS_M= T_sex(is_PMS_sex & T_sex.Sex==0, :);
PMS_F= T_sex(is_PMS_sex & T_sex.Sex==1, :);

[p_sex_td,  tbl_sex_td]  = kruskalwallis([TD_M.ITPC;  TD_F.ITPC],  [zeros(height(TD_M),1);  ones(height(TD_F),1)],  'off');
[p_sex_asd, tbl_sex_asd] = kruskalwallis([ASD_M.ITPC; ASD_F.ITPC], [zeros(height(ASD_M),1); ones(height(ASD_F),1)], 'off');
[p_sex_pms, tbl_sex_pms] = kruskalwallis([PMS_M.ITPC; PMS_F.ITPC], [zeros(height(PMS_M),1); ones(height(PMS_F),1)], 'off');

fprintf('KW TD  M vs F: H=%.3f, p=%.4f (M n=%d, F n=%d)\n', tbl_sex_td{2,5},  p_sex_td,  height(TD_M),  height(TD_F));
fprintf('KW ASD M vs F: H=%.3f, p=%.4f\n', tbl_sex_asd{2,5}, p_sex_asd);
fprintf('KW PMS M vs F: H=%.3f, p=%.4f\n', tbl_sex_pms{2,5}, p_sex_pms);


%% ---- SEX PLOTS ----

figure('Name','Sex Effects','Color','w','Position',[50 50 1300 480]);

subplot(1,3,1); hold on;
T_M = T_sex(T_sex.Sex==0, :);
T_F = T_sex(T_sex.Sex==1, :);
scatter(T_M.Age, T_M.ITPC, 35,col_M,'o','filled','MarkerFaceAlpha',0.35);
scatter(T_F.Age, T_F.ITPC, 35,col_F,'s','filled','MarkerFaceAlpha',0.35);

b_sex = mdl_sex.Coefficients.Estimate;
y_M = b_sex(1) + b_sex(2)*log(age_range_fine);
y_F = b_sex(1) + b_sex(2)*log(age_range_fine) + b_sex(3) + b_sex(4)*log(age_range_fine);
plot(age_range_fine, y_M, '-',  'Color',col_M, 'LineWidth',2.5);
plot(age_range_fine, y_F, '--', 'Color',col_F, 'LineWidth',2.5);
xlabel('Age (years)','FontSize',11); ylabel('40Hz ITPC','FontSize',11);
title('ITPC vs Age by Sex','FontSize',11,'FontWeight','bold');
legend({'Male','Female'},'Location','northwest','FontSize',9,'Box','off');
set(gca,'Box','off','TickDir','out','FontSize',10);
ylim([0 0.65]); xlim([1 40]); hold off;

subplot(1,3,2); hold on;
data_sex_td   = [TD_M.ITPC; TD_F.ITPC];
labels_sex_td = [repmat({'M'},height(TD_M),1); repmat({'F'},height(TD_F),1)];
bp7 = boxplot(data_sex_td,labels_sex_td,'Colors',[col_M;col_F],'Symbol','o','Widths',0.5);
set(bp7,'LineWidth',1.5);
scatter(1+(rand(height(TD_M),1)-0.5)*jitter, TD_M.ITPC, 25,col_M,'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(height(TD_F),1)-0.5)*jitter, TD_F.ITPC, 25,col_F,'filled','MarkerFaceAlpha',0.5);
if p_sex_td < 0.05
    y_sig = max(data_sex_td)+0.03;
    plot([1 2],[y_sig y_sig],'k-','LineWidth',1.2);
    text(1.5,y_sig+0.01,'*','HorizontalAlignment','center','FontSize',14);
end
ylabel('40Hz ITPC','FontSize',11);
title(sprintf('TD: Male vs Female\n(KW p=%.3f)',p_sex_td),'FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',10); ylim([0 0.7]); hold off;

subplot(1,3,3); hold on;
data_sex_all   = [ASD_M.ITPC; ASD_F.ITPC; PMS_M.ITPC; PMS_F.ITPC];
labels_sex_all = [repmat({'ASD M'},height(ASD_M),1); repmat({'ASD F'},height(ASD_F),1); ...
                  repmat({'PMS M'},height(PMS_M),1); repmat({'PMS F'},height(PMS_F),1)];
bp8 = boxplot(data_sex_all,labels_sex_all,'Symbol','o','Widths',0.5);
set(bp8,'LineWidth',1.5);
scatter(1+(rand(height(ASD_M),1)-0.5)*jitter, ASD_M.ITPC, 20,col_ASD*0.8,'filled','MarkerFaceAlpha',0.5);
scatter(2+(rand(height(ASD_F),1)-0.5)*jitter, ASD_F.ITPC, 20,col_ASD,    'filled','MarkerFaceAlpha',0.5);
scatter(3+(rand(height(PMS_M),1)-0.5)*jitter, PMS_M.ITPC, 20,col_PMS*0.8,'filled','MarkerFaceAlpha',0.5);
scatter(4+(rand(height(PMS_F),1)-0.5)*jitter, PMS_F.ITPC, 20,col_PMS,    'filled','MarkerFaceAlpha',0.5);
ylabel('40Hz ITPC','FontSize',11);
title(sprintf('ASD & PMS: M vs F\n(ASD p=%.3f | PMS p=%.3f)',p_sex_asd,p_sex_pms),'FontSize',11,'FontWeight','bold');
set(gca,'Box','off','TickDir','out','FontSize',10); ylim([0 0.65]); hold off;

sgtitle('Sex Effects Post-Hoc (N=131)','FontSize',13,'FontWeight','bold');


%% ============================================================
%  TRIAL COUNT ANALYSIS
%% ============================================================

fprintf('\n=== TRIAL COUNT ANALYSIS ===\n');

seg_grps = {
    'TD',                   T(is_TD, :);
    'PMS Under 18',         T(is_PMS & is_under18, :);
    'ASD w/ID Under 18',    T(is_ASD & is_ID & is_under18, :);
    'ASD w/o ID Under 18',  T(is_ASD & is_noID & is_under18, :);
    'ASD w/o ID Adults',    T(is_ASD & is_noID & is_adult, :);
};

seg_data  = {};
seg_label = {};
for g = 1:size(seg_grps,1)
    name = seg_grps{g,1};
    grp  = seg_grps{g,2};
    fprintf('  %s: n=%d, M=%.2f, SD=%.2f\n', name, height(grp), mean(grp.Segments), std(grp.Segments));
    seg_data{g}  = grp.Segments;
    seg_label{g} = repmat({name}, height(grp), 1);
end

all_segs   = cell2mat(seg_data');
all_labels = vertcat(seg_label{:});
[p_seg_kw, tbl_seg_kw] = kruskalwallis(all_segs, all_labels, 'off');
fprintf('\nKW trial count: H=%.3f, p=%.4f\n', tbl_seg_kw{2,5}, p_seg_kw);

[r_seg_itpc, p_seg_itpc] = corr(T.Segments, T.ITPC, 'Type','Pearson');
fprintf('Pearson Segments vs ITPC: r=%.3f, p=%.4f\n', r_seg_itpc, p_seg_itpc);

mdl_age_only = fitlm(T, 'ITPC ~ logAge');
itpc_resid   = mdl_age_only.Residuals.Raw;
[r_partial, p_partial] = corr(T.Segments, itpc_resid, 'Type','Pearson');
fprintf('Partial r (controlling log Age): r=%.3f, p=%.4f\n', r_partial, p_partial);


%% ============================================================
%  SUMMARY
%% ============================================================

fprintf('\n\n======================================\n');
fprintf('KEY RESULTS SUMMARY\n');
fprintf('======================================\n');
fprintf('Model 0A (KW all groups collapsed): H=%.3f, p=%.4f (ns)\n', tbl_kw0a{2,5}, p_kw0a);
fprintf('Model 0B (OLS all groups): F(6,121)=5.164, R2=0.204, p=0.0001\n');
fprintf('  logAge (TD slope): beta=0.083, p<0.001\n');
fprintf('  logAge*ASD interaction: beta=-0.067, p=0.013\n');
fprintf('  logAge*PMS interaction: beta=-0.053, p=0.121 (ns)\n');
fprintf('Model 1 (TD vs ASD w/o ID): F(4,75)=6.592, R2=0.260, p=0.0001\n');
fprintf('  TD within: beta=0.080, p=0.008\n');
fprintf('  ASD within: beta=0.024, p=0.241\n');
fprintf('  KW adults: p=0.045 | KW children: p=0.451\n');
fprintf('Model 2 (TD+ASD w/ID+PMS <18): F(6,61)=1.343, R2=0.117, p=0.252\n');
fprintf('Model 3S (IQ sensitivity): F(5,69)=1.498, R2=0.098, p=0.202\n');
fprintf('  IQ: p=0.927\n');
fprintf('Sex (N=131): F(3,127)=6.617, R2=0.135, p=0.0003\n');
fprintf('  TD M>F: KW H=7.09, p=0.008\n');
fprintf('Trial count vs ITPC: r=%.3f, p=%.4f\n', r_seg_itpc, p_seg_itpc);
fprintf('Partial r (age-controlled): r=%.3f, p=%.4f\n', r_partial, p_partial);
