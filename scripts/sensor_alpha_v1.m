%% Sensor level analysis for cmm - Juliette Lanskey 2023

% set up variables
clearvars
E = environment_blk;

addpath(genpath('/imaging/rowe/users/jl01/meg/dcm_cmc/scripts/to_publish/'))%cmc paper has some functions (gitlab/jlansk/dcm_cmc_ntad

raw = E.raw; %directory for MEG files
scr = E.scr; % directory of scripts
ana_dir = raw; % where the pre-processed MEG files are

load([scr filesep 'BLsubs.mat']); % 1x56 cell with baseline IDs
load([scr filesep 'AFsubs.mat']); % 1x33 cell with follow-up IDs

pat = find(contains(BLsubs, 'P'));
BLsubs = BLsubs(pat);

task='mmn';
session1 = 'BL';
session2 = 'AF';
lfile= 'bCPffmraeMaffffdtsss.mat'; % file with preprocessing and combined grads

%% load baseline MMN amplitude
[mmn_pa_BL, mmnBL] = mmn_amp(BLsubs, ana_dir, 'BL', task, lfile, 0);

%% load follow-up MMN amplitude
[mmn_pa_AF, mmnAF] =mmn_amp(AFsubs, ana_dir, 'AF', task, lfile, 0); %12 is dev-stn, not needed

%% statistics
load([scr filesep 'BL_mmse_acer']);
BL_acer=BL_mmse_acer;

for ss=1:length(BLsubs)
    idx=find(contains(BL_acer(:,1),BLsubs(ss)));
    Xacer(ss,1) = BL_acer(idx,1);
    Xacer(ss,2) = {BL_acer{idx,2}/100}; % should be the value corresponding to the peb_subs
end
%% run lm model

% Plot without outliers
idx=find(isoutlier(mean(mmn_pa_BL(:,1:10),2))==0); % remove 2 patients with outlying MMSE (
mdl = fitlm([Xacer{idx,2}]*100, (mean(mmn_pa_BL(idx,1:10),2)))

plot(mdl); hold on
scatter([Xacer{idx,2}]*100, (mean(mmn_pa_BL(idx,1:10),2)), "filled", "black"); box off; set(gcf, 'color', 'white'); legend off; %make it prettier!

title(['slope=' num2str(mdl.Coefficients{2,1}) ' p=' num2str(mdl.Coefficients{2,4})]);
ylabel('Mismatch amplitude (fT/m)'); xlabel('MMSE');

if ~exist([scr '/figures/MMSE_blk_regression.png'], 'file')
    exportgraphics(gcf, [scr '/figures/MMSE_blk_regression.png'], 'Resolution', '720');
end

%% plot with outliers to show their info
%one person with low MMSE=18 and one with high MMSE=27 both have low MMN
gscatter([Xacer{:,2}]*100, mean(mmn_pa_BL(:,1:10),2), isoutlier(mean(mmn_pa_BL(:,1:10),2)), [0 0 0; 1 0 0]); box off; legend off; xlabel('MMSE'); ylabel('MMN amplitude'); title('Outliers in red');
%exportgraphics(gcf, [scr '/figures/MMSE_blk_regression_overlayedoutliers.png'], 'Resolution', '720');

hold off; gscatter([Xacer{:,2}]*100, mean(mmn_pa_BL(:,1:10),2), isoutlier(mean(mmn_pa_BL(:,1:10),2)), [0 0 0; 1 0 0]); box off; legend off; xlabel('MMSE'); ylabel('MMN amplitude'); title('Outliers in red');
%exportgraphics(gcf, [scr '/figures/MMSE_blk_scatter_outliers.png'], 'Resolution', '720');

% % The same 2 patients are outliers for MMN amplitude at baseline and follow up
idx_outlier_BL=isoutlier(mean(mmn_pa_BL(:,1:10),2));
idx_outlier_AF=isoutlier(mean(mmn_pa_AF(:,1:10),2));

meanpa(:,1)=mean(mmn_pa_BL(:,1:10),2);
meanpa(size(mmn_pa_BL,1)+1:size(mmn_pa_BL,1)+size(mmn_pa_AF,1),1)=mean(mmn_pa_AF(:,1:10),2);
meanpa(size(mmn_pa_BL,1)+1:size(mmn_pa_BL,1)+size(mmn_pa_AF,1),2)=1;
meanpa(:,3)=[idx_outlier_BL; idx_outlier_AF]
gscatter(meanpa(:,2), meanpa(:,1), meanpa(:,3), [0 0 0; 1 0 0]); % or copy mean_pa_ids into an excel sheet, save and use in R for lines between dots

% For lines between the above scatter, copy mean_pa_ids into excel, save, and run
% R script (mmn_amplitude_bardotline
mmn_amplitude_dotplot(:,1)=BLsubs; mmn_amplitude_dotplot(length(BLsubs)+1:length(BLsubs)+length(AFsubs))=AFsubs;
mmn_amplitude_dotplot(:,2:4) = num2cell(meanpa(:,1:3));


%% exclude the 2 outliers
BLsubs=BLsubs(idx);

%% Longitudinal comparison
Lsubs = AFsubs(ismember(AFsubs, BLsubs));
LBsubs = find(contains(BLsubs, Lsubs));
LAsubs = find(contains(AFsubs, Lsubs));

[h,p,ci,stats]=ttest(mean(mmn_pa_BL(LBsubs,1:10),2), mean(mmn_pa_AF(LAsubs,1:10),2), 'tail','left');
bnanmean=nanmean(mmn_pa_BL(LBsubs,5)); fnanmean=nanmean(mmn_pa_AF(LAsubs,5));

AFstatistics(1,1)=stats.tstat; AFstatistics(2,1)=p; AFstatistics(3,1)=bnanmean; AFstatistics(4,1)=fnanmean;

%% Longitudinal plot
BL_av = squeeze(nanmean(mmnBL(LBsubs,:,:),1))
BL_se= squeeze(nanstd(mmnBL(LBsubs,:,:),1))./sqrt(size(LBsubs,2));

AF_av = squeeze(nanmean(mmnAF(LAsubs,:,:),1))
AF_se = squeeze(nanstd(mmnAF(LAsubs,:,:),1))./sqrt(size(LAsubs,2))

%% longitudinal plot
close all
lh(1)=plot(BL_av(:, 5), 'Color', [0 0 0.6]); hold on; lh(2)=plot(AF_av(:, 5), 'Color', [0.05 0.6 0.65]); hold on;
set(lh(1), 'linewidth', 1); set(lh(2), 'linewidth', 1);
boundedline([1:size(BL_av)],BL_av(:,5),BL_se(:,5), 'transparency', 0.333, 'cmap', [0 0.1 0.7], 'alpha'); boundedline([1:size(AF_av)],AF_av(:,5),AF_se(:,5), 'transparency', 0.333, 'cmap',  [0.05 0.6 0.75], 'alpha'); %[0.05 0.8 0.85]see uisetcolor for colour options
xlim([0 250]); xticks([0 50 100 150 200 250]); xticklabels({'-100','0', '100', '200', '300', '400'}); xlabel('Time (ms)');
xlabel('Time (ms)'); ylabel('Mismatch response (fT/m)');
%  patch([beg2,beg2, fin2, fin2, beg2], [0.1, -1.0, -1.0, 0.1, 0.1], 'black', 'EdgeColor', 'none', 'FaceColor', 'black', 'FaceAlpha', 0.1);
begin2=140; beg2=(begin2+100)/2;
finish2=160 ; fin2=(finish2+100)/2;
ylim([-0.7 0.1]); box off; set(gcf, 'color', 'w');
legend('baseline', 'follow-up'); legend('Location', 'best'); legend('boxoff')
titlestr={'MMN_r1r5'};

if ~exist([scr filesep '/figures/BF_MMN_05.png'], 'file')
    exportgraphics(gcf, [scr filesep '/figures/BF_MMN_05.png'], 'Resolution', 720) %only works in 2020
end

%% all together
BL_av = squeeze(nanmean(nanmean(mmnBL(LBsubs,:,1:10),3),1));
BL_se= squeeze(nanstd(mean(mmnBL(LBsubs,:,1:10),3),1))./sqrt(size(LBsubs,2));

AF_av = squeeze(nanmean(nanmean(mmnAF(LAsubs,:,1:10),3),1))
AF_se = squeeze(nanstd(mean(mmnAF(LAsubs,:,1:10),3),1))./sqrt(size(LAsubs,2))

% need to check this as think i need to average across reps still...

close all
lh(1)=plot(BL_av(:), 'Color', [0 0 0.6]); hold on; lh(2)=plot(AF_av(:), 'Color', [0.05 0.6 0.65]); hold on;
set(lh(1), 'linewidth', 1); set(lh(2), 'linewidth', 1);
boundedline([1:length(BL_av)],BL_av(:),BL_se(:), 'transparency', 0.333, 'cmap', [0 0.1 0.7], 'alpha'); boundedline([1:length(AF_av)],AF_av(:),AF_se(:), 'transparency', 0.333, 'cmap',  [0.05 0.6 0.75], 'alpha'); %[0.05 0.8 0.85]see uisetcolor for colour options

xlim([0 250]); xticks([0 50 100 150 200 250]); xticklabels({'-100','0', '100', '200', '300', '400'}); xlabel('Time (ms)');
xlabel('Time (ms)'); ylabel('Mismatch response (fT/m)');
%  patch([beg2,beg2, fin2, fin2, beg2], [0.1, -1.0, -1.0, 0.1, 0.1], 'black', 'EdgeColor', 'none', 'FaceColor', 'black', 'FaceAlpha', 0.1);
begin2=140; beg2=(begin2+100)/2;
finish2=160 ; fin2=(finish2+100)/2;
%ylim([-0.7 0.1]); 
box off; set(gcf, 'color', 'w');
legend('baseline', 'follow-up'); legend('Location', 'best'); legend('boxoff')
titlestr={'MMN_rs'};

if ~exist([scr filesep '/figures/BF_MMN_rs.png'], 'file')
    exportgraphics(gcf, [scr filesep '/figures/BF_MMN_rs.png'], 'Resolution', 720) %only works in 2020
end