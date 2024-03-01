%% run between group PEB of winning model
%% Set up variables
clearvars
addpath(genpath('/imaging/rowe/users/jl01/meg/dcm_cmc/scripts/to_publish'))
E = environment_blk();

scr = E.scr;
anaB = E.anaB;
anaL = E.anaL;

load([anaB filesep 'Bsubs.mat']); %or BMC_BMA_2 and rcm.mat
age_cov=1;
days_cov=1;
load([anaL 'Lsubs.mat']);

field = {'Mg'};
labels = { 'Mg'};

%% Load DCMs
for ss=1:length(Bsubs)
    GCM{ss,1}=[anaB filesep 'DCM_' Bsubs{ss} '_full.mat'];
end
BGCM=spm_dcm_load(GCM);

%% ACE-R PEB
load([scr 'BL_mmse_acer.mat']);
BL_acer=BL_mmse_acer;

clearvars Xacer
for ss=1:length(Bsubs)
    try
        idx=find(contains(BL_acer(:,1),Bsubs(ss)));
        Xacer(ss,1) = BL_acer(idx,1);
        Xacer(ss,2) = {BL_acer{idx,2}/100};
    catch; continue
    end
end

clearvars peb_acer_subs
count=1;
for ss=1:length(Xacer)
    if size(Xacer{ss,2},2)==1;
        peb_acer_subs{count,1}=Xacer{ss,1};
        Xac(count,2)=Xacer{ss,2};
        
        if strcmp(extractBetween(BGCM{ss}.name, 'DCM_','_'), Xacer{ss,1})
            DCMwinAcer{count}=BGCM{ss};
        else
            disp('subject order for DCMwin and Xptau do not match')
        end
        
        count=count+1;
    else
        excluded{ss}=BGCM{ss}.name
    end
end
peb_acer_subs = peb_acer_subs';

%% specify PEB model
pat = find(contains(peb_acer_subs, 'P'));    
pDCMwinAcer = DCMwinAcer(pat);

pXac = Xac(pat,:); % patients only
pXac(:,1)=ones; %first collumn to ones to model overall mean
pXac(:,2:end)=zscore(pXac(:,2:end)); % Do not invert ace-r score, as may confuse people familiar with acer
Xac_labels={'mean','BL_acer'}

clearvars Mac PEBaceG RCMaceG RMAaceG
Mac=struct(); %a structure to specify the PEB settings
Mac.Q = 'all'; % between-subject variability will be estimated for each connection
Mac.X=pXac; %design matrix
Mac.Xnames = Xac_labels;

[PEBace, Pace] = spm_dcm_peb(pDCMwinAcer', Mac, field);

spm_dcm_peb_review_fig_jl(PEBace, 2, 0.95, 1);
xticklabels({'left IPC', 'right IPC'});
title('BL patient MMSE, threshold = 0.95');

if ~exist([scr '/figures/PEB_mmse.png'], 'file')
    exportgraphics(gcf, [scr '/figures/PEB_mmse.png'], 'Resolution', '720');
end

%% extract PEB-informed dcm values for plots of Mg block vs MMSE
for ss = 1:length(Pace)
    ace_Ep(ss, 1) = Pace{ss}.Ep.Mg(1)
    ace_Ep(ss, 2) = Pace{ss}.Ep.Mg(2)
end
ace_Ep(:, 3) = Mac.X(:,2); %meancorrected mmse
ace_Ep(:, 4) = Xac(pat,2)*100; %mmse
scatter(ace_Ep(:,4), ace_Ep(:,2));

%% Longitudinal PEB analysis
count=1; clearvars LBsubs
for ss=1:length(Lsubs)
    if any(find(contains(Lsubs{ss}, Bsubs)))
        LBsubs(count)=Lsubs(ss);
        count=count+1;
    end
end

clearvars LGCM
for ss=1:length(LBsubs)
    LGCM{ss,1}=[anaB filesep 'DCM_' LBsubs{ss} '_full.mat'];
    LGCM{length(LBsubs)+ss,1}=[anaL filesep 'DCM_' LBsubs{ss} '_full.mat'];
end

LGCM=spm_dcm_load(LGCM);

LX=zeros(2*length(LBsubs),2);
LX(:,1)=ones; %first collumn to ones to model overall mean

if days_cov == 1
    load([scr '/days_af-bl.mat']) %loads into a variable called 'days'
    for ss=1:length(LBsubs)
        idx=find(contains(days(:,1),LBsubs(ss)));
        LX(length(LBsubs)+ss,3) = days{idx,2}/365.2425;
    end
    
    LX(length(LBsubs)+1:end,3) = zscore(LX(length(LBsubs)+1:end,3));
    LX_labels={'mean','BL->AF','days'};
    LX(length(LBsubs)+1:end,2) = 1;
else
    LX(length(LBsubs)+1:end,2) = 1;
    LX_labels={'mean','BL->AF'};
end

LM=struct(); %a structure to specify the PEB settings
LM.Q = 'all'; % between-subject variability will be estimated for each connection
LM.X=LX; %design matrix
LM.Xnames = LX_labels;

%% Run longitudinal PEB with block parameter (Mg)
[PEBL, PL] = spm_dcm_peb(LGCM, LM, field);
spm_dcm_peb_review_fig_jl(PEBL, 2, 0.95, 1)
xticklabels({'left IPC', 'right IPC'});

if ~exist([scr '/figures/PEB_L.png'], 'file')
    exportgraphics(gcf, [scr '/figures/PEB_L.png'], 'Resolution', '720');
end

% Add parameters into L_Ep which is used for bar plot (in R)
for ss = 1:length(PL)
    L_Ep(ss, 1) = PL{ss}.Ep.Mg(1)
    L_Ep(ss, 2) = PL{ss}.Ep.Mg(2)
    L_Ep(ss, 3) = LX(ss,2)
    
     CL_Ep(ss,1) = PL{ss}.Ep.Mg(2);
     CL_Ep(ss,2) = PL{ss}.Cp(2,2); 
     CL_Ep(ss, 3) = LX(ss,2)
end

%% Wilcoxin signed rank test
[p, h, stats] = signrank(L_Ep(1:28,1), L_Ep(29:56,1), 'tail', 'right');
[p, h, stats] = signrank(L_Ep(1:28,2), L_Ep(29:56,2), 'tail', 'right');
[p, h, stats] = signtest(L_Ep(1:28,2), L_Ep(29:56,2), 'tail', 'right');

%% Supplementary materials - Group PEB
X=zeros(length(Bsubs),2);
X(:,1)=ones; %first collumn to ones to model overall mean

con = find(contains(Bsubs, 'C'));
pat = find(contains(Bsubs, 'P'));

X(1:con(end),2)=0;%0; %set controls (index 1:14) to 0
X_labels={'mean','0cons -> 1pats','age'};
X(pat:end,2)=1; %set patients to 1

load('/imaging/rowe/users/jl01/meg/dcm_cmm/reps0-5/BL_Feb23/full/BLageMEG.mat')
for ss=1:length(Bsubs)
    idx=find(contains(BLage(:,1),Bsubs(ss)));
    Xage(ss,1) = BLage(idx,1);
    Xage(ss,2) = BLage(idx,2); % should be the value corresponding to the peb_subs
end

for ss=1:length(Bsubs)
    if strcmp(Xage{ss,1}, Bsubs{ss})
        X(ss,3)=Xage{ss,2}; %age if needed?
    else
        disp('subject orders do not match')
    end
end

X_labels={'mean','cons2pats','age'};
Xexp = X; %save X before age is zscored for plots
X(:,3:end)=zscore(X(:,3:end));

M=struct(); %a structure to specify the PEB settings
M.Q = 'all'; % between-subject variability will be estimated for each connection
M.X=X; %design matrix
M.Xnames = X_labels;

%% Run PEB
[PEB P] = spm_dcm_peb(BGCM, M, {'Mg'});

spm_dcm_peb_review_fig_jl(PEB, 2, 0.95, 1);
xticklabels({'left IPC', 'right IPC'});
exportgraphics(gcf, [scr '/figures/PEB_group_age1.png'], 'Resolution', '720');

spm_dcm_peb_review_fig_jl(PEB, 3, 0.95, 1);
xticklabels({'left IPC', 'right IPC'});
exportgraphics(gcf, [scr '/figures/PEB_group_age2.png'], 'Resolution', '720');

%% Extract PEB-informed DCM posterior values
for ss = 1:length(P)
    GEp(ss, 1) = P{ss}.Ep.Mg(1)
    GEp(ss, 2) = P{ss}.Ep.Mg(2)
end
GEp(:, 3) = M.X(:,2);
GEp(:, 4) = Xexp(:,3);

subplot(2,1,1); gscatter(GEp(:,4), GEp(:,1), GEp(:,3), [1 0.42 0.16; 0 0.32 0.7]); set(gcf, 'color', 'w'); box off; legend('controls', 'patients', 'Location', 'northeastoutside'); legend('boxoff'); title('Left Parietal NMDA blockade');
subplot(2,1,2); gscatter(GEp(:,4), GEp(:,2), GEp(:,3), [1 0.42 0.16; 0 0.32 0.7]); set(gcf, 'color', 'w'); box off; legend('controls', 'patients', 'Location', 'northeastoutside'); legend('boxoff'); title('Right Parietal NMDA blockade');

if ~exist([scr '/figures/PEB_NTAD_group_age.png'], 'file')
    exportgraphics(gcf, [scr '/figures/PEB_NTAD_group_age.png'], 'Resolution', '720');
end

%% suplementary materials - including age in the MMSE PEB
pXac(:,3)=[Xage{pat,2}];
pXac(:,2:end)=zscore(pXac(:,2:end)); 

clearvars Mac PEBaceG RCMaceG RMAaceG
Mac=struct(); %a structure to specify the PEB settings
Mac.Q = 'all'; % between-subject variability will be estimated for each connection
Mac.X=pXac; %design matrix
Mac.Xnames = [Xac_labels {'age'}];
[PEBace, Pace] = spm_dcm_peb(pDCMwinAcer', Mac, field);
spm_dcm_peb_review(PEBace)

for ss = 1:length(Pace)
    ace_Ep(ss, 1) = Pace{ss}.Ep.Mg(1)
    ace_Ep(ss, 2) = Pace{ss}.Ep.Mg(2)
end

ace_Ep(:, 3) = Mac.X(:,2); %meancorrected mmse
ace_Ep(:, 4) = Xac(pat,2)*100; %mmse
scatter(ace_Ep(:,4), ace_Ep(:,2));
box off; set(gcf, 'color', 'w');