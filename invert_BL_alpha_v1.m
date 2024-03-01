%% DCM using the cmc model to test baseline differences (patients vs controls)
%  This version is an update of Lanskey cmc paper /jlansk/dcm_cmc_ntad
%  1. Inverts the full DCM for each individual and check fits

%% Set up environment
%DO NOT ADD SPM TO PATH - MAKE SURE IT ISN'T ADDED
clearvars
E = environment_blk();

%% Checks
%check magnesium block parameters are allowed to vary
open spm_cmm_NMDA_priors.m
open spm_fx_cmm_NMDA.m
open mg_switch.m

%% Set up variables
scr = E.scr;
ana = E.ana;
anaB = E.anaB;
raw = E.raw;
rawlast='/BL/mmn/ffmraeMaffffdtsss.mat'; % file that has had all NTAD preprocessing steps applied) and link to mri

load([scr filesep 'BLsubs']);
subjects=BLsubs; 

%% switches
invert = 1;
plot1 = 1;

%% Set up  model space

%create models from my model space with different A and C matrices
[DCMa] = dcm_cmm_ipc_gen_model_space();

% specify rest of parameters for the full DCM model
for ss=1:length(subjects)
    dcms{ss} = dcm_cmm_ipc_gen(DCMa{end}, raw, anaB, subjects{ss}, rawlast);    
end

dcms=dcms';
GCM=spm_dcm_load(dcms);

%% Invert first-level DCMs
%  subject that do not converge go into the 'failed' variable

if invert == 1
    
    failed={};
    parfor ss=1:length(subjects)
        try
            spm_dcm_erp(GCM{ss});
        catch
            failed{ss}={strcat('failed_', num2str(ss))};
            continue;
        end
    end
    
    fail = find(~cellfun(@isempty,failed));
    failed_subs = subjects(fail);
    save([anaB '/failedsubs.mat'],'failed_subs')
end

%% Plot model fits (line graphs)

if plot1 == 1
    
    try
        load([anaB '/Bsubs.mat'])
        subjects = Bsubs;
    catch; end
    
    dcms = cellstr(spm_select('FPList', [anaB filesep], '^*full.mat'));
    dcms = spm_dcm_load(dcms);
    figure
    title('NTAD patients model fit');
    
    count=0; nan_subs={};
    for ss = 15:length(subjects)
        try
            DCM =spm_dcm_load([anaB filesep 'DCM_' subjects{ss} '_full.mat']);
            DCM = DCM{1};
            if any(isnan([DCM.H{2}(:); DCM.H{1}(:)]))
                count=count+1
                nan_subs{count}=subjects{ss}
            end
            
            reps=1;
            
            subplot(6,14,ss)
            
            for c = 1 %:5:size((DCM.H),2)
                plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
                plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
                
                obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
                prd1(ss,:)=DCM.H{c}(:,1);
                
                ylim([-5 5]);
                
                xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
                xlabel(extractBetween(DCM.name, 'DCM_', '_'))
                set(gcf, 'color', 'w');
                set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
            end
            box off
        catch; continue
        end
    end
    legend({'DEV(pred)', 'DEV(obs)'}); legend('boxoff');
    
    if ~exist([scr '/figures/pats_BL_fits.png'], 'file')
        exportgraphics(gcf, [scr '/figures/pats_BL_fits.png'], 'Resolution', '720');
    end
        
    count=0; nan_subs={};
    for ss = 1:14%length(subjects)
        try
            DCM =spm_dcm_load([anaB filesep 'DCM_' subjects{ss} '_full.mat']);
            DCM = DCM{1};
            if any(isnan([DCM.H{2}(:); DCM.H{1}(:)]))
                count=count+1
                nan_subs{count}=subjects{ss}
            end
            
            reps=1;
            
            subplot(6,14,ss)
            
            for c = 1 %:5:size((DCM.H),2)
                plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
                plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
                
                obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
                prd1(ss,:)=DCM.H{c}(:,1);
                
                ylim([-5 5]);
                
                xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
                xlabel(extractBetween(DCM.name, 'DCM_', '_'))
                set(gcf, 'color', 'w');
                set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
            end
            box off
        catch; continue
        end
    end
    legend({'DEV(pred)', 'DEV(obs)'}); legend('boxoff');
    
    if ~exist([scr '/figures/cons_BL_fits.png'], 'file')
        exportgraphics(gcf, [scr '/figures/cons_BL_fits.png'], 'Resolution', '720');
    end
    
end
%% Remove failed or bad fits
% From the variable 'failed_subs' created during the inversion and from
% seeing which have NANs in the H matrix (as shown through above plots) we
% remove those failed subjects

try
        load([anaB '/Bsubs.mat'])
catch
    
    if ~exist('failed_subs', 'var')
        load([anaB '/failedsubs.mat'])
    end
    
    % run plots for nan_subs
    
    count=0;
    for ss=1:length(subjects)
        if ~any(find(contains(subjects{ss},[failed_subs, nan_subs])))
            count=count+1;
            
            Bsubs{count}=subjects{ss};
        end
    end
    
    save([anaB filesep 'Bsubs'], 'Bsubs');
    
end