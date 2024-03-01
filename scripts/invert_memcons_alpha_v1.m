%% DCM using the cmc model to test difference in AD vs controls
% 'ntad_mmn_preprocessing_v2_jlrep10.m' these files also had source recon
% done ('mmn_source_recon_jl.m') so MRI is already linked to the file (same
% as cmc paper preprocessing)

%% Set up environment
clearvars
E = environment_blk;

%% Checks - make sure you are using the edited versions
open spm_fx_cmm_NMDA.m
open spm_cmm_NMDA_priors.m
open mg_switch.m

%% Set up variables
scr = E.scr;
anaD = E.anaD;
rawD = E.rawD;
load([scr filesep 'mem_cons']);
subjects=mem_cons;

session = 'mem'; %'pla';
anaDs= [anaD session];
rawlast=['/' session '/mmn/ffmraeMaffffdtsss.mat']; % file that has had all NTAD preprocessing steps applied) and link to mri


%% Set up  model space
%create models from my model space with different A and C matrices
[DCMa] = dcm_cmm_ipc_gen_model_space();

% specify rest of parameters for the full DCM model
for ss=1:length(subjects)
    dcms{ss} = dcm_cmm_ipc_gen(DCMa{end}, rawD, anaDs, subjects{ss}, rawlast);
end

dcms=dcms';
GCM=spm_dcm_load(dcms);

%% Invert first-level DCMs
%  subject that do not converge go into the 'failed' variable
   
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
    save([anaDs '/failedsubs.mat'], 'failed_subs')
    

%% Remove failed or bad fits
% From the variable 'failed_subs' created during the inversion and from
% seeing which have NANs in the H matrix (as shown through above plots) we
% remove those failed subjects

anaP = [anaD filesep 'pla/'];
anaM = [anaD filesep 'mem/'];

nan_subs_P={}; i=1;
for ss=1:length(subjects)
    try
        load([anaP '/DCM_' subjects{ss} '_full.mat']);
        if any(isnan(DCM.R{1}))
            nan_subs_P(i)=subjects(ss);
            i=i+1;
        end
    catch;
    end
end

nan_subs_M={}; i=1;
for ss=1:length(subjects)
    try
        load([anaM '/DCM_' subjects{ss} '_full.mat']);
        if any(isnan(DCM.R{1}))
            nan_subs_M(i)=subjects(ss);
            i=i+1;
        end
    catch
    end
end

failed_subs_M=load([anaM 'failedsubs.mat']); failed_subs_M=failed_subs_M.failed_subs;
failed_subs_P=load([anaP 'failedsubs.mat']); failed_subs_P=failed_subs_P.failed_subs;

failedsubs=[failed_subs_M; failed_subs_P]';

i=0;
for ss=1:length(subjects)
    if ~any(find(strcmp(subjects{ss},[nan_subs_P, nan_subs_M, failedsubs])))
        i=i+1;
        
        DCMsubs{i}=subjects{ss};
    end
end

if ~exist([scr '/DCMsubs.mat'], 'file')
    save([scr '/DCMsubs'], 'DCMsubs');
end


%% Plot model fits (line graphs)
    
    % placebo model fits
    figure(1)
        
    for ss = 1:length(DCMsubs)
        try
            DCM =spm_dcm_load([anaP filesep 'DCM_' DCMsubs{ss} '_full.mat']);
            DCM = DCM{1};
            reps=1;
            
            subplot(2,9,ss)
            
            for c = 1 %:5:size((DCM.H),2)
                plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
                plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
                
                obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
                prd1(ss,:)=DCM.H{c}(:,1);
                
                cortemp=corrcoef(DCM.H{c}(:,1) + DCM.R{c}(:,1), DCM.H{c}(:,1));
                corL.s(ss,1)=DCMsubs(ss);
                corL.d(ss,1)=cortemp(2);
                
                ylim([-7 4]); yticks([-6 4]);
                xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
                xlabel(extractBetween(DCM.name, '/DCM_', '_'));
                set(gcf, 'color', 'w');
                set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
            end
            box off
            
        catch; continue
        end
    end
    legend({'DEV(pred)', 'DEV(obs)'}); legend('boxoff');
    
    if ~exist([scr '/figures/placebo_model_fits.png'], 'file')
        exportgraphics(gcf, [scr '/figures/placebo_model_fits.png'], 'Resolution', '720');
    end
    
    %Memantine model fits
    figure(2)
        
    for ss = 1:length(DCMsubs)
        try
            DCM =spm_dcm_load([anaM filesep 'DCM_' DCMsubs{ss} '_full.mat']);%'_full.mat']);%(dcms{ss});
            DCM = DCM{1};
            reps=1;
            
            subplot(2,9,ss)
            
            for c = 1 %:5:size((DCM.H),2)
                plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', [0.4 0.4 0.4], 'Linewidth', 1.5); hold on
                plot(DCM.H{c}(:,1), 'color', [0.3    0.8    0.6510], 'Linewidth', 1.5);
                
                obs1(ss,:)=DCM.H{c}(:,1) + DCM.R{c}(:,1);
                prd1(ss,:)=DCM.H{c}(:,1);
                
                cortemp=corrcoef(DCM.H{c}(:,1) + DCM.R{c}(:,1), DCM.H{c}(:,1));
                corL.s(ss,1)=DCMsubs(ss);
                corL.d(ss,1)=cortemp(2);
                
                ylim([-7 4]); yticks([-6 4]);
                xlim([0 125]); xticks([25 75]); %xticklabels({'0', '200'});
                xlabel(extractBetween(DCM.name, '/DCM_', '_'));
                set(gcf, 'color', 'w');
                set(gcf, 'Position', [100 + 400*(reps-1) 100 400 800]);
            end
            box off            
        catch; continue
        end
    end
    legend({'DEV(pred)', 'DEV(obs)'}); legend('boxoff');
    
    
    if ~exist([scr '/figures/placebo_model_fits.png'], 'file')
        exportgraphics(gcf, [scr '/figures/memantine_model_fits.png'], 'Resolution', '720');
    end
    
