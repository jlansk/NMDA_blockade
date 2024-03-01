%% Assess the memantine effect 
% Runs a series of PEBs to identify where memantine acts and direction
% of effect

%% Set up environment
clear all
E = environment_blk;

%% Set up variables
scr = E.scr;
anaD = E.anaD;
anaP = [anaD filesep 'pla/'];
anaM = [anaD filesep 'mem/'];

load([anaD filesep 'DCMsubs.mat']);
memcons=DCMsubs;

%% Load the DCMs into a GCM file - placebo first (anaP) then memantine (anaM)
for ss=1:length(memcons)
    GCM{ss,1}=[anaP filesep 'DCM_' memcons{ss} '_full.mat']; 
    GCM{length(memcons)+ss,1}=[anaM filesep 'DCM_' memcons{ss} '_full.mat'];
end

GCM=spm_dcm_load(GCM);


%% Baseline Group comparison (PEB analysis)
% specify PEB model
X=zeros(2*length(memcons),2);
X(:,1)=ones; %first collumn to ones to model overall mean

X(1:length(memcons),2)= 0;
X(length(memcons)+1:end,2) = 1;

X_labels={'mean','pla->mem'};

M=struct(); %a structure to specify the PEB settings
M.Q = 'all'; % between-subject variability will be estimated for each connection
M.X=X; %design matrix
M.Xnames = X_labels;

field{1} = {'T(:,1)'};
field{2} = {'T(:,2)'}; 
field{3} = {'T(:,3)'};
field{4} = {'Mg'};
field{5} = {'H'};
field{6}=[field{1} field{2} field{3} field{4} field{5}]

labels = {'AMPA-T',  'GABA-T', 'NMDA-T', 'NMDA-Blk', 'H', 'All (T H Blk)'};

%--------------------------------------------------------------------------
for c = 1:length(field)
    [PEB{c}, Pmem{c}] = spm_dcm_peb(GCM, M, field{c});
    F(c)   = PEB{c}.F;
end

%% Plot Bayesian model comparison of all PEB models
%--------------------------------------------------------------------------
plot_bmc_F_jl(F, labels);

if ~exist([scr '/figures/memantine_comparison.png'], 'file')
    exportgraphics(gcf, [scr '/figures/memantine_comparison_all.png'], 'Resolution', '720');
end

win=find(F==max(F));
spm_dcm_peb_review_fig_jl(PEB{win}, 2, 0.95, 1);

if ~exist([scr '/figures/memantine_PEB_95.png'], 'file')
    exportgraphics(gcf, [scr '/figures/memantine_PEB_95.png'], 'Resolution', '720');
end

for ss = 1:length(Pmem{win})
    mem_Ep(ss, 1) = Pmem{win}{ss}.Ep.Mg(1);
    mem_Ep(ss,3) = Pmem{win}{ss}.Cp(1,1); %precision
end
mem_Ep(:, 2) = M.X(:,2); %group

scatter(mem_Ep(:,2), mem_Ep(:,1));

mem_Ep_exp = [num2cell(mem_Ep)]
mem_Ep_exp(:, 3) = [DCMsubs'; DCMsubs'];


%% Wilcoxin signed rank test
[p, h, stats] = signrank(mem_Ep(18:34,1), mem_Ep(1:17,1), 'tail', 'right');

