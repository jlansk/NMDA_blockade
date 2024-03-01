function [DCMa] = dcm_cmm_2r_gen_model_space()
%creates A and C matrices from model space

LA1=1;
RA1=2;

Nareas=2;

da=1; %we only have one model

DCMa{da}.A{1} = zeros(Nareas, Nareas);   % forward connections

DCMa{da}.A{2} = zeros(Nareas, Nareas);    % backward connections

DCMa{da}.A{3} = zeros(Nareas,Nareas);    % reciprocal connections i.e. binary constraints on extrinsic connections
DCMa{da}.A{3}(LA1, LA1) = 1;
DCMa{da}.A{3}(RA1,RA1) = 1;

DCMa{da}.C = [1; 1];

%% Names

DCMa{da}.name = 'full'

end
