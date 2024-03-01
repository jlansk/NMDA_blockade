function E = environment_blk()

% This function sets the environment for the cmm_NMDA_blk analyses
scr =     '/imaging/rowe/users/jl01/meg/dcm_cmm/reps0-5/scripts/Mg/to_publish/';
addpath(genpath(scr));
addpath(genpath('/imaging/rowe/users/jl01/meg/dcm_cmm/reps0-5/'), '-begin') %MAKE SURE THIS IS TOP OF PATH so correct cmm_NMDA.m functions are used

ana='/imaging/rowe/users/jl01/meg/dcm_cmm/reps0-5/';

% part 1 memcons:
anaD = [ana 'mem_cons/publish'];
rawD='/imaging/rowe/users/jl01/meg/source_mmn/data/mem_cons';

% part 2 ntad:
raw='/imaging/rowe/users/jl01/meg/dcm_cmc/meg_data/';
anaB= [ana '/BL_Feb23/publish/'];
anaL = [ana '/AF_Feb23/publish/'];

% fill in fields of environment structure
%----------------------------------------------------------------------
E.scr = scr;
E.ana = ana;
E.anaD = anaD;
E.rawD = rawD;
E.anaB = anaB;
E.anaL = anaL;
E.raw = raw;

end
