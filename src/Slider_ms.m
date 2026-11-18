% Earthquake Cycle THMC 1.1.0: current THMC* manuscript configuration.
% Yuhan Wang, ORCID 0000-0003-0490-5977. See CITATION.cff.
release_root=fileparts(fileparts(mfilename('fullpath')));
addpath(release_root);
result=run_thmc('THMC*');
