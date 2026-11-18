function result=run_thmc(case_name,options)
%RUN_THMC Run the current manuscript code. Default: THMC* and Dirichlet.
% result=run_thmc('THMC*',struct('boundary_type','robin'));
% Yuhan Wang, ORCID 0000-0003-0490-5977. See CITATION.cff.
if nargin<1 || isempty(case_name),case_name='THMC*';end
if nargin<2,options=struct();end
validateattributes(options,{'struct'},{'scalar'});
root=fileparts(mfilename('fullpath'));
aliases={'M','cycles_0';'THM','cycles_001';'THMC','cycles_011'; ...
    'THM*','cycles_101';'THMC*','cycles_111';'THM_no_dilatancy','dense_000'};
case_name=char(case_name);
for i=1:size(aliases,1)
    if strcmpi(case_name,aliases{i,1}),case_name=aliases{i,2};break,end
end
allowed={'cycles_0','cycles_001','cycles_011','cycles_101','cycles_111', ...
    'dense_000','dense_001','dense_011'};
assert(ismember(case_name,allowed),'Unknown case: %s',case_name);
if ~isfield(options,'boundary_type'),options.boundary_type='dirichlet';end
options.boundary_type=validatestring(options.boundary_type,{'dirichlet','robin'});
if ~isfield(options,'output_directory')
    parent=fullfile(root,'output');if ~isfolder(parent),mkdir(parent);end
    options.output_directory=tempname(parent);
end
options.output_directory=char(options.output_directory);
assert(~isfolder(options.output_directory) && ~isfile(options.output_directory), ...
    'Output exists; choose a new output_directory.');
if isfield(options,'end_seconds'),validateattributes(options.end_seconds,{'numeric'},{'real','finite','scalar','positive'});end
if isfield(options,'dt_limit'),validateattributes(options.dt_limit,{'numeric'},{'real','finite','scalar','positive'});end
previous_path=path;restore_path=onCleanup(@()path(previous_path)); %#ok<NASGU>
addpath(fullfile(root,'src'));
if strcmp(case_name,'cycles_111')
    cfg=jsondecode(fileread(fullfile(root,'THMC_star.json')));
    names=fieldnames(options);for i=1:numel(names),cfg.(names{i})=options.(names{i});end
    assert(strcmp(char(cfg.case_code),'111'),'Use the named manuscript cases to change physical processes.');
    result=run_damaged_case(cfg);
else
    all_options=jsondecode(fileread(fullfile(root,'manuscript_cases.json')));
    cfg=all_options.(case_name);
    allowed_options={'boundary_type','output_directory','end_seconds','dt_limit','max_steps', ...
        'boundary_layer_thickness','robin','field_interval','stop_after_first_event', ...
        'event_max_dt','max_fields'};
    names=fieldnames(options);
    for i=1:numel(names)
        assert(ismember(names{i},allowed_options),'Unsupported archived-case option: %s',names{i});
        cfg.(names{i})=options.(names{i});
    end
    if ~isfield(cfg,'max_steps'),cfg.max_steps=Inf;end
    result=feval(['run_' case_name],cfg);
    if ~isfield(result,'completed')
        saved=load(fullfile(options.output_directory,'final_state.mat'),'completed');
        result.completed=saved.completed || result.first_event_complete;
    end
end
result.case_name=case_name;result.boundary_type=options.boundary_type;
result.output_directory=options.output_directory;result.version=thmc_version();
save(fullfile(options.output_directory,'release_run.mat'),'result','options');
fprintf('Earthquake Cycle THMC %s: %s, %s boundaries.\n',result.version,case_name,result.boundary_type);
end
