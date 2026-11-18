function report=inelast_stream_fields(action,filename,row,index)
% Preserve every legacy exported field in bounded buffers and v7.3 MAT storage.
% Numerical state and scalar export cadence are unchanged.
persistent state
if nargin<3,row=struct();end
if nargin<4,index=[];end
switch action
 case 'append'
  if isempty(state)
   assert(~isfile(filename),'Field output already exists.');
   names=fieldnames(row);capacity=4096;buffer=struct();
   for j=1:numel(names),buffer.(names{j})=zeros(capacity,numel(row.(names{j})));end
   stream_schema=struct('format','every original exported row, double precision','fields',{names},'buffer_capacity',capacity); %#ok<NASGU>
   save(filename,'stream_schema','-v7.3');
   state=struct('filename',filename,'file',matfile(filename,'Writable',true),'names',{names}, ...
    'capacity',capacity,'buffer',buffer,'count',0,'written',0);
  end
  assert(strcmp(state.filename,filename),'One field stream per process.');
  assert(index==state.written+state.count+1,'Nonconsecutive field export.');
  state.count=state.count+1;
  for j=1:numel(state.names)
   name=state.names{j};state.buffer.(name)(state.count,:)=row.(name)(:).';
  end
  if state.count==state.capacity,state=flush(state);end
  report=struct('rows',state.written+state.count);
 case 'finish'
  if isempty(state),report=struct('rows',0);return,end
  assert(strcmp(state.filename,filename));state=flush(state);
  report=struct('rows',state.written,'fields',{state.names});state=[];
 otherwise
  error('Unsupported field stream action.');
end
end
function state=flush(state)
if state.count==0,return,end
rows=state.written+(1:state.count);
for j=1:numel(state.names)
 name=state.names{j};state.file.(name)(rows,1:size(state.buffer.(name),2))=state.buffer.(name)(1:state.count,:);
end
state.written=state.written+state.count;state.count=0;
end
