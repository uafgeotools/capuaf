% script to check whether all green's function files are present

% INPUT PARAMETERS
smodel = 'scak';                    % structural model
fkdir = '/store/wf/FK_synthetics/'; % path to green's function
% Also check whether the directory exists or not
% ls /store/wf/FK_synthetics/smodel (replace $smodel)

mindep = 0;     % minimum depth
maxdep = 200;   % maximum depth
depinc = 1;     % depth increment
mindis = 1;     % minimum epicentral distance
maxdis = 500;  % maximum epicentral distance

ddir = strcat(fkdir,smodel,'/');
D = dir(ddir);

mod = dlmread(strcat(ddir,smodel));
emptydir = cumsum([0;mod(:,1)]);  % layer interface = no green's function generated

% OUTPUT log file
log_file = strcat('./greens_checklsit_',smodel);
fid = fopen(log_file,'w');

% File existense log is generated at $log_file in the current directory
for ii=mindep:depinc:maxdep
    depdir = strcat(ddir,smodel,'_',num2str(ii));
    D1 = dir(depdir);
    if length(D1)==12*(maxdis-mindis+1)+2
        disp(sprintf('%s_%d is complete',smodel,ii));
        fprintf(fid,'%s_%d is complete\n',smodel,ii);
    elseif max(eq(ii,emptydir))
        disp(sprintf('%s_%d is at layer interface',smodel,ii));
        fprintf(fid,'%s_%d is at layer interface\n',smodel,ii);
    else
        disp(sprintf('%s_%d is incomplete',smodel,ii));
        fprintf(fid,'%s_%d is incomplete\n',smodel,ii);
    end
end
    
% FUTURE WORK:  ALSO TELL WHAT DEPTH and DISTANCE's green's function are
% absent (time-issue...takes too long to surf through ALL filenames)