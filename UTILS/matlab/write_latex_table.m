function [] = write_latex_table(elat,elon,dep,M,eid,Nstn)
%By default it will write origin time (not event id) - see iwriteotime
% REQUIRED UPDATES:
% UPDATE for other format of table format (if MT needs to be dispalyed)
% UPDATE to use vargin 


[gamma,delta,M0,strike,dip,rake] = CMT2TT(M,0);
Mw = m02mw(1,M0);

iwriteotime = 1;

otime = eid2otime(eid); 

if iwriteotime
    for ii=1:length(elat)
        evnt_id = datestr(otime(ii),31); % use origin time instead of event id
        if nargin==6
            fprintf('%d & %s & %3.2f & %3.2f & %3.0f & %3.0f & %3.0f & %2.1f & %3.1f & %d \\\\ \n',...
                ii,evnt_id,elat(ii),elon(ii),strike(ii),dip(ii),rake(ii),Mw(ii),dep(ii),Nstn(ii));
        else
            fprintf('%d & %s & %3.2f & %3.2f & %3.0f & %3.0f & %3.0f & %2.1f & %3.1f \\\\ \n',...
                ii,evnt_id,elat(ii),elon(ii),strike(ii),dip(ii),rake(ii),Mw(ii),dep(ii));
            if 0
                fprintf('%s %3.2f  %3.2f  %3.0f  %3.0f  %3.0f  %2.1f  %3.1f \n',...
                    evnt_id{1},elat(ii),elon(ii),strike(ii),dip(ii),rake(ii),Mw(ii),dep(ii))
            end
        end
    end
    
else
    for ii=1:length(elat)
        evnt_id = otime2eid(otime(ii)); % use event id
        if nargin==6
            fprintf('%d & %s & %3.2f & %3.2f & %3.0f & %3.0f & %3.0f & %2.1f & %3.1f & %d \\\\ \n',...
                ii,evnt_id{1},elat(ii),elon(ii),strike(ii),dip(ii),rake(ii),Mw(ii),dep(ii),Nstn(ii));
        else
            fprintf('%d & %s & %3.2f & %3.2f & %3.0f & %3.0f & %3.0f & %2.1f & %3.1f \\\\ \n',...
                ii,evnt_id{1},elat(ii),elon(ii),strike(ii),dip(ii),rake(ii),Mw(ii),dep(ii));
            if 0
                fprintf('%s %3.2f  %3.2f  %3.0f  %3.0f  %3.0f  %2.1f  %3.1f \n',...
                    evnt_id{1},elat(ii),elon(ii),strike(ii),dip(ii),rake(ii),Mw(ii),dep(ii))
            end
        end
    end
end
