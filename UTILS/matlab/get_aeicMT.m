% Returns the moment tensor computed by AEIC using first P-motion polarity
%
% Note : Works for southern Alaska during MOOS array. For other, change
% ax3 and oran. 
%
% Vipul Silwal
% Oct 21, 2013

function [Mid] = get_aeicMT(event_id)

iregion = 1;
switch iregion
    case 1
        % SCAK
        axmoos = [-154 -146 58 62.5];
        ax3 = [axmoos -10 700];
        Mwran = [4 10];
    case 2
        % MFSZ
        axmoos = [-151 -147.5 63.5 65.5];
        ax3 = [axmoos -10 700];
        Mwran = [2 10];
end
oran=[];
[otime,slat,slon,sdep,M,M0,Mw,eid] = read_mech_AEICfp(oran,ax3,Mwran);
%[~,isort] = sort(eid);
[~,isort] = sort(sdep,'ascend');
%[~,isort] = sort(Mw,'descend');
display_eq_list(isort,otime,slon,slat,sdep,Mw,eid);

x=find(strcmp(eid(:),event_id)==1);
Mid = M(:,x)

if 1==0
    event_id = '20130714091656927';
    [Mid]= get_aeicMT(event_id);
end

