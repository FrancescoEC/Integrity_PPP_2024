function [sv,stat]=satpos(time,obs,nav,ephopt,sv)

global glc
% persistent Index_SAT
% 
% if isempty Index_SAT
%     Index_SAT
% end

teph=obs.time; sat=obs.sat;

switch ephopt
    case glc.EPHOPT_BRDC
        %% MODMENZ
        [sv,stat,tk]=ephpos(time,teph,sat,nav,-1,sv); return;
    case glc.EPHOPT_PREC
        [sv,stat]=peph2pos(time,sat,nav,1,sv); return;
    otherwise
        stat=0;return;
end

%% ModMENZ (ADD ERROR)

