function [L,P,Lc,Pc]=corr_meas(rtk,obs,nav,dantr,dants,phw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR test not surpport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global glc

k=2;%% Mod

lam=nav.lam(obs.sat,:);
L=zeros(glc.NFREQ,1);P=zeros(glc.NFREQ,1);Lc=0;Pc=0;

for i=1:glc.NFREQ
    L(i)=0; P(i)=0;
    if lam(i)==0||obs.L(i)==0||obs.P(i)==0,continue;end
    
    %antenna phase center and phase windup correction
    L(i)=obs.L(i)*lam(i)-dants(i)-dantr(i)-phw*lam(i);
    P(i)=obs.P(i)       -dants(i)-dantr(i);
    
end

% DCB correction 
[cbias,~]=getdcb(nav,obs,rtk.opt);
for i=1:glc.NFREQ
    if P(i)~=0,P(i)=P(i)-cbias(i);end
end
C1= lam(2)^2/(lam(2)^2-lam(1)^2);
C2=-lam(1)^2/(lam(2)^2-lam(1)^2);

%IFLC measurements
sat=obs.sat;
[sys,~]=satsys(sat);
%     if sys==glc.SYS_GAL,k=3;else,k=2;end
% if glc.NFREQ>=3&&(sys==glc.SYS_GPS),k=3;end
% if glc.NFREQ>=3&&(sys==glc.SYS_GAL),k=3;end
if L(1)~=0&&L(k)~=0,Lc=C1*L(1)+C2*L(k);end
if P(1)~=0&&P(k)~=0,Pc=C1*P(1)+C2*P(k);end

return


