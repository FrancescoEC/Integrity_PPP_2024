%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [num,den]=flicker_tf(const,fcs,ndec);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Computation of the approximate filter such that the PSD
%
%            PSD(f) = const*(fcs/f) for f < fcs and PSD(f) = const for f >= fcs
%
%  is shaped through [num(s),den(s)] in a ndec decades bandwidth
%  below the fcs frequency.
%
%  Based on NBS Technical Note 604, by J.A. Barnes and Stephen Jarvis, Jr
%  
%
%            Usage: [num,den]=flicker_tf(const,fcs,ndec);
%
%             Input :   
%                           const = PSD plateau value (high frequency flat component)
%                             fcs = corner frequency of 1/f behaviour
%                            ndec = number of frequency decades for the 1/f behaviour
%
%            Output : [num,den] of shaping filter
%
%
%
% Ver 1.0                 Jan 18 2002                   G.C.
%
%=======================================================================
% 
function [num,den]=flicker_tf(const,fcs,ndec); 
%
%
f0=fcs;
%
Aflat=sqrt(const);
%
nn=ndec;
%
%
%
beta      = 1/9;
%
%tau2      = 1/(2*pi*f0)/beta^(nn-1);
tau2      = 1/(pi*f0)/beta^(nn-1);
%
tau1ptau2 = 3*tau2;
%
%
betan=zeros(nn,1);
%
num=[0 1];
den=[0 1];
%
for i=1:nn
   betan(i) = beta^(i-1);
   num = conv([tau2*betan(i) 1],num);
   den = conv([tau1ptau2*betan(i) 1],den);
end;
num=num*Aflat*den(2)/num(2);
%
%
%
%=======================================================================
