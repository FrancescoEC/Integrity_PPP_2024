function [tsv_corr,dtsv_corr]=sv_clock_corr(Af0,Af1,Af2,toc,t_raw)

%sv_clock_corr computes satellite clock error, using the parameters
%broadcasted in the Navigation Message
%
%in input:
%         correction parameters Af0,Af1,Af2
%         corrections reference epoch toc
%         raw epoch to correct t_raw
%in output:
%         subtractive correction to t_raw  

dt = check_t(t_raw-toc);
tsv_corr=Af0+Af1.*(dt)+Af2.*(dt).^2;
dtsv_corr=Af1+2*Af2.*(dt);