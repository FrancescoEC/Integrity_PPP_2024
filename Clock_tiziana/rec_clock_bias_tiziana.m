function [delta_t_rec,offset_f_rec,clk_bias,clk_drift,time]=rec_clock_bias_tiziana(x0,y0,d0,delta,nsim,seed_c)

global glc

Ts=delta;
Tf=(nsim-1)*delta;       %Tf=num_iterazioni-1
time=0:Ts:Tf;
flag_clock=1;

A=[1,Ts;0,1];
B=[Ts,Ts^2;0,Ts];

%Calcolo del contributo costante (forzante)

U_cost=[0.5*d0*(Ts^2);d0*Ts]*ones(1,length(time));

% Calcolo del contributo randomico (forzante)

[x,y]=clock_noise(Ts,Tf,flag_clock,seed_c);

wnc=U_cost;
wnc=wnc';
IC=[x0;y0];

yy = ltitr(A,B,wnc,IC);
delta_t_rec=yy(:,1)+x;
offset_f_rec=yy(:,2)+y;

clk_bias = delta_t_rec*glc.CLIGHT;
clk_drift = offset_f_rec*glc.CLIGHT;

end