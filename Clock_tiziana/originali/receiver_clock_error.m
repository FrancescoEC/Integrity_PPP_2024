clc
clear all
close all

global LIGHT_SPEED

LIGHT_SPEED=2.99792458*(10^8);


Ts=0.1;
Tf=10000;
time=0:Ts:Tf;

flag_clock=1;           % flag_clock=1 --> US-OCXO
                        % flag_clock=2 --> Cesium High Performance
if(flag_clock==1)
x0=0.0067817;
y0=8.617/LIGHT_SPEED;
d0=0;
else
x0=0.1935958862e-3;
y0=0.3637978807e-11;
d0=0;
end 


A=[1,Ts;0,1];
B=[Ts,Ts^2;0,Ts];

%Calcolo del contributo costante (forzante)

U_cost=[0.5*d0*(Ts^2);d0*Ts]*ones(1,length(time));

% Calcolo del contributo randomico (forzante)

[x,y]=clock_noise(Ts,Tf,flag_clock);
%[x,y]=clock_noise_tot(Ts,Tf,flag_clock);

wnc=U_cost;
wnc=wnc';
IC=[x0;y0];

yy = ltitr(A,B,wnc,IC);
delta_t_rec=yy(:,1)+x;
offset_f_rec=yy(:,2)+y;

clk_bias = delta_t_rec*LIGHT_SPEED;
clk_drift = offset_f_rec*LIGHT_SPEED;

figure
plot(time,delta_t_rec);grid;
title ('time error');
xlabel('time (s)');
ylabel(' time error (s) ');

figure
plot(time,offset_f_rec);grid;
title ('frequency error');
xlabel('time (s)');
ylabel(' frequency error ');

% da confrontare con la figura a pag. 179 tesi Ambrosini
% figure
% plot(time,clk_bias);grid;
% title ('clock bias');
% xlabel('time (s)');
% ylabel(' clock bias (m) ');
% 
% figure
% plot(time,clk_drift);grid;
% title ('clock drift');
% xlabel('time (s)');
% ylabel(' clock drift (m/s) ');

