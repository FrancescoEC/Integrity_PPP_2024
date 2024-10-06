function [x,y]=clock_noise(Ts,Tf,flag_clock,seed_cbias)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                            
%    Phase and Frequency noise computation (Random walk FM + Flicker FM + White noise FM)
%    
%    INPUT: Ts               --> sampling time;
%           Tf               --> time duration;
%           flag_clock       --> tipe of clock : 1 = US-OCXO
%                                                2 = Cesium
%
%    OUTPUT: x               --> clock_time_noise
%            y               --> clock_frequency_noise
%
%    Ref: "Relativity in the Global Positioning System", Neil Ashby, 2003
%         "Precise Time and Time Interval Clocks, Time Frames and
%         Frequency", James R. Clynch, 2003
%         "Techniques for Frequency Stability Analysis", W. J. Riley,
%         Symmetricom, 2003
% 
%    Author : Tiziana Ingrassia 1/08/2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time=0:Ts:Tf;
z=0.7;

if (flag_clock==1)       %US-OCXO
    
% RANDOM WALK FM (PSD is 40 db/dec, so filter must be 20 db/dec with double threshold frequency)

FF=logspace(-4,1,5000); 
f0=1/86400;
f1=0.1;
t0=1/(2*pi*f0);
t1=1/(2*pi*f1);
num=[t1 1]*t0/t1;
den=[t0 1];

[a,b,c,d]=tf2ss(num,den);                 % transfer function to state-space conversion
Filter1_ss=ss(a,b,c,d);                   % creates an object SYS representing the continuous-time state-space model

Filter1_d=c2d(Filter1_ss,Ts,'tustin');    % converts continuous-time dynamic system to discrete time
[az,bz,cz,dz]=ssdata(Filter1_d);          % returns the A,B,C,D matrices of the discrete state-space model FL_SYSD

TF_Rwfm=freqs(num,den,2*pi*FF);           % TF_RWFM  Laplace-transform (s-domain) frequency response

% FLICKER (PSD is 20 db/dec, so filter must be 10 db/dec with single threshold frequency)

G=2e-25;                                % Allan deviation@1sec US-OCXO Tab. Oscillator Options - Summary Performance data (p.1) 
                                          % (Ref:."EndRun Technologies "Disciplined Oscillator Options for GPS-Synchronized Time&Frequency Standards"")

[FL_num,FL_den]=flicker_tf(G,1,1);        % transfer function of the flicker noise 

[a2,b2,c2,d2]=tf2ss(FL_num,FL_den);       % transfer function to state-space conversion
Filter2_ss=ss(a2,b2,c2,d2);               % creates an object SYS representing the continuous-time state-space model

Filter2_d=c2d(Filter2_ss,Ts,'tustin');    % converts continuous-time dynamic system to discrete time
[afz,bfz,cfz,dfz]=ssdata(Filter2_d);      % returns the A,B,C,D matrices of the discrete state-space model FL_SYSD

TF_flicker=freqs(FL_num,FL_den,2*pi*FF);

% Integrator
num_int=[0 1];
den_int=[1 0];
int_tf=freqs(num_int,den_int,2*pi*FF);

[a3,b3,c3,d3]=tf2ss(num_int,den_int);
Filter3_ss=ss(a3,b3,c3,d3);

else if (flag_clock==2)
    
% RANDOM WALK FM : Filtro 40 db/dec con doppia freq di taglio
FF=logspace(-3,2,5000); 

% Cesium
f0=1e-9;
f1=1e-7;

t0=1/(2*pi*f0);
t1=1/(2*pi*f1);
num=[t1 1]*t0/t1;
den=[t0 1];

[a,b,c,d]=tf2ss(num,den);
Filter1_ss=ss(a,b,c,d);

Filter1_d=c2d(Filter1_ss,Ts,'tustin');    % converts continuous-time dynamic system to discrete time
[az,bz,cz,dz]=ssdata(Filter1_d);          % returns the A,B,C,D matrices of the discrete state-space model FL_SYSD

TF_Rwfm=freqs(num,den,2*pi*FF);           % RWFM

% FLICKER 
% Cesium
G=5e-23;                                  % Allan deviation@1sec Cesium Clocks- High Performance Tab.1: Specifications
                                          % (Ref:."Cesium Clock and Hydrogen Maser Compared", Application note, Oscilloquartz
[FL_num,FL_den]=flicker_tf(G,1e-3,4);     % transfer function of the flicker noise 

[a2,b2,c2,d2]=tf2ss(FL_num,FL_den);
Filter2_ss=ss(a2,b2,c2,d2);

Filter2_d=c2d(Filter2_ss,Ts,'tustin');    % converts continuous-time dynamic system to discrete time
[afz,bfz,cfz,dfz]=ssdata(Filter2_d);      % returns the A,B,C,D matrices of the discrete state-space model FL_SYSD

TF_flicker=freqs(FL_num,FL_den,2*pi*FF);

figure
loglog(FF,(abs(TF_flicker).^2));grid

% Integrator
num_int=[0 1];
den_int=[1 0];
int_tf=freqs(num_int,den_int,2*pi*FF);

[a3,b3,c3,d3]=tf2ss(num_int,den_int);
Filter3_ss=ss(a3,b3,c3,d3);
    
    else
        
        FF=logspace(-3,2,5000);

        % Rubidium
        f0=1e-7;
        f1=1e-4;

        t0=1/(2*pi*f0);
        t1=1/(2*pi*f1);
        num=[t1 1]*t0/t1;
        den=[t0 1];

        [a,b,c,d]=tf2ss(num,den);
        Filter1_ss=ss(a,b,c,d);

        Filter1_d=c2d(Filter1_ss,Ts,'tustin');    % converts continuous-time dynamic system to discrete time
        [az,bz,cz,dz]=ssdata(Filter1_d);          % returns the A,B,C,D matrices of the discrete state-space model FL_SYSD

        TF_Rwfm=freqs(num,den,2*pi*FF);           % RWFM

        % FLICKER 
        % Rubidium
        G=2e-22;                                 % Allan deviation@1sec US-OCXO Tab. Oscillator Options - Summary Performance data (p.1) 
                                                 % (Ref:."EndRun Technologies "Disciplined Oscillator Options for GPS-Synchronized Time&Frequency Standards"")
        
        [FL_num,FL_den]=flicker_tf(G,1e-2,2);     % transfer function of the flicker noise 

        [a2,b2,c2,d2]=tf2ss(FL_num,FL_den);
        Filter2_ss=ss(a2,b2,c2,d2);

        Filter2_d=c2d(Filter2_ss,Ts,'tustin');    % converts continuous-time dynamic system to discrete time
        [afz,bfz,cfz,dfz]=ssdata(Filter2_d);      % returns the A,B,C,D matrices of the discrete state-space model FL_SYSD

        TF_flicker=freqs(FL_num,FL_den,2*pi*FF);

        figure
        loglog(FF,(abs(TF_flicker).^2));grid

        % Integrator
        num_int=[0 1];
        den_int=[1 0];
        int_tf=freqs(num_int,den_int,2*pi*FF);

        [a3,b3,c3,d3]=tf2ss(num_int,den_int);
        Filter3_ss=ss(a3,b3,c3,d3);

    end
end

% run simulink 
options=simset('SrcWorkspace','current');
sim('prova_simu_clock',[],options);        % tf


end

