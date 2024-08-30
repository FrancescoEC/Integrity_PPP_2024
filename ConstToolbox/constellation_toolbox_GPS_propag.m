function [x_eci,v_eci,TIME,n_sat] = constellation_toolbox_GPS_propag(start_time,final_time,dt)
%function [x_eci,v_eci] = constellation_toolbox_GPS_propag(start_time,final_time,dt)
% Si tenga presente che in questa funzione vengono usate le routine della sola 
%constellationtoolbox per cui trasformazioni come ECI-ECEF possono essere
%poco accurate. E' stata inserita per la rapidità della lettura dello YUMA.
% Il TIME va rivisto poichè è solo il vettore dei tempi della simulazione
% comeune per tutti i satelliti

%addpath(genpath('CONSTELLATION TOOLBOX'));
%addpath(genpath('Propagatore_anani_per_simulatore'));
%addpath(genpath('/LIBRARY/TIME'));

%Oss ho messo 2010 e mi ha detto che il tempo di simulazione era tropo
%lontano dai dati del database e quindi mi chiamava un warning
%corrispondente ad una previsione troppo distante nel tempo.

[GPS_week_start, GPS_sec_start, GPS_day_start, rollover_flag_start] = utc2gps(start_time);
[GPS_week_final, GPS_sec_final, GPS_day_final, rollover_flag_final] = utc2gps(final_time);

[gps_alm_file, glo_alm_file]=find_alm(GPS_week_start);
[gps_alm] = readyuma(gps_alm_file);
[gps_ephem] = alm2geph(gps_alm);
[t_out, prn, x, v] = propgeph(gps_ephem ,[GPS_week_start,GPS_sec_start] , [GPS_week_final,GPS_sec_final],dt);

%Per ottenere le misure dei singoli satelliti è necessario fare un reshape
%basato sul numero di satelliti contenuto dall'almanacco parri alla
%dimensione delle righe del messaggio n. Ossia ogni n valori passiamo alle
%posizioni e velocità del dt successivo.
n_sat=length(gps_ephem(:,1));
n=n_sat;
m=length(x)/n;
X_1_ecef=reshape(x(:,1),n,m);% riga per riga abbiamo la storia del singolo satellite
X_2_ecef=reshape(x(:,2),n,m);
X_3_ecef=reshape(x(:,3),n,m);

V_1_ecef=reshape(v(:,1),n,m);% riga per riga abbiamo la storia del singolo satellite
V_2_ecef=reshape(v(:,2),n,m);
V_3_ecef=reshape(v(:,3),n,m);

t_out_org_1=reshape(t_out(:,1),n,m);
t_out_org_2=reshape(t_out(:,2),n,m);

%return

x_eci=cell(n,1);
v_eci=cell(n,1);

TIME=0:dt:((m-1)*dt);


for i=1:n
%Convertiamo a questo punto le misure dal riferimento ECEF ad ECI
       [x_eci_app, v_eci_app] = ecef2eci([t_out_org_1(i,:)',t_out_org_2(i,:)'], [X_1_ecef(i,:)',X_2_ecef(i,:)',X_3_ecef(i,:)'], [V_1_ecef(i,:)',V_2_ecef(i,:)',V_3_ecef(i,:)']); 
        x_eci{i}=x_eci_app';
        v_eci{i}=v_eci_app';
       
%questa è quella di default del costellation è sbagliata e accetta
%direttamente le matrici di pos,vel, e tempi di simulazione in settimana
%gps       
            
end