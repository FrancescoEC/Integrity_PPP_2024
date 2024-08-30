function [position_ecef,velocity_ecef,tx,dt_sv,dt_dot_sv]=ComputeBRDC(eph_sel,tx_raw,system_type)
%
% orbital_propagator_r implements the updating of GPS ephemerides and the transformation in ECEF frame 
% 
%in input:
%         - an ephemerides matrix, structured as defined in RINEX format
%           eph (can be 1 eph set for each SV or more)
%         - a vector containing prn of satellites of interest sv (in increasing order)
%         - vector of epochs of the desired orbit position (in which the ECEF coordinates are computed) toss
%
%in output:
%         - satellite positions in ECEF frame at epoch toss (meters) position_ECEF
%         - satellite velocity in ECEF frame at epoch toss (m/sec) velocity ECEF


%light speed [m/sec]
v_light = 299792458;
% costante di gravitazione geocentrica [m^3/sec^2]
if system_type == 'G'
    mu = 3.986005*10^14;
    rot_e = 7.2921151467e-5;
elseif system_type == 'E'
    mu = 3.986004418*10^14;
    rot_e = 7.2921151467e-5;
elseif system_type == 'C'
    mu = 3.986004418*10^14;
    rot_e = 7.2921150e-5;
end
%number of satellites
k=size(eph_sel,2);


%prn=eph_sel(1,:);
M0=eph_sel(10,:)';
roota=eph_sel(14,:)';
Dn=eph_sel(9,:)';
ecc=eph_sel(12,:)';
omega=eph_sel(21,:)';
Cuc=eph_sel(11,:)';
Cus=eph_sel(13,:)';
Crc=eph_sel(20,:)';
Crs=eph_sel(8,:)';
i0=eph_sel(19,:)';
i_dot0=eph_sel(23,:)';
Cic=eph_sel(16,:)';
Cis=eph_sel(18,:)';
OMEGA0=eph_sel(17,:)';
OMEGA_dot=eph_sel(22,:)';
toe=eph_sel(15,:)';
toc=eph_sel(3,:)';
Af0=eph_sel(4,:)';
Af1=eph_sel(5,:)';
Af2=eph_sel(6,:)';
tgd=eph_sel(29,:)';


%satellite clock correction
[dt_sv,dt_dot_sv]=sv_clock_corr(Af0,Af1,Af2,toc,tx_raw);
% dt_sv=dt_sv-tgd;
%approximate relativistic correction
dtr = (-4.442807633e-10)*ecc.*roota.*sin(M0);

if system_type == 'C'
    %approximate relativistic correction
    dtr = -2*sqrt(mu)/v_light^2*ecc.*roota.*sin(M0);
end

%epoch of transmission
tx = tx_raw-dt_sv-dtr;

% semi-major axis
A=roota.^2;

% calcolo del moto medio
n0=sqrt(mu./A.^3);

% intervallo tra epoca di interesse tx ed epoca di riferimento delle effemeridi
dt=tx-toe;
dt=check_t(dt);

% aggiornamento del moto medio all'epoca toss
n=n0+Dn;

% aggiornamento dell'anomalia media all'epoca toss
M=M0+n.*dt;

% derivative of Mean Anomaly M
M_dot=n;

AnEcc=nan(k,1);
for w=1:k
    % calcolo dell'anomalia eccentrica all'epoca toss
    E(1)=M(w);
    E(2)=E(1)-((E(1)-ecc(w)*sin(E(1))-M(w))/(1-ecc(w)*cos(E(1))));
    bb=1 ;
    while abs(E(bb+1)-E(bb))>10^-8 %% da -7 a -8
        bb=bb+1;
        E(bb+1)=E(bb)-((E(bb)-ecc(w)*sin(E(bb))-M(w))/(1-ecc(w)*cos(E(bb))));
    end
    AnEcc(w)=E(bb+1);
end

% derivative of Eccentric Anomaly AnEcc
AnEcc_dot=M_dot./(1-ecc.*cos(AnEcc));
        
% calcolo dell'anomalia vera all'epoca toss
%AnVera=2*atan(sqrt((1+ecc)./(1-ecc)).*tan(AnEcc/2));
AnVera=atan2((sqrt(1-ecc.^2)).*sin(AnEcc)./(1-ecc.*cos(AnEcc)),(cos(AnEcc)-ecc)./(1-ecc.*cos(AnEcc)));
    
% derivative of True Anomaly AnVera
AnVera_dot=sin(AnEcc).*AnEcc_dot.*(1+ecc.*cos(AnVera))./((1-ecc.*cos(AnEcc)).*sin(AnVera));

% calcolo dell'argomento di latitudine all'epoca toss
u=omega+AnVera;

% argment of latitude correction
du=Cuc.*cos(2*u)+Cus.*sin(2*u);
% radius correction
dr=Crc.*cos(2*u)+Crs.*sin(2*u);
% inclination correction
di=Cic.*cos(2*u)+Cis.*sin(2*u);

% derivative of argument of latitude
u_dot=AnVera_dot+2*(Cus.*cos(2*u)-Cuc.*sin(2*u)).*AnVera_dot;
         
% calcolo del raggio vettore all'epoca toss
r=A.*(1-ecc.*(cos(AnEcc)))+dr;
    
% derivative of Range earth center - satellite
r_dot=A.*ecc.*sin(AnEcc).*AnEcc_dot+2*(Crs.*cos(2*u)-Crc.*sin(2*u)).*AnVera_dot;
         
% aggiornamento dell'inclinazione all'epoca toss
i=i0+i_dot0.*dt+di;
    
% derivative of orbit inclination
i_dot=i_dot0+2*(Cis.*cos(2*u)-Cic.*sin(2*u)).*AnVera_dot;
         
% calcolo della longitudine del nodo ascendente all'epoca toss
LAMBDA=OMEGA0+(OMEGA_dot-rot_e).*dt-rot_e*toe;
    
% derivative of Longitude of Ascending Node
LAMBDA_dot=OMEGA_dot-rot_e;    
                                                                                  
% corrected argment of latitude
u=u+du;

% coordinate dei satelliti nel sistema di riferimento orbitale
Xorb=r.*cos(u) ;
Yorb=r.*sin(u) ;
Zorb=zeros(size(Xorb)) ;

% derivative of satellite orbital coordinates
Xorb_dot=r_dot.*cos(u)-Yorb.*u_dot;
Yorb_dot=r_dot.*sin(u)+Xorb.*u_dot;

% coordinate dei satelliti nel sistema di coordinate ECEF 
Xecef=nan(k,1); Yecef=nan(k,1); Zecef=nan(k,1);
for dd=1:k
    Rx=matrix_rot(i(dd),'O','x');
    Rz=matrix_rot(LAMBDA(dd),'O','z');
    Xecef(dd)=Rz(1,:)*Rx*[Xorb(dd);Yorb(dd);Zorb(dd)];
    Yecef(dd)=Rz(2,:)*Rx*[Xorb(dd);Yorb(dd);Zorb(dd)];
    Zecef(dd)=Rz(3,:)*Rx*[Xorb(dd);Yorb(dd);Zorb(dd)];    
end

% derivative of satellite ECEF coordinates
Xecef_dot=Xorb_dot.*cos(LAMBDA)-Yorb_dot.*cos(i).*sin(LAMBDA)+Yorb.*sin(i).*sin(LAMBDA).*i_dot-Yecef.*LAMBDA_dot;
Yecef_dot=Xorb_dot.*sin(LAMBDA)+Yorb_dot.*cos(i).*cos(LAMBDA)-Yorb.*sin(i).*cos(LAMBDA).*i_dot+Xecef.*LAMBDA_dot;
Zecef_dot=Yorb_dot.*sin(i)+Yorb.*cos(i).*i_dot;

position_ecef=[Xecef'; Yecef'; Zecef'];
velocity_ecef=[Xecef_dot'; Yecef_dot'; Zecef_dot'];

%relativistic effect correction
dtr = (-4.442807633e-10)*ecc.*roota.*sin(AnEcc);

%
dt_sv = dt_sv + dtr;
