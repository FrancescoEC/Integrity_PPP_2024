clear all
close all
clc

h_sat=1200000;
minElevAngle=5;
GAL_input= readyuma('GAL_Alm_w220_yuma_defaults_36SVs.yuma');
GPS_input= readyuma('GPS_Alm_w220_yuma_defaults_32SVs_1.yuma');

File_name_GPS='LEO4_GPS.yuma';
File_name_GAL='Gal_ref.yuma';
Traslazione_GAL_Matteo_PRN=0;

Window=120;
StopIndex=100;
%% Fprint
%ID_Free={'PRN-01','PRN-02','PRN-03','PRN-04'};

% ID_Free_n=10;
%ID_Gal=[1 2 3 11 12 18 20 27 29 31 32];
ID_Gal=[];
%ID_GPS=[4 5 7 8 9 14 16 20 27 30];
ID_GPS=[8,9,20,27];

Index_tot=1:9;
subset=3;
% Week=220;

mission.StartDate = datetime(2023,28,06,0,0,0);
mission.Duration = hours(6);
mission.scenario = satelliteScenario(mission.StartDate, ...
    mission.StartDate + mission.Duration, 60);
Station_Position=repmat([45.8140,8.6121,0],2.5*3600,1);
figure
geoplot(Station_Position(:,1), Station_Position(:,2), "b-");
% geolimits([30 50],[-110 -50]);

%Station.obj = (mission.scenario,Station_Position, CoordinateFrame="geographic", Name="Station");
Station.obj=groundStation(mission.scenario,Station_Position(1,1), Station_Position(1,2),MinElevationAngle=minElevAngle);
MinElevationAngle=minElevAngle;
Station.obj.MarkerColor = "green";

%% Constellation - Multiorbit

% Inclination=[80,65,50,10];
% NumberOfSat=[42,81,49,48];
% Inclination=[65,50,10];
% NumberOfSat=[81,49,48];
% Inclination=[65,50];
% NumberOfSat=[81,49];

Inclination=[65,50,10];
NumberOfSat=[81,49,48];

% Inclination=[65,54];
% NumberOfSat=[81,49];

% Inclination=65;
% NumberOfSat=81;
% Planes=[6,9,7,8];
%Planes=[9,7,8];

% Essential
Planes=[9,7,8];

% Inclination=[65,54];
% NumberOfSat=[81,49];

% Inclination=65;
% NumberOfSat=81;
% Planes=[6,9,7,8];
%Planes=[9,7,8];

% Essential
% Planes=[3];

Appo={"Sat:in1","Sat:in2","Sat:in3","Sat:in4"}
Prec=0;

for ii=1:length(Inclination)

constellation(ii).Radius = 6378137+h_sat; % meters
constellation(ii).Inclination = Inclination(ii); % deg
constellation(ii).TotalSatellites =NumberOfSat(ii) ;
constellation(ii).PRN_range=Prec+1:NumberOfSat(ii)+Prec;
constellation(ii).GeometryPlanes = Planes(ii);
constellation(ii).Phasing = 1;
constellation(ii).ArgLat = 15; % deg
constellation(ii).obj = walkerDelta(mission.scenario, ...
constellation(ii).Radius, ...
constellation(ii).Inclination, ...
constellation(ii).TotalSatellites, ...
constellation(ii).GeometryPlanes, ...
constellation(ii).Phasing, ...
ArgumentOfLatitude=constellation(ii).ArgLat, ...
Name=Appo{ii});
Prec=Prec+NumberOfSat(ii);

% sensor(ii).HalfAngle = 70; % deg
% sensor(ii).Names = constellation(ii).obj.Name + " satellite";
% sensor(ii).obj = conicalSensor(constellation(ii).obj, MaxViewAngle=sensor(ii).HalfAngle*2, MountingLocation=[0 0 0], Name=sensor(ii).Names);
% sensor(ii).FOV.obj = fieldOfView(sensor(ii).obj);
% accessAnalysis(ii).obj = access(sensor(ii).obj, Station.obj);
accessAnalysis(ii).obj =access(constellation(ii).obj,Station.obj); 
[accessAnalysis(ii).Status, mission.TimeSteps] = accessStatus(accessAnalysis(ii).obj);
accessAnalysis(ii).Intervals = accessIntervals(accessAnalysis(ii).obj);
accessAnalysis(ii).Intervals = sortrows(accessAnalysis(ii).Intervals,"StartTime");
disp(accessAnalysis(ii).Intervals)

accessAnalysis(ii).SystemWideAccessStatus = any(accessAnalysis(ii).Status, 1);
figure
stairs(mission.TimeSteps, accessAnalysis(ii).SystemWideAccessStatus);
ylim([-.2, 1.2]);
xlabel("Time");
ylabel("System-Wide Access Status");
yticks([0,1]);

accessAnalysis(ii).SystemWideAccessDuration = sum(accessAnalysis(ii).SystemWideAccessStatus) * mission.scenario.SampleTime; % seconds
accessAnalysis(ii).SystemWideAccessPercentage = (accessAnalysis(ii).SystemWideAccessDuration/seconds(mission.Duration))*100;
disp("System Wide Access = " + accessAnalysis(ii).SystemWideAccessPercentage + " %")

accessAnalysis(ii).SystemWideAccessStatus = any(accessAnalysis(ii).Status, 1);
stairs(mission.TimeSteps, accessAnalysis(ii).SystemWideAccessStatus);
ylim([-.2, 1.2]);
xlabel("Time");
ylabel("System-Wide Access Status");
yticks([0,1]);

accessAnalysis(ii).SystemWideAccessDuration = sum(accessAnalysis(ii).SystemWideAccessStatus) * mission.scenario.SampleTime; % seconds
accessAnalysis(ii).SystemWideAccessPercentage = (accessAnalysis(ii).SystemWideAccessDuration/seconds(mission.Duration))*100;
disp("System Wide Access = " + accessAnalysis(ii).SystemWideAccessPercentage + " %")

hide(constellation(ii).obj.Orbit);
constellation(ii).obj.ShowLabel = false;
Station.obj.ShowLabel = false;

end

%% Skyplot
% LEO_PNT(1).Name='MultiInclination';
% LEO_PNT(2).Name='MultiInclination';
% LEO_PNT(3).Name='MultiInclination';
% LEO_PNT(4).Name='MultiInclination';

% queryTime=mission.StartDate;
% acStats = accessStatus(accessAnalysis(1).obj,queryTime);
% [az,el] = aer(Station.obj,constellation(1).obj(acStats),queryTime);
% skyplot(az,el,constellation(1).obj(acStats).Name)

%% All Constellation
LEO_PNT_T.azTotal=[];
LEO_PNT_T.elTotal=[];
LEO_PNT_T.acStatsTotal=[];
LEO_PNT_T.timeHistoryTotal=[];
LEO_PNT_T.NameTotal={};
for ii=1:length(Inclination)
LEO_PNT(ii).Component=constellation(ii).obj(1).Name;
[LEO_PNT(ii).acStatsAllTime,LEO_PNT(ii).timeHistory] = accessStatus(accessAnalysis(ii).obj);
[LEO_PNT(ii).azHistory,LEO_PNT(ii).elHistory] = aer(Station.obj,constellation(ii).obj);
LEO_PNT(ii).numTimeStamps = numel(LEO_PNT(ii).timeHistory);
LEO_PNT(ii).elHistory(LEO_PNT(ii).acStatsAllTime == 0) = NaN;
LEO_PNT(ii).azHistoryTranspose = LEO_PNT(ii).azHistory';
LEO_PNT(ii).elHistoryTranspose = LEO_PNT(ii).elHistory';

LEO_PNT_T.azTotal=[LEO_PNT_T.azTotal,LEO_PNT(ii).azHistoryTranspose];
LEO_PNT_T.elTotal=[LEO_PNT_T.elTotal,LEO_PNT(ii).elHistoryTranspose];
LEO_PNT_T.acStatsTotal=[LEO_PNT_T.acStatsTotal,LEO_PNT(ii).acStatsAllTime.'];
LEO_PNT_T.timeHistoryTotal=[LEO_PNT_T.timeHistoryTotal,LEO_PNT(ii).timeHistory.'];
LEO_PNT_T.NameTotal=[LEO_PNT_T.NameTotal,constellation(ii).obj(1:end).Name];
end

% Find the PRN index of each satellite
satNames = LEO_PNT_T.NameTotal;
prnIndex = 1:sum(NumberOfSat); double(string(satNames(:,5:end)));

% To better visualize each GPS satellite, scale the status with the PRN
% index
acStatsAllTime = double(LEO_PNT_T.acStatsTotal);
acStatsAllTime(acStatsAllTime == 0) = NaN;

% Plot the satellite visibility chart
figure
colors = colororder;
firstLineColor = colors(1,:);
plot(LEO_PNT_T.timeHistoryTotal(:,1),prnIndex.*acStatsAllTime, ...
    Color=firstLineColor,LineWidth=1)
xlim([LEO_PNT_T.timeHistoryTotal(1,1) LEO_PNT_T.timeHistoryTotal(end,1)])
ylim([min(prnIndex)-1 max(prnIndex)+1])
xlabel("Time")
ylabel("Satellite PRN Index")
title("Satellite Visibility Chart")
yticks(prnIndex)
grid on

%% Setting 

figure
plot(LEO_PNT_T.timeHistoryTotal(:,1),sum((~isnan(acStatsAllTime.'))), ...
    Color=firstLineColor,LineWidth=1)
xlim([LEO_PNT_T.timeHistoryTotal(1,1) LEO_PNT_T.timeHistoryTotal(end,1)])
xlabel("Time")
ylabel("Satellite PRN Index")
title("Satellite Visibility Chart")
yticks(prnIndex)
grid on


Index_tot=1:sum(NumberOfSat);
Appo=LEO_PNT_T.elTotal(1:Window,:);
Appo(isnan(Appo))=0;
% Set_Satelliite=Index_tot(any(acStatsAllTime));
Score=sum(Appo)./(1+sum(Appo~=0));
% Score=sum(Appo>20);
%Score=sum(Appo>0);
figure,plot(Score);
[a,b]=sort(Score,'descend');
Mask=LEO_PNT_T.elTotal(1:Window,:)>0; 
if StopIndex<sum(NumberOfSat)
Mask(:,b(StopIndex:end))=0; %b(1:40)
end
% Mask(isnan(Mask))=0;
Appo_2=Appo(1:Window,:).*Mask;
surf(Appo_2);

figure
colors = colororder;
firstLineColor = colors(1,:);
if StopIndex<sum(NumberOfSat)
plot(LEO_PNT_T.timeHistoryTotal(1:Window,1),Index_tot(b(1:StopIndex-1)).*acStatsAllTime(1:Window,b(1:StopIndex-1)), ...
    Color=firstLineColor,LineWidth=1)
else
plot(LEO_PNT_T.timeHistoryTotal(1:Window,1),Index_tot(b).*acStatsAllTime(1:Window,Index_tot(b)), ...
    Color=firstLineColor,LineWidth=1)

end
xlim([LEO_PNT_T.timeHistoryTotal(1,1) LEO_PNT_T.timeHistoryTotal(Window,1)])
ylim([min(prnIndex)-1 max(prnIndex)+1])
xlabel("Time")
ylabel("Satellite PRN Index")
title("Satellite Visibility Chart")
yticks(prnIndex)
grid on



% satNames = char(sat(:).Name');
% prnIndex = double(string(satNames(:,5:end)));
% 
% % To better visualize each GPS satellite, scale the status with the PRN
% % index
% acStatsAllTime = double(acStatsAllTime);
% acStatsAllTime(acStatsAllTime == 0) = NaN;
% 
% % Plot the satellite visibility chart
% colors = colororder;
% firstLineColor = colors(1,:);
% plot(timeHistory,prnIndex.*acStatsAllTime, ...
%     Color=firstLineColor,LineWidth=1)
% xlim([timeHistory(1) timeHistory(end)])
% ylim([min(prnIndex)-1 max(prnIndex)+1])
% xlabel("Time")
% ylabel("Satellite PRN Index")
% title("Satellite Visibility Chart")
% yticks(prnIndex)
% grid on
% 
% sp = skyplot([],[]);
% for tIdx = 1:LEO_PNT(ii).numTimeStamps
%     for ii=1:4
%     set(sp, ...
%         AzimuthData=LEO_PNT(ii).azHistoryTranspose(tIdx,:), ...
%         ElevationData=LEO_PNT(ii).elHistoryTranspose(tIdx,:), ...
%         LabelData=constellation(ii).obj.Name,ColorOrder=(ii*40/256)*ones(1,3));
%     % For slower animation, use 'drawnow' instead of 'drawnow limitrate'
%     drawnow 
%     pause(0.1)
%     end
% end
% 
% % Find the PRN index of each satellite
% satNames = char(sat(:).Name');
% prnIndex = double(string(satNames(:,5:end)));
% 
% % To better visualize each GPS satellite, scale the status with the PRN
% % index
% acStatsAllTime = double(acStatsAllTime);
% acStatsAllTime(acStatsAllTime == 0) = NaN;
% 
% % Plot the satellite visibility chart
% colors = colororder;
% firstLineColor = colors(1,:);
% plot(timeHistory,prnIndex.*acStatsAllTime, ...
%     Color=firstLineColor,LineWidth=1)
% xlim([timeHistory(1) timeHistory(end)])
% ylim([min(prnIndex)-1 max(prnIndex)+1])
% xlabel("Time")
% ylabel("Satellite PRN Index")
% title("Satellite Visibility Chart")
% yticks(prnIndex)
% grid on

%% Ephemeris

if subset==1
Set_Satelliite=Index_tot(any(acStatsAllTime(1:Window,:)));
elseif subset==2
Set_Satelliite=Index_tot(any(acStatsAllTime(1:Window,:)));
QQQ=sum(~isnan(acStatsAllTime(1:Window,Set_Satelliite))*60);
Set_Satelliite=Set_Satelliite(QQQ>500);
else
   Set_Satelliite=Index_tot(b);
end
Ephem=zeros(9,length(Set_Satelliite));
for ii=1:length(Set_Satelliite)

    for jj=1:length(Inclination)

        if Set_Satelliite(ii)<=constellation(jj).PRN_range(end) && Set_Satelliite(ii)>=constellation(jj).PRN_range(1)
            Index=find(constellation(jj).PRN_range==Set_Satelliite(ii));
            [Pos,Vel]=states(constellation(jj).obj(Index),'CoordinateFrame','inertial');
            [a,ecc,incl,RAAN,argp,nu,truelon,arglat,lonper] = ijk2keplerian(Pos(:,1),Vel(:,1));
            Ephem(:,ii)=[a,ecc,deg2rad(incl),deg2rad(RAAN),deg2rad(argp),deg2rad(nu),deg2rad(truelon),...
                deg2rad(arglat),deg2rad(lonper)].';

        end

    end

    

end


fileID = fopen(File_name_GAL,'w');
II=1;
for ii=1:36

if all(ii~=ID_Gal) || isempty(ID_Gal)
fprintf(fileID,'******** Week %d almanac for PRN-%d********\n',GAL_input(ii,end),ii);
fprintf(fileID,'ID:                             %2d \n',ii);
fprintf(fileID,'Health:                         0 \n');
fprintf(fileID,'Eccentricity:                   %1.14e\n',Ephem(2,II));
fprintf(fileID,'Time of Applicability(s):       %6d\n',GAL_input(ii,4));
fprintf(fileID,'Orbital Inclination(rad):       %1.14e\n',Ephem(3,II));
fprintf(fileID,'Rate of Right Ascen(r/s):       0.00000000000000e+00\n');
fprintf(fileID,'SQRT(A)  (m 1/2):               %1.14e\n',sqrt(Ephem(1,II)));
Appo2=Ephem(4,II);%+pi;

if abs(Appo2)<deg2rad(180)
fprintf(fileID,'Right Ascen at Week(rad):       %1.14e\n',Appo2);
else
fprintf(fileID,'Right Ascen at Week(rad):       %1.14e\n',Appo2-2*pi);    
end
fprintf(fileID,'Argument of Perigee(rad):       0.00000000000000e+00\n');
if abs(Ephem(8,II))<deg2rad(180)
fprintf(fileID,'Mean Anom(rad):                 %1.14e\n',Ephem(8,II));
else
fprintf(fileID,'Mean Anom(rad):                 %1.14e\n',Ephem(8,II)-2*pi);    
end
fprintf(fileID,'Af0(s):                         0.00000000000000e+00\n');
fprintf(fileID,'Af1(s/s):                       0.00000000000000e+00\n');
fprintf(fileID,'week:                           %3d\n',GAL_input(ii,end));
fprintf(fileID,'Af1(s/s):                       0.00000000000000e+00\n');
if II>=size(Ephem,2),break,end

II=II+1
else

fprintf(fileID,'******** Week %d almanac for PRN-%2d********\n',GAL_input(ii,end),ii+Traslazione_GAL_Matteo_PRN);
fprintf(fileID,'ID:                             %1d \n',GAL_input(ii,1)+Traslazione_GAL_Matteo_PRN);
fprintf(fileID,'Health:                         0 \n');
fprintf(fileID,'Eccentricity:                   %1.14e\n',GAL_input(ii,3));
fprintf(fileID,'Time of Applicability(s):       %6d\n',GAL_input(ii,4));
fprintf(fileID,'Orbital Inclination(rad):       %1.14e\n',GAL_input(ii,5));
fprintf(fileID,'Rate of Right Ascen(r/s):       0.00000000000000e+00\n');
fprintf(fileID,'SQRT(A)  (m 1/2):               %1.14e\n',GAL_input(ii,7));
fprintf(fileID,'Right Ascen at Week(rad):       %1.14e\n',GAL_input(ii,8)+pi);
fprintf(fileID,'Argument of Perigee(rad):       %1.14e\n',GAL_input(ii,9));
fprintf(fileID,'Mean Anom(rad):                 %1.14e\n',GAL_input(ii,10));
fprintf(fileID,'Af0(s):                         %1.14e\n',GAL_input(ii,11));
fprintf(fileID,'Af1(s/s):                       %1.14e\n',GAL_input(ii,12));
fprintf(fileID,'week:                           %3d\n',GAL_input(ii,end));
end

end


fileID = fopen(File_name_GPS,'w');
for ii=1:32

if all(ii~=ID_GPS) || isempty(ID_GPS)
fprintf(fileID,'******** Week %d almanac for PRN-%d ********\n',GPS_input(ii,end),ii);
fprintf(fileID,'ID:                             %1d \n',ii);
fprintf(fileID,'Health:                         0 \n');
fprintf(fileID,'Eccentricity:                   %1.14e\n',Ephem(2,II));
fprintf(fileID,'Time of Applicability(s):       %6d\n',GPS_input(ii,4));
fprintf(fileID,'Orbital Inclination(rad):       %1.14e\n',Ephem(3,II));
fprintf(fileID,'Rate of Right Ascen(r/s):       0.00000000000000e+00\n');
fprintf(fileID,'SQRT(A)  (m 1/2):               %1.14e\n',sqrt(Ephem(1,II)));
Appo2=Ephem(4,II);%+pi;
if abs(Appo2)<=deg2rad(180)
fprintf(fileID,'Right Ascen at Week(rad):       %1.14e\n',Appo2);
else
fprintf(fileID,'Right Ascen at Week(rad):       %1.14e\n',Appo2-2*pi);    
end
fprintf(fileID,'Argument of Perigee(rad):       0.00000000000000e+00\n');
if abs(Ephem(8,II))<=deg2rad(180)
fprintf(fileID,'Mean Anom(rad):                 %1.14e\n',Ephem(8,II));
else
fprintf(fileID,'Mean Anom(rad):                 %1.14e\n',Ephem(8,II)-2*pi);    
end
fprintf(fileID,'Af0(s):                         0.00000000000000e+00\n');
fprintf(fileID,'Af1(s/s):                       0.00000000000000e+00\n');
fprintf(fileID,'week:                           %3d\n',GPS_input(ii,end));

if II>=size(Ephem,2),break,end

II=II+1

else

fprintf(fileID,'******** Week %d almanac for PRN-%2d********\n',GPS_input(ii,end),ii);
fprintf(fileID,'ID:                             %1d \n',GPS_input(ii,1));
fprintf(fileID,'Health:                         0 \n');
fprintf(fileID,'Eccentricity:                   %1.14e\n',GPS_input(ii,3));
fprintf(fileID,'Time of Applicability(s):       %6d\n',GPS_input(ii,4));
fprintf(fileID,'Orbital Inclination(rad):       %1.14e\n',GPS_input(ii,5));
fprintf(fileID,'Rate of Right Ascen(r/s):       0.00000000000000e+00\n');
fprintf(fileID,'SQRT(A)  (m 1/2):               %1.14e\n',GPS_input(ii,7));
fprintf(fileID,'Right Ascen at Week(rad):       %1.14e\n',GPS_input(ii,8));
fprintf(fileID,'Argument of Perigee(rad):       %1.14e\n',GPS_input(ii,9));
fprintf(fileID,'Mean Anom(rad):                 %1.14e\n',GPS_input(ii,10));
fprintf(fileID,'Af0(s):                         %1.14e\n',GPS_input(ii,11));
fprintf(fileID,'Af1(s/s):                       %1.14e\n',GPS_input(ii,12));
fprintf(fileID,'week:                           %3d\n',GPS_input(ii,end));
end

end

fclose all

% fprintf(fileID,'%6.2f %12.8f\n',A);
% fclose(fileID);
% 
% ******** Week 219 almanac for PRN-01 ********
% ID:                             01
% Health:                         0
% Eccentricity:                   0.00000000000000e+00
% Time of Applicability(s):       589824
% Orbital Inclination(rad):       0.00000000000000e+00
% Rate of Right Ascen(r/s):       0.00000000000000e+00
% SQRT(A)  (m 1/2):               2.81691320000000e+03
% Right Ascen at Week(rad):       0.00000000000000e+00
% Argument of Perigee(rad):       0.00000000000000e+00
% Mean Anom(rad):                 0.00000000000000e+00
% Af0(s):                         0.00000000000000e+00
% Af1(s/s):                       0.00000000000000e+00
% week:                           219


mission.viewer.PlaybackSpeedMultiplier = 200;
play(mission.scenario);