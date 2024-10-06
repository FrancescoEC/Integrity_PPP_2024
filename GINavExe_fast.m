%--------------------------------------------------------------------------
%
%                    %%%%%  %%%%%  %   %    %    %   %
%                    %        %    %%  %
%                    %  %%%   %    % % %   %%%    % %
%                    %   %    %    %  %%
%                    %%%%%  %%%%%  %   %  %   %    %
%
%
%--------------------------------------------------------------------------
%                            GINav v0.1.0
%
% Copyright (C) 2020-2025, by Kai Chen, All rights reserved.
%--------------------------------------------------------------------------
% This program is free software,it's distributed in the hope that provide a
% useful tool for carrying out  GNSS/INS or GNSS-related research,but
% WITHOUT ANY WARRANTY
%--------------------------------------------------------------------------

clc
clear all %#ok
close all
clear global
% global motionV2 
% global satdataV1A2 OptionSatnav
addpath(genpath(pwd))
global_variable;

MODES = 0;%'HW';
% Test_name='Skydel';
Test_name='LEOPNT_clock10Mhz';


%% Global Overwrite
    glc.BRDperPrec=1;
    glc.EphErrOverW=0.1000;
    glc.OUTNAME='Test_PPP_LOW_MEO_MN_GAL_ONLY_EPH_NON_NON_REPORT_9sat';

%LEOPNT_clock10Mhz
ConfFile='GINav_PPP_LEO_MN_GAL_ONLY_DUMMY_REPORT.ini';
% execute GINavCfg to configure input file
root=pwd;
Parth_conf=fullfile(root,'conf','PPP');
Parth_data=fullfile(root,'data',Test_name);

%% File Obs
%file_obs='MOSA00____R_20231790000_01H_01S_MO.rnx';
%file_obs='SEPT00____R_20231790000_59M_01S_MO.rnx';
file_obs='rinex-obs_V1_A1-static_vehicle_MO.rnx';

file_nav='SEPT00____R_20231790000_01H_MN.rnx';
file_sp3 ='GRG0OPSFIN_20231790000_01D_05M_ORB.SP3';
file_clk= 'GRG0OPSFIN_20231790000_01D_30S_CLK.CLK';
File_Error_BRD='BRDC-MGEX_Final_SP3CLKREF_2023_356.mat';
File_Error_HAS='HASSIS-MGEX_Final_SP3CLKREF_2023_356.mat';

opt=decode_cfg(gls.default_opt,fullfile(Parth_conf,ConfFile));
Allo=dir(Parth_data);
readSpirent=1;

for ii=1:length(Allo)
    if length(Allo(ii).name)>4
        switch Allo(ii).name(end-3:end)
            case ''
                file.path='';
            case '.rnx'
                if strcmp(Allo(ii).name(end-6:end),'_MO.rnx')
                    file.obsr=fullfile(Parth_data,file_obs);
                else
                    file.beph=fullfile(Parth_data,file_nav);
                end
            case '.SP3'
                file.sp3=fullfile(Parth_data,file_sp3);
            case '.CLK'
                file.clk=fullfile(Parth_data,file_clk);
            case '.csv'
                if strcmp(Allo(ii).name(1:8),'sat_data') && readSpirent
                    satdataV1A2 = satnav(fullfile(Parth_data,'sat_data_V1A1.csv'));
                    %Time = importTime(fullfile(Parth_data,'sat_data_V1A1.csv'));
                    %epoch2time(2023 06 28 00 00 20.0000000)
                    InitTime=epoch2time([2023 06 28 00 00 00]);
                    SkipTime=epoch2time([2023 06 28 00 00 00]);
                else
                    motionV2 = motion(fullfile(Parth_data,'motion_V1.csv'));      
                end

               case '.23O'
                    file.obsr=fullfile(Parth_data,Allo(ii).name);
  
                case '.23P'
                    file.beph=fullfile(Parth_data,Allo(ii).name);
  
% %                 if strcmp(Allo(ii).name(end-6:end),'_MO.rnx')
% %                     file.obsr=fullfile(Parth_data,Allo(ii).name);
% %                 else
% %                     file.beph=fullfile(Parth_data,'MOSA00____R_20231790000_01H_MN.rnx');
% %                 end
% 
% %                 end

        end

    end
end

file.obsb='';
file.atx='';
file.dcb={'','',''};
file.dcb_mgex='';
file.erp='';
file.blq ='';
file.imu ='';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEO_Section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SearchLEO
%glc.ID_BAN=[4,7,9,16,20,27,30,61,93];
%glc.ID_BAN=[4,7,9,16,20,27,30,72,85,93];
glc.ID_BAN=[4,7,9,16,20,27,30];

if ~isempty(SetGAL_LEO)
    SetGAL_LEO=SetGAL_LEO+59;
end
glc.ID_LEO=[SetGPS_LEO(:);SetGAL_LEO(:)];
glc.ID_GAL=SetGAL(:)+59;
glc.ID_GPS=SetGPS(:);
glc.ConfFile=ConfFile;



% if gui_flag==0
%     return;
% end
% load(fullfile(Parth_data,'satdataV1A1.mat'));

% vect=zeros(size(satdataV1A1,1),1);
% start=epoch2time([2023 06 28 0 0 0]);
% motionV1=table2struct(satdataV1A1,'ToScalar',true);
% motionV2=motionV1;
% motionV2.TimeGPS=vect;
% % for ii=1:length(vect)
% %     motionV2.TimeGPS(ii,1)=start.time+motionV1.Time_ms(ii,1)/1000;
% % end

opt.Galileo_REF=0;
opt.LEO_Aug=0;
opt.LEO_pseudoonly=0;
opt.LEOsingleFreq=0;
opt.Slant_TEC=1;
opt.LEO_IONO_MODEL=0;
opt.CumdtrFlag=1;
opt.MAXLEO=1000;

if opt.LEO_Aug
gls.EPH_LEO_ERR=1;
else
gls.EPH_LEO_ERR=0;
end


if MODES==0
opt.ClkSim=2;
OptionSatnav=1;
OptionRinexDummy=1;
OptionClock=0;
OptionSatNavEph=1;
else
opt.ClkSim=2;
OptionSatnav=1;
OptionRinexDummy=0;
OptionClock=0;
OptionSatNavEph=0;
end





% motionV2(:,1)=vect;
opt.L5=1;


%% %% GLOBAL OVERWRITE %%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New glc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

glc.Alloc=20000;
glc.Prealloc=86400;

%% Increment of variances 
gls.STD_BRDCCLK=0;
glc.FlagBRD2FINE_PPP=1;
gls.EPH_Fact=10;
glc.URA_FINE=0.2^2;
glc.URA_BRD=5;

glc.FlagBRD2FINE_SPP=0;
glc.STD_LEO=10;



if MODES==0
%% Noise Characterization
glc.Activate_Noise_PRange=1*[0.5;0.5;0]; % Per Frequency
glc.Activate_Noise_LRange=1*[0.01;0.01;0];% Per Frequency
glc.GAL_FactP=0.9; % Scale factor GPS vs Galileo 
glc.GAL_FactL=0.9; % Scale factor GPS vs Galileo
glc.LEO_FactP=1; % Scale factor GPS vs Galileo 
glc.LEO_FactL=1; % Scale factor GPS vs Galileo
glc.Activate_ErrorPRange=1; 
glc.Activate_ErrorLrRange=1;
glc.Activate_EPHErr=0;
else
glc.Activate_ErrorPRange=0; 
glc.Activate_ErrorLrRange=0;
glc.Activate_EPHErr=0;
end

if MODES==0
%% Intersystem clock error 
glc.Activate_ClockError=1*[0.1;0;0.05;0;0]; % Per System (5)
glc.Activate_ClockError_rate=1*[64;0;64;0;0]; % Per System (5)
glc.Activate_ClockError_LEO=0;
else
glc.Activate_ClockError=0*[0.1;0;0.05;0;0]; % Per System (5)
glc.Activate_ClockError_rate=0*[64;0;64;0;0]; % Per System (5)
glc.Activate_ClockError_LEO=0;
end


%% Ephemeris Error 

glc.FileBRD=[];
glc.FileHAS=[];

glc.EPHERR=[];
glc.EPHERR_LEO=[];
glc.Scale_GNSS=0;
glc.Scale_LEO=0;
glc.HASBRDC_GNSS=0;
glc.HASBRDC_LEO=0;
glc.ADDEPHERR_LEO_MODEL=0;
glc.Shaffle_GNSS=0;
glc.Shaffle_LEO=0;
glc.ref_sat_err='G04';
glc.SPPBRD=1;

glc.ADDEPHERR_FromFile_SPP=0; % Valid Only for GPS and Galileo  (To be implemented) 
glc.ADDEPHERR_SPP_GPS=[0;0;0];
glc.ADDEPHERR_SPP_GAL=[0;0;0];
glc.ADDEPHERR_SPP_LEO=[0;0;0];

glc.ADDEPHERR_FromFile_PPP=0; % Valid Only for GPS and Galileo   
glc.ADDEPHERR_PPP_GPS=[0;0;0];
glc.ADDEPHERR_PPP_GAL=[0;0;0]; % SISRE Model Coefficient
glc.ADDEPHERR_PPP_LEO=glc.ADDEPHERR_LEO_MODEL*[0.1,0.01,0.001]/2; % SISRE Model Coefficient
glc.ADDEPHERR_PPP_LEO_Component_Factor=1*[1,1,2]; % SISRE Model Coefficient

% In case of LEO-HAS
%glc.ADDEPHERR_PPP_LEO=2*[0.05,0.001,0.0001]; % SISRE Model Coefficient


%% 
glc.Relativity=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New Gls
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% New gls
gls.Log.SPP=[];
gls.Log.PPP=[];


% execute positioning
exepos(opt,file);
%global gls
file_obs_1=file_obs;
file_obs_1(end-3:end)=[];

Name=[file_obs_1,'_PPP_K',glc.OUTNAME,'.pos'];
% Name=[file_obs_1,'_SPP',glc.OUTNAME,'.pos'];

File_1=readGINavsol(fullfile('result',Name));
plot_pva_err_station(File_1,[4402957.8071   666835.8982  4550862.1942],4,3000)

%plot_pva_err_station(File_1,[6378137.0000        0.0000        0.0000],4)


% figure, plot(gls.Log.SPP.sim,gls.Log.SPP.Nsat)

figure, plot(gls.Log.PPP.sim,gls.Log.PPP.Nsat,'LineWidth',3)
if opt.LEO_Aug==0
title('N-SAT Galileo Only')
else
title('N-SAT Galileo+LEO ')
end
xlabel('Simulation Time [s]')
xlim([0,3000])

figure, plot((gls.Log.PPP.RESIDUALS_SATNAV_3(151:300,:).'))
grid on
xlabel('Simulation Time [s]')
ylabel('[m]')
title('Tropo Error')
xlim([0,3000])


figure, plot((gls.Log.PPP.RESIDUALS_SATNAV_4(151:300,:).'))
grid on
xlabel('Simulation Time [s]')
title('Iono Residuals')
ylabel('[m]')
xlim([0,3000])

figure, plot((gls.Log.PPP.RESIDUALS_SATNAV_2(151:300,:).'))
grid on
xlabel('Simulation Time [s]')
title('Clock Residuals')
ylabel('[m]')
xlim([0,3000])

%% Ephemeris ERROR

figure,
subplot 311
plot(glc.POS_ERR_ECEF(glc.ID_GAL,1:3000,1).');
grid on
xlabel('Simulation Time [s]')
title('Eph Error (GALILEO)')
ylabel('XECEF [m]')
subplot 312
plot(glc.POS_ERR_ECEF(glc.ID_GAL,1:3000,2).');
grid on
xlabel('Simulation Time [s]')
ylabel('YECEF [m]')
subplot 313
plot(glc.POS_ERR_ECEF(glc.ID_GAL,1:3000,3).');
grid on
xlabel('Simulation Time [s]')
ylabel('ZECEF [m]')

figure,
subplot 311
plot(glc.POS_ERR_ECEF(glc.ID_LEO,1:3000,1).');
grid on
xlabel('Simulation Time [s]')
title('Eph Error (LEO)')
ylabel('XECEF [m]')
subplot 312
plot(glc.POS_ERR_ECEF(glc.ID_LEO,1:3000,2).');
grid on
xlabel('Simulation Time [s]')
ylabel('YECEF [m]')
subplot 313
plot(glc.POS_ERR_ECEF(glc.ID_LEO,1:3000,3).');
grid on
xlabel('Simulation Time [s]')
ylabel('ZECEF [m]')

mkdir(glc.OUTNAME);
save(fullfile(glc.OUTNAME,[glc.OUTNAME,num2str(now*10000),'.mat']));
saveallfigures(fullfile(pwd,glc.OUTNAME));


% GPS SVs List = [4   7   9  16  19  20  27  30 ]
% GAL SVs List = [1   2   3   8  11  12  13  14  18  19  20  22  25  26  29  31  32  33  34  35  36 ]
% LEOPNT-GPS SVs List = [1   2   3   5   6   8  10  11  12  13  14  15  17  18  21  22  23  24  25  26  28  29  31  32 ]
% LEOPNT-GAL SVs List = [4   5   6   7   9  10  15  16  17  21  23  24  27  28  30 ]
%--------------------------------------------------------------------------

