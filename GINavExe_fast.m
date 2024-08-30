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
% Test_name='Skydel';
Test_name='200LEO_T';


%% Global Overwrite
    glc.BRDperPrec=1;
    glc.EphErrOverW=0.1000;
    glc.OUTNAME='Test_PPP_SAtnav_LEO200_LEO_dummy';


ConfFile='GINav_PPP_LEO_Satnav.ini';
% execute GINavCfg to configure input file
root=pwd;
Parth_conf=fullfile(root,'conf','PPP');
Parth_data=fullfile(root,'data',Test_name);

%% File Obs
%file_obs='MOSA00____R_20231790000_01H_01S_MO.rnx';
file_obs='rinex-obs_V1_A1_MO.rnx';
% file_nav='MOSA179Z.23P';
file_nav='MOSA1791.23P'; % edit Heiko
file_sp3 ='GRG0OPSFIN_20231790000_01D_05M_ORB.SP3';
file_clk= 'GRG0OPSFIN_20231790000_01D_30S_CLK.CLK';

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
if ~isempty(SetGAL_LEO)
    SetGAL_LEO=SetGAL_LEO+60;
end
glc.ID_LEO=[SetGPS_LEO;SetGAL_LEO];

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


opt.LEO_Aug=1;
OptionSatnav=1;
OptionRinexDummy=1;
opt.LEO_pseudoonly=1;
opt.LEOsingleFreq=1;
opt.LEO_thr_pos=3000000;

if opt.LEO_Aug
gls.EPH_LEO_ERR=1;
else
gls.EPH_LEO_ERR=0;
end


% motionV2(:,1)=vect;
opt.L5=1;
% execute positioning
exepos(opt,file);
%global gls
File_1=readGINavsol(fullfile('result','rinex-obs_V1_MO_SPPTest_SPP_SAtnav.pos'));
plot_pva_err_station(File_1,[4402957.80710000	666835.898200000	4550862.19420000],4)
figure, plot(gls.Log.PPP.sim,gls.Log.PPP.Nsat)
save('Sim_rinex-obs_V1_MO_SPPTest_SPP_SAtnav')

% Edit Heiko (HE)
figure
plot(gls.Log.PPP.STATE(16:end,1:gls.ti))

%--------------------------------------------------------------------------

