function exepos(opt,file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute positioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright(c) 2020-2025, by Kai Chen, All rights reserved.
%8/12/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global glc gls motionV2 satdataV1A2 ClockBias InitTime
global File_Error_BRD File_Error_HAS
tic

rtk=gls.rtk;

% read input file
[filepath,name,ext] = fileparts(file.obsr);
mkdir(fullfile(filepath,glc.ConfFile));
matFiles = dir(fullfile(fullfile(filepath,glc.ConfFile), '*.mat'));
if isempty(matFiles)
    [obsr,obsb,nav,imu]=read_infile(opt,file);
    save(fullfile(fullfile(filepath,glc.ConfFile),'obsr.mat'),'obsr','-v7.3');
    save(fullfile(fullfile(filepath,glc.ConfFile),'obsb.mat'),'obsb','-v7.3');
    save(fullfile(fullfile(filepath,glc.ConfFile),'nav.mat'),'nav','-v7.3');
    save(fullfile(fullfile(filepath,glc.ConfFile),'imu.mat'),'imu','-v7.3');

else

    load(fullfile(filepath,glc.ConfFile,'obsr.mat'));
    load(fullfile(filepath,glc.ConfFile,'obsb.mat'));
    load(fullfile(filepath,glc.ConfFile,'nav.mat'));
    load(fullfile(filepath,glc.ConfFile,'imu.mat'));

end

% read external clock File
if  isfile(fullfile(filepath,'mosaic_pvt_ecef.mat'))
    load(fullfile(filepath,'mosaic_pvt_ecef.mat'));
else
    opt.ClkSim=2;
end
if opt.ClkSim==1
    ClockBias.Bias=0*rx_pvat_ecef.RxClkBias*10^-3*glc.CLIGHT;
    ClockBias.TimeTag=rx_pvat_ecef.TimeTag+10*365*86400+604800;
elseif opt.ClkSim==0
    ClockBias.Bias=rx_pvat_ecef.RxClkBias*10^-3*glc.CLIGHT;
    ClockBias.TimeTag=rx_pvat_ecef.TimeTag+10*365*86400+604800;
else
    ClockBias.Bias=zeros(glc.Prealloc,1);
    ClockBias.TimeTag=InitTime.time:InitTime.time+glc.Prealloc;
end

% Read Ephemeris Error
cons_ID='GE';
Window_Err=0:86400-1;
Flag_BRD=0;
if isfile(fullfile(filepath,'EPH_errors',File_Error_BRD))
    load(fullfile(filepath,'EPH_errors',File_Error_BRD));
    POS_ERR_ECEF_BRD=POS_ERR_ECEF;
    CLK_ERR_SEC_BRD=CLK_ERR_SEC;
    Flag_BRD=1;
end
Flag_HAS=0;
if isfile(fullfile(filepath,'EPH_errors',File_Error_HAS))
    load(fullfile(filepath,'EPH_errors',File_Error_HAS));
    POS_ERR_ECEF_HAS=POS_ERR_ECEF;
    CLK_ERR_SEC_HAS=CLK_ERR_SEC;
    Flag_HAS=1;
end

% Scale Ephemeris Error from Ephemeris Error

glc.POS_ERR_ECEF=zeros(glc.MAXSAT,length(Window_Err),3);
glc.CLK_ERR=zeros(glc.MAXSAT,length(Window_Err),1);
glc.POS_ERR_ECEF_SPP=zeros(glc.MAXSAT,length(Window_Err),3);
glc.CLK_ERR_SPP=zeros(glc.MAXSAT,length(Window_Err),1);

cons_sat=[32,36];
rng(123456);
Err_rand=randn(glc.MAXSAT,1);

if Flag_BRD==1 || Flag_HAS==1
    for cons = 1:2

        for sat = 1:cons_sat(cons)

            satid    = [cons_ID(cons) num2str(sat,'%02.f')];

            switch cons_ID(cons)

                case 'E'
                    % find set of ephemeris for the current galileo satellite
                    ind_sat  =  sat +59;
                case 'G'
                    % find set of ephemeris for the current GPS satellite
                    ind_sat  =  sat ;
            end
            Factor_2=1;
            if ~isfield(POS_ERR_ECEF,satid) || ~isfield(CLK_ERR_SEC_BRD,satid) || isempty(ind_sat)
                disp(['Relased with',satid])
                satid=glc.ref_sat_err;
                Factor_2=Err_rand(ind_sat);

            end
%             Flag_BRD_2=0;
            Flag_HAS_2=0;
            %% Case LEO
            if any(ind_sat==glc.ID_LEO)
                if glc.HASBRDC_LEO==1
                    Flag_HAS_2=1;
                end
                Factor=glc.Scale_LEO*Factor_2;
            else
                if glc.HASBRDC_GNSS==1
%                     Flag_BRD_2=1;
%                 else
                    Flag_HAS_2=1;
                end
                Factor=glc.Scale_GNSS*Factor_2;
            end

            if Flag_BRD

                if glc.SPPBRD==0

                    if all(isnan(POS_ERR_ECEF_BRD.(satid).X)) || all(isnan(CLK_ERR_SEC_BRD.(satid).clk))
                        disp(['Relased with',satid])
                        satid=glc.ref_sat_err;
                        Factor=Err_rand(ind_sat)*Factor;
                    end
                    InterpPOS_X=interp1(POS_ERR_ECEF_BRD.(satid).t-POS_ERR_ECEF_BRD.(satid).t(1),...
                        POS_ERR_ECEF_BRD.(satid).X,Window_Err);
                    InterpPOS_X(isnan(InterpPOS_X))=nanmean(InterpPOS_X);
                    glc.POS_ERR_ECEF(ind_sat,:,1)=Factor*InterpPOS_X;
                    InterpPOS_Y=interp1(POS_ERR_ECEF_BRD.(satid).t-POS_ERR_ECEF_BRD.(satid).t(1),...
                        POS_ERR_ECEF_BRD.(satid).Y,Window_Err);
                    InterpPOS_Y(isnan(InterpPOS_Y))=nanmean(InterpPOS_Y);
                    glc.POS_ERR_ECEF(ind_sat,:,2)=Factor*InterpPOS_Y;
                    InterpPOS_Z=interp1(POS_ERR_ECEF_BRD.(satid).t-POS_ERR_ECEF_BRD.(satid).t(1),...
                        POS_ERR_ECEF_BRD.(satid).Z,Window_Err);
                    InterpPOS_Z(isnan(InterpPOS_Z))=nanmean(InterpPOS_Z);
                    glc.POS_ERR_ECEF(ind_sat,:,3)=Factor*InterpPOS_Z;
                    InterpCLK=interp1(CLK_ERR_SEC_BRD.(satid).t-CLK_ERR_SEC_BRD.(satid).t(1),...
                        CLK_ERR_SEC_BRD.(satid).clk,Window_Err);
                    glc.CLK_ERR(ind_sat,:,1)=Factor*InterpCLK*glc.CLIGHT;

                else
                    if all(isnan(POS_ERR_ECEF_BRD.(satid).X)) || all(isnan(CLK_ERR_SEC_BRD.(satid).clk))
                        disp(['Relased with',satid])
                        satid=glc.ref_sat_err;
                        Factor=Err_rand(ind_sat)*Factor;
                    end
                    InterpPOS_X=interp1(POS_ERR_ECEF_BRD.(satid).t-POS_ERR_ECEF_BRD.(satid).t(1),...
                        POS_ERR_ECEF_BRD.(satid).X,Window_Err);
                    InterpPOS_X(isnan(InterpPOS_X))=nanmean(InterpPOS_X);
                    glc.POS_ERR_ECEF_SPP(ind_sat,:,1)=Factor*InterpPOS_X;
                    glc.POS_ERR_ECEF(ind_sat,:,1)=Factor*InterpPOS_X;
                    InterpPOS_Y=interp1(POS_ERR_ECEF_BRD.(satid).t-POS_ERR_ECEF_BRD.(satid).t(1),...
                        POS_ERR_ECEF_BRD.(satid).Y,Window_Err);
                    InterpPOS_Y(isnan(InterpPOS_Y))=nanmean(InterpPOS_Y);
                    glc.POS_ERR_ECEF_SPP(ind_sat,:,2)=Factor*InterpPOS_Y;
                    glc.POS_ERR_ECEF(ind_sat,:,2)=Factor*InterpPOS_Y;        
                    InterpPOS_Z=interp1(POS_ERR_ECEF_BRD.(satid).t-POS_ERR_ECEF_BRD.(satid).t(1),...
                        POS_ERR_ECEF_BRD.(satid).Z,Window_Err);
                    InterpPOS_Z(isnan(InterpPOS_Z))=nanmean(InterpPOS_Z);
                    glc.POS_ERR_ECEF_SPP(ind_sat,:,3)=Factor*InterpPOS_Z;
                    glc.POS_ERR_ECEF(ind_sat,:,3)=Factor*InterpPOS_Z;

                    InterpCLK=interp1(CLK_ERR_SEC_BRD.(satid).t-CLK_ERR_SEC_BRD.(satid).t(1),...
                        CLK_ERR_SEC_BRD.(satid).clk,Window_Err);
                    glc.CLK_ERR_SPP(ind_sat,:,1)=Factor*InterpCLK*glc.CLIGHT;
                    glc.CLK_ERR(ind_sat,:,1)=Factor*InterpCLK*glc.CLIGHT;                    

                end
            end
            if Flag_HAS && Flag_HAS_2
                if all(isnan(POS_ERR_ECEF_HAS.(satid).X))
                    disp(['Relased with',satid])
                    satid=glc.ref_sat_err;
                    Factor=Err_rand(ind_sat)*Factor;

                end
                InterpPOS_X=interp1(POS_ERR_ECEF_HAS.(satid).t-POS_ERR_ECEF_HAS.(satid).t(1),...
                    POS_ERR_ECEF_HAS.(satid).X,Window_Err);
                InterpPOS_X(isnan(InterpPOS_X))=nanmean(InterpPOS_X);
                glc.POS_ERR_ECEF(ind_sat,:,1)=Factor*InterpPOS_X;
                InterpPOS_Y=interp1(POS_ERR_ECEF_HAS.(satid).t-POS_ERR_ECEF_HAS.(satid).t(1),...
                    POS_ERR_ECEF_HAS.(satid).Y,Window_Err);
                InterpPOS_Y(isnan(InterpPOS_Y))=nanmean(InterpPOS_Y);
                glc.POS_ERR_ECEF(ind_sat,:,2)=Factor*InterpPOS_Y;
                InterpPOS_Z=interp1(POS_ERR_ECEF_HAS.(satid).t-POS_ERR_ECEF_HAS.(satid).t(1),...
                    POS_ERR_ECEF_HAS.(satid).Z,Window_Err);
                InterpPOS_Z(isnan(InterpPOS_Z))=nanmean(InterpPOS_Z);
                glc.POS_ERR_ECEF(ind_sat,:,3)=Factor*InterpPOS_Z;
                InterpCLK=interp1(CLK_ERR_SEC_HAS.(satid).t-CLK_ERR_SEC_HAS.(satid).t(1),...
                    CLK_ERR_SEC_HAS.(satid).clk,Window_Err);
                glc.CLK_ERR(ind_sat,:,1)=Factor*InterpCLK*glc.CLIGHT;
            end

        end

    end

end




% initlize output file
rtk=initoutfile(rtk,opt,file,obsr);

% high efficiency by converting struct to matrix
obsr=adjobs(obsr,opt);
obsb=adjobs(obsb,opt);
nav =adjnav(nav,opt);

% set base position for relative positioning mode
if opt.mode>=glc.PMODE_DGNSS&&opt.mode<=glc.PMODE_STATIC
    rtk=baserefpos(rtk,opt,obsb,nav);
end

% set anttena parameter for satellite and reciever
if opt.mode~=glc.PMODE_SPP
    nav=L2_L5pcv_copy(nav);
    if isfield(obsr,'sta'),stas(1)=obsr.sta;end
    if isfield(obsb,'sta'),stas(2)=obsb.sta;end
    time0.time=obsr.data(1,1);time0.sec=obsr.data(1,2);
    [nav,opt]=setpcv(time0,opt,nav,nav.ant_para,nav.ant_para,stas);
end
% InitTime=obsr.data(1,1);
% process all data
if opt.ins.mode==glc.GIMODE_OFF
    % gnss
    gnss_processor(rtk,opt,obsr,obsb,nav);
elseif opt.ins.mode==glc.GIMODE_LC||opt.ins.mode==glc.GIMODE_TC
    % gnss/ins integration
    gi_processor(rtk,opt,obsr,obsb,nav,imu);
end

toc

return

