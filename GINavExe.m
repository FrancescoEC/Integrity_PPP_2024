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
addpath(genpath(pwd))
global_variable;

% execute GINavCfg to configure input file
[opt,file,gui_flag]=GINavCfg;

if gui_flag==0
    return;
end 
load('C:\Users\menzifr\Documents\3_My_SW\LEO_SPS\GINav-main(1)\GINav-main\data\LEO\satdataV1A1.mat');

vect=zeros(size(satdataV1A1,1),1);
start=epoch2time([2023 06 28 0 0 0]);
motionV1=table2struct(satdataV1A1,'ToScalar',true);
motionV2=motionV1;
motionV2.TimeGPS=vect;
opt.LEO_Aug=1;
OptionSatnav=0;
opt.LEO_pseudoonly=0;

for ii=1:length(vect)
    motionV2.TimeGPS(ii,1)=start.time+motionV1.Time_ms(ii,1)/1000;
end

% motionV2(:,1)=vect;
opt.L5=1;
% execute positioning
exepos(opt,file);

%--------------------------------------------------------------------------

