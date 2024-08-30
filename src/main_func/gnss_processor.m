function gnss_processor(rtk,opt,obsr,obsb,nav)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start gnss processor to generate navigation solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Copyright(c) 2020-2025, by Kai Chen, All rights reserved.
%8/12/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global glc gls OptionSatnav OptionRinexDummy
ti=0;

hbar=waitbar(0,'Processing...','Name','GINav', 'CreateCancelBtn', 'delete(gcbf);');
H=get(0,'ScreenSize'); w=600; h=450; x=H(3)/2-w/2; y=H(4)/2-h/2; 
hfig=figure;set(gcf,'Position',[x y w h]);

% edit Heiko HE
PL_fig = figure; set(PL_fig,'Position',[x+w y w h]);
subplot(3,1,1); ylabel('PL [m]'); title('X');
subplot(3,1,2); ylabel('PL [m]'); title('Y');
subplot(3,1,3); ylabel('PL [m]'); title('Z'); xlabel('Epoch')
p1 = plot(0,0,'b*');
hold on
p2 = plot(0,0,'ro');
legend([p1,p2],{'PL','Error'},'AutoUpdate','off');

bias_fig = figure; set(bias_fig,'Position',[x-w y w h]);

figure(hfig);

% initialize rtk sturct
rtk=initrtk(rtk,opt);

% set time span
tspan=timespan(rtk,obsr);
if tspan<=0,error('Time span is zero!!!');end


%% Create DUmmy 

% if OptionSatnav
%     for ii=1:length(obsr)
%     if obsr(ii).sat>=60
%     obsr(ii).P(obsr(ii).P~=0)=obsr(ii).P(obsr(ii).P~=0)+2*10^5;
%     obsr(ii).L(obsr(ii).L~=0)=obsr(ii).L(obsr(ii).L~=0)+2*10^5;
%     else
%     obsr(ii).P(obsr(ii).P~=0)=obsr(ii).P(obsr(ii).P~=0)+2*10^5+15;
%     obsr(ii).L(obsr(ii).L~=0)=obsr(ii).L(obsr(ii).L~=0)+2*10^5+15;
% 
%     end
% 
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ATTENTION THIS BRANCH APPLY ONLY IF YOU HAVE ALL LEO ON GPS  
%% IF YOU WANT TO MIX GPS AND GALILEO YOU HAVE TO MODIFY THE SCRIPT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Noise 
rng(1111);
gls.Activate_Noise_PRange=zeros(glc.MAXSAT*glc.NFREQ,glc.Prealloc);
if any(glc.Activate_Noise_PRange)
    for ii=1:length(glc.Activate_Noise_PRange)
        gls.Activate_Noise_PRange((ii-1)*glc.MAXSAT+1:ii*glc.MAXSAT,:)=glc.Activate_Noise_PRange(ii)*randn(glc.MAXSAT,glc.Prealloc);
    end
end
gls.Activate_Noise_LRange=zeros(glc.MAXSAT*glc.NFREQ,glc.Prealloc);
if any(glc.Activate_Noise_LRange)
    for ii=1:length(glc.Activate_Noise_LRange)
        gls.Activate_Noise_LRange((ii-1)*glc.MAXSAT+1:ii*glc.MAXSAT,:)=glc.Activate_Noise_LRange(ii)*randn(glc.MAXSAT,glc.Prealloc);
    end
end

%% Clock
if any(glc.Activate_ClockError)
    for ii=1:length(glc.Activate_ClockError)
        gls.Activate_ClockError(ii,:)=ones(1,glc.Prealloc);
    end
end


if any(glc.Activate_ClockError)
    for ii=1:length(glc.Activate_ClockError)
        gls.Activate_ClockError(ii,:)=glc.Activate_ClockError(ii)*ones(1,glc.Prealloc);
    end
end

%% Ephemeris(PPP) - SISRE
ERR_SISRE=zeros(glc.MAXSAT,glc.Prealloc);

%% GPS
rng(1112);
if any(glc.ADDEPHERR_PPP_LEO)

    Coeff_1=rand(glc.MAXSAT,glc.Prealloc);
    Coeff_2=rand(glc.MAXSAT,glc.Prealloc);
    Coeff_3=rand(glc.MAXSAT,glc.Prealloc);
    Phase=rand(glc.MAXSAT,glc.Prealloc);

    for ii=1:length(glc.ID_LEO)
            ERR_SISRE(glc.ID_LEO(ii),:)=periodic_parabola(glc.ADDEPHERR_PPP_LEO(1)*Coeff_1(ii,:),...
            glc.ADDEPHERR_PPP_LEO(2)*Coeff_2(ii,:)/200, glc.ADDEPHERR_PPP_LEO(3)*Coeff_3(ii,:)/500, Phase,300, glc.Prealloc);
    end
end

rng(1111);

%% GPS
if any(glc.ADDEPHERR_PPP_GPS)

    Coeff_1=rand(glc.NSATGPS,glc.Prealloc);
    Coeff_2=rand(glc.NSATGPS,glc.Prealloc);
    Coeff_3=rand(glc.NSATGPS,glc.Prealloc);
    Phase=rand(glc.NSATGPS,glc.Prealloc);

    for ii=1:glc.NSATGPS
        if any(glc.ID_LEO==ii)
            continue
        else
        ERR_SISRE(ii,:)=periodic_parabola(glc.ADDEPHERR_PPP_GPS(1)*Coeff_1(ii,:),...
            glc.ADDEPHERR_PPP_GPS(2)*Coeff_2(ii,:)/200, glc.ADDEPHERR_PPP_GPS(3)*Coeff_3(ii,:)/500, Phase,300, glc.Prealloc);
        end
    end
end

%% GAL
if any(glc.ADDEPHERR_PPP_GAL)

    Coeff_1=rand(glc.NSATGPS,glc.Prealloc);
    Coeff_2=rand(glc.NSATGPS,glc.Prealloc);
    Coeff_3=rand(glc.NSATGPS,glc.Prealloc);
    Phase=rand(glc.NSATGPS,glc.Prealloc);

    for ii=glc.NSATGPS+glc.NSATGLO+1:glc.NSATGPS+glc.NSATGLO+glc.NSATGAL
        if any(glc.ID_LEO==ii)
            continue
        else
        ERR_SISRE(ii,:)=periodic_parabola(glc.ADDEPHERR_PPP_GAL(1)*Coeff_1(ii,:),...
            glc.ADDEPHERR_PPP_GAL(2)*Coeff_2(ii,:)/200, glc.ADDEPHERR_PPP_GAL(3)*Coeff_3(ii,:)/500, Phase,300, glc.Prealloc);
        end
    end
end

while 1

    %% Counter 
    gls.ti=ti;

    if ti>tspan,break;end
    
    % search rover obs
    [obsr_,nobs,obsr]=searchobsr(obsr);
    if nobs<0
        str=sprintf('Processing... %.1f%%',100);
        waitbar(ti/tspan,hbar,str);
        break;
    end
    % exclude rover obs
    [obsr_,nobs]=exclude_sat(obsr_,rtk);
    if nobs==0,continue;end

    if opt.mode>=glc.PMODE_DGNSS&&opt.mode<=glc.PMODE_STATIC
        % search base obs
        [obsb_,nobs]=searchobsb(obsb,obsr_(1).time);
        % exclude base obs
        if nobs~=0,[obsb_,~]=exclude_sat(obsb_,rtk);end
    else
        obsb_=NaN;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PARSE Errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MENx
if OptionRinexDummy==1 && glc.Activate_ErrorPRange==1
    for ii=1:length(obsr_)


   if any(obsr_(ii).sat==glc.ID_LEO)
      obsr_(ii).P(obsr_(ii).P~=0)=obsr_(ii).P(obsr_(ii).P~=0)+glc.LEO_FactP*gls.Activate_Noise_PRange(Index(obsr_(ii).P~=0),ti+1).';
    obsr_(ii).L(obsr_(ii).L~=0)=obsr_(ii).L(obsr_(ii).L~=0)+glc.LEO_FactL*gls.Activate_Noise_LRange(Index(obsr_(ii).P~=0),ti+1).';
      continue
   end


    Index=[obsr_(ii).sat,obsr_(ii).sat+glc.MAXSAT,obsr_(ii).sat+2*glc.MAXSAT];
    if obsr_(ii).sat>=60
    obsr_(ii).P(obsr_(ii).P~=0)=obsr_(ii).P(obsr_(ii).P~=0)+glc.GAL_FactP*gls.Activate_Noise_PRange(Index(obsr_(ii).P~=0),ti+1).';
    obsr_(ii).L(obsr_(ii).L~=0)=obsr_(ii).L(obsr_(ii).L~=0)+glc.GAL_FactL*gls.Activate_Noise_LRange(Index(obsr_(ii).P~=0),ti+1).';
    else
    obsr_(ii).P(obsr_(ii).P~=0)=obsr_(ii).P(obsr_(ii).P~=0)+gls.Activate_Noise_PRange(Index(obsr_(ii).P~=0),ti+1).';
    obsr_(ii).L(obsr_(ii).L~=0)=obsr_(ii).L(obsr_(ii).L~=0)+gls.Activate_Noise_LRange(Index(obsr_(ii).P~=0),ti+1).';
    end
    end
end

if OptionSatnav==1 && glc.Activate_EPHErr==1

    %% SISRE

    for ii=1:length(obsr_)
    Index=[obsr_(ii).sat,obsr_(ii).sat+glc.MAXSAT,obsr_(ii).sat+2*glc.MAXSAT];
%     if obsr_(ii).sat>=60
        %% LEO GAL
        if any(obsr_(ii).sat==glc.ID_LEO)
        APPO=ERR_SISRE(obsr_(ii).sat,ti+1)*obsr_(ii).P~=0;
        obsr_(ii).P(obsr_(ii).P~=0)=obsr_(ii).P(obsr_(ii).P~=0)+APPO;
        end
%     else
%         if any(obsr_(ii).sat==glc.ID_LEO)
%         APPO=ERR_SISRE(obsr_(ii).sat,ti+1)*obsr_(ii).P~=0;
%         obsr_(ii).P(obsr_(ii).P~=0)=obsr_(ii).P(obsr_(ii).P~=0)+APPO;
%         end
    end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Single Frequency LEO Filte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if opt.LEOsingleFreq==1

        for ii=1:length(obsr_)
            if any(obsr_(ii).sat==glc.ID_LEO)
                obsr_(ii).P(2:3)=0;
                obsr_(ii).L(2:3)=0;
                obsr_(ii).D(2:3)=0;

            end

        end

end


    % solve an epoch
    [rtk,~]=gnss_solver(rtk,obsr_,obsb_,nav);
    
    if rtk.sol.stat~=0
        % write solution to output file
        outsol(rtk);
        % kinematic plot
        plot_trajectory_kine(hfig,rtk);

        % Heiko (HE): Protection level plot
        plot_protection_level(PL_fig,rtk);

        % Heiko (HE): Bias plot
        cm = colormap;
        cm = cm(round(linspace(1,256,150)),:);

        figure(bias_fig);
        hold on;
        for i = 1:150
        plot(ti,gls.Log.PPP.STATE(i+15),'marker','*','Color',cm(i,:));
        end
        
        figure(hfig); % reset current figure
    else
        [week,sow]=time2gpst(obsr_(1).time);
        fprintf('Warning:GPS week = %d sow = %.3f,GNSS unavailable!\n',week,sow);
    end
    
    % update progress bar
    ti=ti+1;
    str=sprintf('Processing... %.1f%%',100*ti/tspan);
    waitbar(ti/tspan,hbar,str);
     
end

close(hbar);

return

