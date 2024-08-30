function [rtk,stat]=pppos(rtk,obs,nav,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% precise point positioning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:  rtk - rtk control struct
%        obs - observations
%        nav - navigation message
%output: rtk - rtk control struct
%        stat - state (0:error 1:ok)

% Version 27-06-2024: Heiko Engwerda (HE): included residual based FDE and
% protection level computation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  glc gls %RESIDUALS_PPP motionV2 STATE EPH
persistent II 

FlagBRD2FINE=glc.FlagBRD2FINE_PPP;

if isempty(gls.Log.PPP)
    gls.Log.PPP.RESIDUALS=NaN(glc.MAXSAT*2*3,glc.Alloc);
    gls.Log.PPP.STATE=NaN(rtk.nx,glc.Alloc);
    gls.Log.PPP.EPH=NaN(glc.MAXSAT*3,glc.Alloc);
    gls.Log.PPP.Nsat=NaN(1,glc.Alloc);
    gls.Log.PPP.Tropo=NaN(glc.MAXSAT,glc.Alloc);
    gls.Log.PPP.Iono=NaN(glc.MAXSAT,glc.Alloc);
    gls.Log.PPP.time=NaN(1,glc.Alloc);
    gls.Log.PPP.sim=NaN(1,glc.Alloc);

end
if isempty(II)
    II=1;
end

nobs=size(obs,1); MAXITER=8; iter=1; stat=0;
exc=zeros(nobs,1); dr=zeros(3,1);

% initialize rtk.sat.fix
for i=1:glc.MAXSAT
    for j=1:rtk.opt.nf
        rtk.sat(i).fix(j)=0;
    end
end

% debug tracing
% fprintf('\n \n \n');
% fprintf('time= %d',obs(1).time.time);fprintf('\n');
% fprintf('before time update');
% printx_ppp(rtk.x,rtk);
% printP(rtk.P,rtk);

% time update
rtk=udstate_ppp(rtk,obs,nav);

% debug tracing
% fprintf('after time update');
% printx_ppp(rtk.x,rtk);
% printP(rtk.P,rtk);

% cumpute satellite position,clock bias,velocity,clock drift
sv=satposs(obs,nav,rtk.opt.sateph,rtk,FlagBRD2FINE);

% exclude measurements of eclipsing satellite (block IIA)
% if rtk.opt.posopt(4)==1
%     sv=testeclipse(obs,nav,sv);
% end

% if opt.LEO_Aug==0
%     for ii=1:length(sv)
%         rt=norm(sv(ii).pos)-glc.RE_WGS84;
%         if abs(rt)<1e7
%             sv(ii).pos=zeros(3,1);
%             sv(ii).svh=-1;
%         end
%     end
% end

% compute earth tides correction
if rtk.opt.tidecorr==1
    dr=tidedisp(gpst2utc(obs(1).time),rtk.sol.pos,7,nav.erp,nav.otlp);
end

% HE: initialize sliding window buffers
sliding_window = 120; % sliding window size ideally according to correlation time (2 minutes for instance)
Buffer_V_inv = zeros(50,50,sliding_window); % max 50 satellites and sliding window of "sliding_window" length
Buffer_y = zeros(50,sliding_window); % max 50 satellites and sliding window of "sliding_window" length
PRN_list = [];

while iter<=MAXITER
    
    % calculate RESIDUALS_PPP,measurement matrix,measurement noise matrix
    [v,H,R,~,exc,nv,rtk,PRN]=ppp_res(0,rtk.x,rtk,obs,nav,sv,dr,exc);
    if nv==0,break;end
    
    
    % debug tracing
%     printv(v);
%     printH(H);
%     printR(R);

    % measurement update
    [X,P,stat_tmp]=Kfilter_h(v,H,R,rtk.x,rtk.P);
    if stat_tmp==0
        [week,sow]=time2gpst(obs(1).time);
        fprintf('Warning:GPS week = %d sow = %.3f,filter error!\n',week,sow);
        stat=0;break;
    end

    % Edit Heiko (HE) -------------------------------------------------
    [PL1,PL2,PL3,PL4,PRN_list] = calculate_integrity(X,H,P,R,v,Buffer_V_inv,Buffer_y,PRN(PRN>0),PRN_list);
    gls.Log.PPP.PL = PL4;
    % End Edit Heiko (HE) ----------------------------------------------

    % debug tracing
%     fprintf('after measurement update');
%     printx_ppp(X,rtk);
%     printP(P,rtk);
    
    % calculate posteriori RESIDUALS_PPP,validate the solution
    [v_post,H_post,R_post,azel_post,exc,stat,rtk,PRN,v_tot]=ppp_res(iter,X,rtk,obs,nav,sv,dr,exc);
         if stat_tmp==1
            gls.Log.PPP.STATE(1:length(X),II)=X;
            Sat_vis=vertcat(obs.sat);
            Appo=reshape(vertcat(sv.pos),3,length(Sat_vis));
            gls.Log.PPP.EPH(Sat_vis,II)=Appo(1,:).';
            gls.Log.PPP.EPH(Sat_vis+glc.MAXSAT+1,II)=Appo(2,:).';
            gls.Log.PPP.EPH(Sat_vis+glc.MAXSAT+glc.MAXSAT+1,II)=Appo(3,:).';
%           lll=vertcat(obs.sat);
%             PRN=PRN(PRN~=0);
            gls.Log.PPP.Nsat(1,II)=round(sum(v_tot~=0)/4);
            gls.Log.PPP.RESIDUALS_PPP(PRN(PRN~=0),II)=v_tot(PRN~=0);
            gls.Log.PPP.time(1,II)=obs(1).time.time;
            gls.Log.PPP.sim(1,II)=II;
%             RESIDUALS_PPP(90+PRN,II)=v_post(2:2:2*length(PRN));
%             if opt.ionoopt~=2
%                 switch opt.nf
%                     case 1
%                     RESIDUALS_PPP(PRN,II)=v_post(1:2:2*length(PRN));
%                     RESIDUALS_PPP(90+PRN,II)=v_post(2:2:2*length(PRN));
%                     case 2
%                     RESIDUALS_PPP(PRN,II)=v_post(1:2:2*length(PRN));
%                     RESIDUALS_PPP(90+PRN,II)=v_post(2:2:2*length(PRN));
%                     RESIDUALS_PPP(180+PRN,II)=v_post(2*length(PRN)+1:2:4*length(PRN));
%                     RESIDUALS_PPP(270+PRN,II)=v_post(2*length(PRN)+1:2:4*length(PRN));
%                     case 3
%                     RESIDUALS_PPP(PRN,II)=v_post(1:2:2*length(PRN));
%                     RESIDUALS_PPP(90+PRN,II)=v_post(2:2:2*length(PRN));
%                     RESIDUALS_PPP(180+PRN,II)=v_post(2*length(PRN)+1:2:4*length(PRN));
%                     RESIDUALS_PPP(270+PRN,II)=v_post(2*length(PRN)+1:2:4*length(PRN));
%                     RESIDUALS_PPP(360+PRN,II)=v_post(4*length(PRN)+1:2:6*length(PRN));
%                     RESIDUALS_PPP(450+PRN,II)=v_post(4*length(PRN)+1:2:6*length(PRN));                  
%                     
%                 end
% 
%             end
            II=II+1;
         end




    iter=iter+1;
    
    if stat==1
        rtk.x=X; rtk.P=P;
        break;
    end
    
end

if iter>MAXITER
    [week,sow]=time2gpst(obs(1).time);
    fprintf('Warning:GPS week=%d sow=%.3f, ppp iteration overflows',week,sow);
end

if stat==1
    % PPP ambiguity resolution (Not supported for the time being)
    
    % update solution
    rtk=update_stat(rtk,obs,glc.SOLQ_PPP);

    % hold fixed ambiguity (Not supported for the time being)
    
end

return

