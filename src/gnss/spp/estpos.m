function [rtk,sat_,stat]=estpos(rtk,obs,nav,sv,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimate reciever position and clock bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global  glc gls %RESIDUALS motionV2 EPH
persistent II 
if isempty(gls.Log.SPP)
    gls.Log.SPP.RESIDUALS=NaN(glc.MAXSAT*2*3,glc.Alloc);
    gls.Log.SPP.STATE=NaN(rtk.nx,glc.Alloc);
    gls.Log.SPP.EPH=NaN(glc.MAXSAT*3,glc.Alloc);
    gls.Log.SPP.Nsat=NaN(1,glc.Alloc);
    gls.Log.SPP.Tropo=NaN(glc.MAXSAT,glc.Alloc);
    gls.Log.SPP.Iono=NaN(glc.MAXSAT,glc.Alloc);
    gls.Log.SPP.time=NaN(1,glc.Alloc);
    gls.Log.SPP.sim=NaN(1,glc.Alloc);
end
if isempty(II)
    II=1;
end

NX=3+glc.NSYS; MAXITER=10; iter=1;
time=obs(1).time; xr0=[rtk.sol.pos';zeros(glc.NSYS,1)];



while iter<=MAXITER
    
    % compute residual,measurement model,weight matrix
    [v,H,P,vsat,azel,resp,nv,ns,tropoerr_vect,ionoerr_vect]=rescode(iter,obs,nav,sv,xr0,opt);
%     if obs(1).time.time>1.687910420000000e+09+200 && opt.LEO_Aug %vsat(1)==1 && opt.LEO_Aug
% 
%         disp('first_LEO');
% 
%     end

    sat_.vsat=vsat; 
    sat_.azel=azel;
    sat_.resp=resp;
    
    if nv<NX, stat=0;return;end
    
    % least square algrithm
    [dx,Q]=least_square(v,H,P);
    
    xr0=xr0+dx;
    iter=iter+1;
    
    if dot(dx,dx)<1e-4
        rtk.sol.time=timeadd(obs(1).time,-xr0(4)/glc.CLIGHT);
        rtk.sol.ns  =ns;
        rtk.sol.ratio=0;
        rtk.sol.pos =xr0(1:3)';
        rtk.sol.vel =zeros(1,3);
        rtk.sol.posP(1)=Q(1,1);rtk.sol.posP(2)=Q(2,2);rtk.sol.posP(3)=Q(3,3);
        rtk.sol.posP(4)=Q(1,2);rtk.sol.posP(5)=Q(2,3);rtk.sol.posP(6)=Q(1,3);
        rtk.sol.dtr(1) =xr0(4)/glc.CLIGHT;
        rtk.sol.dtr(2) =xr0(5)/glc.CLIGHT;
        rtk.sol.dtr(3) =xr0(6)/glc.CLIGHT;
        rtk.sol.dtr(4) =xr0(7)/glc.CLIGHT;
        rtk.sol.dtr(5) =xr0(8)/glc.CLIGHT;
          
        % validate solution
        if opt.L5==0
        stat=valsol(time,v,P,vsat,azel,opt);
        else
        stat=1;
        end
        if stat==1
%             lll=vertcat(obs.sat);
%             RESIDUALS(lll,II)=resp;
%             Sat_vis=vertcat(obs.sat);
%             Appo=reshape(vertcat(sv.pos),3,length(Sat_vis));
%             EPH(Sat_vis,II)=Appo(1,:).';
%             EPH(Sat_vis+90,II)=Appo(2,:).';
%             EPH(Sat_vis+190,II)=Appo(3,:).';
%             II=II+1;
%             rtk.sol.stat=glc.SOLQ_SPP; 

            gls.Log.SPP.STATE(1:length(xr0),II)=xr0;
            Sat_vis=vertcat(obs.sat);
            Appo=reshape(vertcat(sv.pos),3,length(Sat_vis));
            gls.Log.SPP.EPH(Sat_vis,II)=Appo(1,:).';
            gls.Log.SPP.EPH(Sat_vis+glc.MAXSAT+1,II)=Appo(2,:).';
            gls.Log.SPP.EPH(Sat_vis+glc.MAXSAT+glc.MAXSAT+1,II)=Appo(3,:).';
%           lll=vertcat(obs.sat);
%             PRN=PRN(PRN~=0);
            gls.Log.SPP.Nsat(1,II)=round(sum(vsat~=0)/2);
            gls.Log.SPP.RESIDUALS(Sat_vis(vsat~=0),II)=v(1:sum((vsat~=0)));
            gls.Log.SPP.tropoerr_vect(Sat_vis(vsat~=0),II)=tropoerr_vect(vsat~=0);
            gls.Log.SPP.ionoerr_vect(Sat_vis(vsat~=0),II)=ionoerr_vect(vsat~=0);
            gls.Log.SPP.time(1,II)=obs(1).time.time;
            gls.Log.SPP.sim(1,II)=II;
            II=II+1;
            rtk.sol.stat=glc.SOLQ_SPP;
            return;
        end
        return;
    end  
    
end

if iter>MAXITER
    stat=0;
    [week,sow]=time2gpst(time);
    fprintf('Warning:GPS week = %d sow = %.3f,SPP iteration divergent!\n',week,sow);
    return;
end

return



