function rtk=udsol_rtkins(rtk,stat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update solution status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global glc
ins=rtk.ins;

if stat==glc.SOLQ_FIX
    xa=ins.xa; Pa=ins.Pa; 
    [rr,Cen]=blh2xyz(xa(7:9)); T=Dblh2Dxyz(xa(7:9));
    rtk.sol.pos=rr';
    rtk.sol.vel=(Cen*xa(4:6))';
    rtk.sol.att=xa(1:3)'/glc.D2R;
    
    if rtk.sol.att(3)>=0
        rtk.sol.att(3)=360-rtk.sol.att(3);
    else
        rtk.sol.att(3)=-rtk.sol.att(3);
    end
    
    posvar=Pa(7:9,7:9); posvar=T*posvar*T';        %blh_var to xyz_var
    velvar=Pa(4:6,4:6); velvar=Cen*velvar*Cen';    %enu_var to xyz_var
    attvar=Pa(1:3,1:3); attvar=attvar./glc.D2R^2;  %rad_var to deg_var
    
    rtk.sol.posP(1)=posvar(1,1);
    rtk.sol.posP(2)=posvar(2,2);
    rtk.sol.posP(3)=posvar(3,3);
    rtk.sol.posP(4)=posvar(1,2);
    rtk.sol.posP(5)=posvar(2,3);
    rtk.sol.posP(6)=posvar(1,3);
    
    rtk.sol.velP(1)=velvar(1,1);
    rtk.sol.velP(2)=velvar(2,2);
    rtk.sol.velP(3)=velvar(3,3);
    rtk.sol.velP(4)=velvar(1,2);
    rtk.sol.velP(5)=velvar(2,3);
    rtk.sol.velP(6)=velvar(1,3);
    
    rtk.sol.attP(1)=attvar(1,1);
    rtk.sol.attP(2)=attvar(2,2);
    rtk.sol.attP(3)=attvar(3,3);
    rtk.sol.attP(4)=attvar(1,2);
    rtk.sol.attP(5)=attvar(2,3);
    rtk.sol.attP(6)=attvar(1,3);
else
    x=ins.x; P=ins.P;
    rtk.nfix=0;
    [rr,Cen]=blh2xyz(x(7:9)); T=Dblh2Dxyz(x(7:9));
    rtk.sol.pos=rr';
    rtk.sol.vel=(Cen*x(4:6))';
    rtk.sol.att=x(1:3)'/glc.D2R;
    
    if rtk.sol.att(3)>=0
        rtk.sol.att(3)=360-rtk.sol.att(3);
    else
        rtk.sol.att(3)=-rtk.sol.att(3);
    end
    
    posvar=P(7:9,7:9); posvar=T*posvar*T';        %blh_var to xyz_var
    velvar=P(4:6,4:6); velvar=Cen*velvar*Cen';    %enu_var to xyz_var
    attvar=P(1:3,1:3); attvar=attvar./glc.D2R^2;  %rad_var to deg_var
    
    rtk.sol.posP(1)=posvar(1,1);
    rtk.sol.posP(2)=posvar(2,2);
    rtk.sol.posP(3)=posvar(3,3);
    rtk.sol.posP(4)=posvar(1,2);
    rtk.sol.posP(5)=posvar(2,3);
    rtk.sol.posP(6)=posvar(1,3);
    
    rtk.sol.velP(1)=velvar(1,1);
    rtk.sol.velP(2)=velvar(2,2);
    rtk.sol.velP(3)=velvar(3,3);
    rtk.sol.velP(4)=velvar(1,2);
    rtk.sol.velP(5)=velvar(2,3);
    rtk.sol.velP(6)=velvar(1,3);
    
    rtk.sol.attP(1)=attvar(1,1);
    rtk.sol.attP(2)=attvar(2,2);
    rtk.sol.attP(3)=attvar(3,3);
    rtk.sol.attP(4)=attvar(1,2);
    rtk.sol.attP(5)=attvar(2,3);
    rtk.sol.attP(6)=attvar(1,3);
end

return

