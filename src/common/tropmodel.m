function corr=tropmodel(el,pos,humi)

temp0=15.0; %temparature at sea level

if (pos(3)<-100 || pos(3)>5000 || el<=0) 
    corr=0;
else
    % standard atmosphere
    if pos(3)<0,hgt=0;else,hgt=pos(3);end
    pres=1013.25*(1.0-2.2557e-5*hgt)^5.2568;
    temp=temp0-6.5^-3*hgt+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));
    
    %saastamoninen model
    z=pi/2.0-el;
    trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos(1))-0.00028*hgt/1000)/cos(z);
    trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
    corr=trph+trpw;
end

return

