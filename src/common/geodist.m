function [r,LOS]=geodist(rs,rr)

global glc OptionSatNavEph;
r=-1; LOS=zeros(1,3);

if norm(rs)<glc.RE_WGS84, return; end

delta=rs-rr;
r=norm(delta);
LOS=delta'/r;
if ~OptionSatNavEph
r=r+glc.OMGE*(rs(1)*rr(2)-rs(2)*rr(1))/glc.CLIGHT;
end

return