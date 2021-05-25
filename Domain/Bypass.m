%------------------------------ PolyTop ----------------------------------%
% Ref: 
%-------------------------------------------------------------------------%
function [x] = Bypass(Demand,Arg)
d = 0.6; 
L1 = 6; L2 = 3; L5 = L1 - L2 +(d/2);
H = 0.8;
BdBox = [0 2*L1 0 L5+H];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds();
    case('BdBox'); x = BdBox;
    case('FixedPoints');  x = FixedPoints(Arg{:});
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox) 
d = 0.6; 
L1 = 6; L2 = 3; L5 = L1 - L2 +(d/2); L3 = 1.9; L4 = d/4;
H = 0.8; e = 0.5;
d1 = dRectangle(P,0,L2+d+L3,0,H);
d2 = dRectangle(P,L1+e,2*L1,0,H);
d3 = dRectangle(P,L2,L2+d,H,L5+H);
d4 = dRectangle(P,2*L1-L2-d,2*L1-L2,H,H+L5);
d5 = dRectangle(P,L2+d,2*L1-L2-d,H+L4,H+L5);
Dist = dUnion(dUnion(dUnion(dUnion(d1,d2),d3),d4),d5);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds()
x = {[],[]};
%------------------------------------------------------------- FIXED POINTS
function [x] = FixedPoints(P,R_P)
d = 0.6; 
L1 = 6;L2 = 3; L5 = L1 - L2 +(d/2); L3 = 1.9;  L4 = d/4;
H = 0.8; e = 0.5;
PFix = [       L2, H   ,2;
          2*L1-L2, H   ,1;
             L2+d, H   ,1;
        2*L1-L2-d, H   ,2;
             L2+d, H+L4,4;
        2*L1-L2-d, H+L4,3];
switch nargin
  case 1
    PP = P;
  case 2
    PP = [P;R_P];
    PFix(:,3) = 0;
end
for i = 1:size(PFix,1)
  [B,I] = sort(sqrt((PP(:,1)-PFix(i,1)).^2+(PP(:,2)-PFix(i,2)).^2));
  if PFix(i,3)==0
    for j = 2:4
      n = PP(I(j),:) - PFix(i,1:2); n = n/norm(n);
      PP(I(j),:) = PP(I(j),:)-n*1/sqrt(2)*(B(j)-B(1));
    end
  elseif PFix(i,3)==1
    PP(I(1),:) = PFix(i,1:2)+1/sqrt(2)*[-B(1), B(1)];
    PP(I(2),:) = PFix(i,1:2)+1/sqrt(2)*[-B(1),-B(1)];
    PP(I(3),:) = PFix(i,1:2)+1/sqrt(2)*[ B(1),-B(1)];
  elseif PFix(i,3)==2
    PP(I(1),:) = PFix(i,1:2)+1/sqrt(2)*[ B(1), B(1)];
    PP(I(2),:) = PFix(i,1:2)+1/sqrt(2)*[-B(1),-B(1)];
    PP(I(3),:) = PFix(i,1:2)+1/sqrt(2)*[ B(1),-B(1)];
  elseif PFix(i,3)==3
    PP(I(1),:) = PFix(i,1:2)+1/sqrt(2)*[ B(1), B(1)];
    PP(I(2),:) = PFix(i,1:2)+1/sqrt(2)*[-B(1), B(1)];
    PP(I(3),:) = PFix(i,1:2)+1/sqrt(2)*[ B(1),-B(1)];
  elseif PFix(i,3)==4
    PP(I(1),:) = PFix(i,1:2)+1/sqrt(2)*[ B(1), B(1)];
    PP(I(2),:) = PFix(i,1:2)+1/sqrt(2)*[-B(1), B(1)];
    PP(I(3),:) = PFix(i,1:2)+1/sqrt(2)*[-B(1),-B(1)];
  end
end
P = PP(1:size(P,1),:); R_P = PP(1+size(P,1):end,:);
x = {P,R_P};
%-------------------------------------------------------------------%