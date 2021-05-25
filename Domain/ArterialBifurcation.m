%----------------------------------------------------------- MAIN FUNCTIONS
function [x] = ArterialBifurcation(Demand,Arg)
 BdBox = [0 14/10 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds();
    case('BdBox'); x = BdBox;
    case('FixedPoints');  x = FixedPoints(Arg{:});
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
d1 = dRectangle(P,2/10,12/10,0,1);
d2 = dRectangle(P,0,2/10,1/6,2/6);
d3 = dRectangle(P,0,2/10,4/6,5/6);
d4 = dRectangle(P,12/10,14/10,1/6,2/6);
Dist = dUnion(dUnion(dUnion(d1,d2),d3),d4);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds()
  x = {[],[]};
%------------------------------------------------------------- FIXED POINTS
function [x] = FixedPoints(P,R_P)
PFix = [  2/10, 1/6, 3;
          2/10, 2/6, 2;
          2/10, 4/6, 3;
          2/10, 5/6, 2;
         12/10, 1/6, 4;
         12/10, 2/6, 1];
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