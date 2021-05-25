function [wQ,xQ,N_p,S,S_grad,P_grad] = ShapeFncVEM(vx,vy)
%---------Forming projections for linear VEM --------------------------%
% Author: HENG CHI

% Input: vx: vector of x cooridnates of the vertices
%        vy: vector of y coordinates of the vertices          
%-------------------------------------------------------------------------%

%-----forming an integration rule on E that is exact for constant functions
[xC,A] = PolyCentroid(vx,vy);
xQ=xC;
wQ=A;
h_E = sqrt(A);
%-----forming the M,R,S matrices for projection \Pi_E^0
 M=A.*eye(2);
 nn = length(vx);
 Enodelist=[vx,vy;vx(1),vy(1)];
 dL=Enodelist(2:end,:)-Enodelist(1:end-1,:);
 n=[dL(:,2) -dL(:,1)];
 IntX=zeros(nn,1);
 IntY=zeros(nn,1);
 for edge=1:nn
   if edge==nn
      IntX(edge,:)=IntX(edge,:)+n(edge,1).*1/2;
      IntY(edge,:)=IntY(edge,:)+n(edge,2).*1/2;
      IntX(1,:)=IntX(1,:)+n(edge,1).*1/2;
      IntY(1,:)=IntY(1,:)+n(edge,2).*1/2;
   else
       IntX(edge,:)=IntX(edge,:)+n(edge,1).*1/2;
       IntY(edge,:)=IntY(edge,:)+n(edge,2).*1/2;
       IntX(edge+1,:)=IntX(edge+1,:)+n(edge,1).*1/2;
       IntY(edge+1,:)=IntY(edge+1,:)+n(edge,2).*1/2;
   end    
  end 
  R=[IntX,IntY];
  S=R/M;    
%-----forming the M^{nabla}, R^{nabla} and S^{nabla} matrices for projection \Pi_E^{nabla}    
 temp=[nn;(sum(vx)-nn*xC(:,1))/h_E;(sum(vy)-nn*xC(:,2))/h_E];
 M_grad=[temp,[0,0;eye(2)]];
 R_grad=[ones(nn,1),1/h_E*IntX,1/h_E*IntY];
 S_grad=R_grad/M_grad;
%-----forming the P^{nabla} matrix for projection \Pi_E^{nabla}      
 G_grad=[ones(1,nn);(vx-xC(:,1))'/h_E;(vy-xC(:,2))'/h_E];
 P_grad=S_grad*G_grad;
%---------------------------------------------------------------
% N: evaluated in the center of the element
m_k = BasisPk(xQ,1,xC,h_E);
N_p=zeros(nn,1,length(wQ));
for q=1:length(wQ)
    N_p(:,1,1)=S_grad*m_k(:,q);
end
% N: evaluated in the coordinates of the node of the element
% m_k = BasisPk([vx vy],1,xC,h_E);
% N_v=zeros(nn,1,nn);
% for q=1:nn
%     N_v(:,1,q)=S_grad*m_k(:,q);
% end
%-----------------------------------------------------------------
function m_k = BasisPk(X,k,Xc,h)
nQua = size(X,1);
x=(X(:,1)-Xc(1))/h; y=(X(:,2)-Xc(2))/h;
switch k
    case(1)
        m_k = zeros(3,nQua);
        for q=1:nQua
            m_k(:,q) = [1,x(q),y(q)];
        end
    case(2) 
        m_k = zeros(6,nQua);
        for q=1:nQua
            m_k(:,q) = [1,x(q),y(q),x(q)^2,x(q)*y(q),y(q)^2];
        end
end

function [xC,A] = PolyCentroid(vx,vy)
nv = length(vx);
vxS=vx([2:nv 1]); vyS=vy([2:nv 1]); %Shifted vertices
temp = vx.*vyS - vy.*vxS; 
A = 0.5*sum(temp);
xC = 1/(6*A)*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
