

function [f,dfdE] = SensitivityElemOpt(vem,fluid,C)  
f = 0; 
dfdE = zeros(length(vem.ElemOpt),1); 
visc_nn = fluid.visc_nn;
visc   = fluid.mi;
for iele= 1:length(vem.ElemOpt)
     
kappal = vem.kappa(vem.ElemOpt(iele)); 
NNodeElem = size(vem.Element{vem.ElemOpt(iele)},2);
velocity = zeros(2,NNodeElem);
for ilnode=1:NNodeElem
    velocity(1,ilnode) = C(vem.ASSMtrx(ilnode,vem.ElemOpt(iele)));
    velocity(2,ilnode) = C(vem.ASSMtrx(ilnode+NNodeElem,vem.ElemOpt(iele)));
end

u1 =velocity(1,:)';
u2 =velocity(2,:)';

U = [u1;
     u2]; 

S = vem.ShapeFnc{vem.ElemOpt(iele)}.S;

dudx = S(:,1)'*u1 ; 
dudy = S(:,2)'*u1 ;
dvdx = S(:,1)'*u2 ; 
dvdy = S(:,2)'*u2 ;

%     
 if visc_nn == 0  
      eta = visc; 
 else 
      %  ****  Viscosity: Carreau-Yasuda model   ****
        [eta, ~] = Visc_Carreau(dudx,dudy,dvdx,dvdy);
 end

  [Kk,Kd] = Mat_kk_kd(vem.ElemOpt(iele),vem,eta);
  
   % Computing objective function (f) and sensitivity (dfdE)
  f = f + 0.5*U'*(Kd + kappal*Kk)*U; 
  
  dfdE(iele,1) = 0.5*vem.dkappa(vem.ElemOpt(iele))*U'*Kk*U; 
  
end

function [Kk,Kd] = Mat_kk_kd(iele,vem,eta)
%
NNodeElem = size(vem.Element{iele},2);

w0 = zeros(NNodeElem,NNodeElem);


 S = vem.ShapeFnc{iele}.S;
 N_p = vem.ShapeFnc{iele}.N_p;
 A =   vem.ShapeFnc{iele}.W;
 P_grad = vem.ShapeFnc{iele}.P_grad;
 
 Mm = N_p*N_p'*A ; 
 Kk = [ Mm    w0;
        w0    Mm];
   
 eDof = vem.Element{iele}(:); nn = size(eDof,1);
 m = A.*eye(3);
 Cmu = eta*[2 0 0; 0 2 0; 0 0 1];
 s = zeros(3,2*nn);
 s(1,   1:1:nn)   = S(:,1)'; 
 s(2,nn+1:1:2*nn) = S(:,2)'; 
 s(3,   1:1:nn)   = S(:,2)'; 
 s(3,nn+1:1:2*nn) = S(:,1)';

 pgrad = [P_grad zeros(nn); zeros(nn) P_grad];

 ss = s'*m*Cmu*s;
 Ip = (eye(2*nn)-pgrad)*(eye(2*nn)-pgrad');
 Kd = ss  + Ip; fact = 1.0;
 Kd = fact*(1/2)*(Kd + Kd');




   

