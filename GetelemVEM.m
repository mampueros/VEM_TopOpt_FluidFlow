function [elemJ,elemRv,vem] = GetelemVEM(iele,vem,fluid,C,kappa)
%
eDof = vem.Element{iele}(:); nn = size(eDof,1);

Rho = fluid.ro;
NNodeElem = size(vem.Element{iele},2);
ElemNDoFPress = 1;

iuldof = zeros(1,NNodeElem); ivldof = zeros(size(iuldof));
for i = 1:NNodeElem
   iuldof(i) = i;
   ivldof(i) = i+NNodeElem;
end

ipldof = zeros(1,ElemNDoFPress);
for i = 1:ElemNDoFPress
   ipldof(i) = i+2*NNodeElem;
end
velocity = zeros(2,NNodeElem);
for ilnode=1:NNodeElem
    velocity(1,ilnode) = C(vem.ASSMtrx(ilnode,iele));
    velocity(2,ilnode) = C(vem.ASSMtrx(ilnode+NNodeElem,iele));
end

pressure = zeros(1,ElemNDoFPress);
for ipress = 1:ElemNDoFPress
    pressure(ipress) = C(vem.ASSMtrx(2*NNodeElem+ipress,iele)); 
end

visc_nn = fluid.visc_nn;
visc   = fluid.mi;
 
w0 = zeros(NNodeElem,NNodeElem);
u1 =velocity(1,:)';
u2 =velocity(2,:)';

U = [u1;
     u2 ;
    pressure]; 
%%
 S = vem.ShapeFnc{iele}.S;
 N_p = vem.ShapeFnc{iele}.N_p;
 A =   vem.ShapeFnc{iele}.W;
 P_grad = vem.ShapeFnc{iele}.P_grad;
 S_grad = vem.ShapeFnc{iele}.S_grad;

dudx = S(:,1)'*u1 ; 
dudy = S(:,2)'*u1 ;
dvdx = S(:,1)'*u2 ; 
dvdy = S(:,2)'*u2 ;

 if visc_nn == 0  
      eta = visc; d_eta = 0;
 else 
      %  ****  Viscosity: Carreau-Yasuda model  ****
        [eta, d_eta,vem] = Visc_Carreau(dudx,dudy,dvdx,dvdy,vem);
         
        vem.Gamma2(iele,:) = vem.gamma2; 
 end


 uu = N_p'*u1 ;
 vv = N_p'*u2 ;
 
 % Local matrix Kk
 Mm = N_p*N_p'*A ; 
 
 Kk = kappa*[ Mm    w0;
              w0    Mm];
          
 % Local matrix Kd
         
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
 
% Assembling: linear

% Local matrix: C
 ci = Rho*N_p*A*(uu*S(:,1)' + vv*S(:,2)') ;

 C = [ ci    w0;
       w0    ci];

% Local vector: Q
 h_E = sqrt(A);
 Qu = h_E*S_grad(:,2); 
 Qv = h_E*S_grad(:,3);

 Q = [ Qu ;
       Qv ] ;
   

  matl = [C + Kd + Kk   -Q;
            -Q'          0];
      
% Assembling: non-linear

  C11 = A*Rho*N_p*dudx*N_p';
  C12 = A*Rho*N_p*dudy*N_p';
  C21 = A*Rho*N_p*dvdx*N_p';
  C22 = A*Rho*N_p*dvdy*N_p';
   
  Pw0 = zeros(NNodeElem,1);

  Cnl = [ C11  C12  Pw0;
          C21  C22  Pw0;
          Pw0' Pw0'  0];
               
 % Mg
 
 Mg11 = 2*S(:,1)*S(:,1)' + S(:,2)*S(:,2)';
 Mg12 = S(:,2)*S(:,1)';
 Mg21 = S(:,1)*S(:,2)';
 Mg22 = 2*S(:,2)*S(:,2)' + S(:,1)*S(:,1)';
 
 
 Mg = [Mg11 Mg12; 
       Mg21  Mg22];
   
 Ueta = [u1; 
         u2];
 
 term = d_eta*ones(2*NNodeElem,1)'*Mg*Ueta*Kd;
 PW0 = zeros(2*NNodeElem,1);
 
 term2 = [term  PW0;
           PW0'  0];
 
 matnl = Cnl + term2;
         
% Computing Jacobian matrix

 elemJ = matl + matnl;
 
 % Computing residual vector
 
 elemRv = matl*U;
 

 

     