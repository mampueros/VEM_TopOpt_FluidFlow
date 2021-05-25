function [C,vem] = VEM_NStokes(vem,fluid)
tol = 1e-6;
itermax = 30;
C = zeros(vem.NDoF,1);
%
iter = 0;
erro = 1.0 ; 
%
dat_err_vem = zeros(itermax,4);
while  ( (erro > tol) && (iter < itermax) )    
   iter = iter + 1;
   % Evaluate Jacobian Matrix and Residual Vector
   [J,vem,Rv] = FormJ(vem,fluid,C);  
    
   DC = J\-Rv;
   C = DC + C;
  
   rvnorm = norm(Rv);
   error = rvnorm;
   errod = norm(DC);
   erro = errod + error;
   
   dat_err_vem(iter,:) = [iter error errod erro]; 
   fprintf('iter:%i \t  erro: | dc| =%12.5e |R| =%12.5e |dc| +|R| =%12.5e\n',...
                      iter,errod,error,erro);
end
vem.dat_err_vem = dat_err_vem; 

if (rvnorm > tol)
    disp('Convergence Not Achieved ! ');
    return;
end

function [J,vem,Rv] = FormJ(vem,fluid,C)  
index = 0 ;
 Rv = zeros(vem.NDoF,1);
 
 vem.Gamma2 = zeros(vem.NElem,1);
for iele= 1:vem.NElem
  kappal = vem.kappa(iele);  
  [elemJ,elemRv,vem] = GetelemVEM(iele,vem,fluid,C,kappal);
  % Nodal boundary conditions
  if vem.NNodeBC > 0  
        [elemJ,elemRv] = NodalBC(iele,vem,elemJ,elemRv,C) ;
  end
  % Neumann boundary conditions
  if vem.NElemNBC > 0  
        [elemJ,elemRv] = NeumanBC(iele,vem,elemJ,elemRv);
  end
  NDofElem = vem.ElemNDof(iele);
  vem.k(index+1:index+NDofElem^2) = reshape(elemJ,NDofElem^2,1);
  index = index + NDofElem^2;
  % Form Rv (global)
  lm = vem.eDof{iele};  
  Rv(lm) = Rv(lm) + elemRv;
end

J = sparse(vem.i,vem.j,vem.k,vem.NDoF,vem.NDoF);


function [elemJ,elemRv] = NodalBC(iele,vem,elemJ,elemRv,C)  
% 
NNodeElem = size(vem.Element{iele},2);
 
nl = NNodeElem;
NDofElem = vem.ElemNDof(iele) ;
%  **** Velocity nodal boundary condition  ****  
%
    enode = vem.Element{iele}(:); 
    NodeBCl = vem.NodeBC(enode,1:2) ; 
    vmbclu = find(NodeBCl(:,1));
    ul = C(vem.GlobalDoFID(1,enode));
    vl = C(vem.GlobalDoFID(2,enode));
    if (vmbclu)        
       vmbclv = vmbclu + nl;
       fl = vem.FBC(enode,1:2);
       elemJ(vmbclu,1:NDofElem) = 0.0 ;
       elemJ(vmbclv,1:NDofElem) = 0.0 ; 
       elemRv(vmbclu) = ul(vmbclu) - fl(vmbclu,1) ;
       elemRv(vmbclv) = vl(vmbclu) - fl(vmbclu,2);  
       for i = 1:size(vmbclu)
          iu = vmbclu(i) ;
          iv = vmbclv(i) ; 
          elemJ(iu,iu) = 1.0 ;  
          elemJ(iv,iv) = 1.0 ; 
       end 
    end

function [elemJ,elemRv] = NeumanBC(iele,vem,elemJ,elemRv)  
% 
NDofElem = vem.ElemNDof(iele) ;
%  **** Pressure element boundary condition  ****
   ind = find(vem.NeumanBC(:,1)==iele);     
    if (ind)
        if vem.NeumanBC(ind,3) == 1
           flp = vem.ElemFBC(iele,1);
           elemJ(NDofElem,1:NDofElem) = 0.0 ;
           elemJ(NDofElem, NDofElem) = 1.0 ;
           elemRv(NDofElem) = flp ; 
        end
    end

   

