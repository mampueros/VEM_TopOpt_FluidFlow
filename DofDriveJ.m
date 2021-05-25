function [vem] = DofDriveJ(vem)
%                                                                     
%.... program to establish equation numbers   		                    
% 
%   output
%         GlobalDoFID :    nodal dof 
%         ASSMtrx     :    element dof
%

MaxElemNDoF = 2*max(cellfun(@numel,vem.Element)) + 1;   MaxNodeNDoF = 2;
vem.ElemPressDof = zeros(vem.NElem,1);

vem.ASSMtrx = zeros(MaxElemNDoF,vem.NElem);
vem.GlobalDoFID = zeros(MaxNodeNDoF,vem.NNode);
vem.ElemNDof = 2*cellfun(@numel,vem.Element)+1;
%
  ieq=0;
  for iele =1:vem.NElem
      nl = size(vem.Element{iele},2);
     for k =1: MaxNodeNDoF       
       for  j =1:nl 
             ignode = vem.Element{iele}(j);
             ildof = vem.GlobalDoFID(k,ignode);
             if(ildof == 0) 
               ieq = ieq +1;
               vem.GlobalDoFID(k,ignode) = ieq;
             end  
       end % j
     end % k
%
     ieq = ieq +1;
     vem.ElemPressDof(iele) = ieq;
  end % iele
     NDoF = ieq;
     vem.NDoF  = NDoF ; 
 for iele=1:vem.NElem
     nl = size(vem.Element{iele},2);
     ip = 0 ;
     for k=1:MaxNodeNDoF   
        for j=1:nl 
            ip = ip+1;
            ignode = vem.Element{iele}(j);   
            ildof = vem.GlobalDoFID(k,ignode);
            vem.ASSMtrx(ip,iele) = ildof;
        end 
     end 
     ip = ip +1;
     vem.ASSMtrx(ip,iele) = vem.ElemPressDof(iele);
 end % nel
 %
 
  %  create sparse matrix 
  vem.i = zeros(sum(vem.ElemNDof.^2),1);
  vem.j = vem.i; vem.k = vem.i;
  index = 0;
  for iele= 1:vem.NElem
    NDofElem = vem.ElemNDof(iele);
    eDof = vem.ASSMtrx(1:NDofElem,iele);
    vem.eDof{iele} = eDof;
    II = repmat(eDof ,1,NDofElem);
    JJ = repmat(eDof',NDofElem,1);
    vem.i(index+1:index+NDofElem^2) = reshape(II ,NDofElem^2,1);
    vem.j(index+1:index+NDofElem^2) = reshape(JJ ,NDofElem^2,1);
    index = index + NDofElem^2;
  end
   