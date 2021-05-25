function [vem]=ByPassBoundCond(vem,fluid,BdBox)
A = BdBox(2)-BdBox(1); 
B = BdBox(4)-BdBox(3); 
l1 = 6; l2 = 3;  H = 0.8; e = 0.5;
Umax1 = fluid.Umax_1 ;
Umax2 = fluid.Umax_2 ;
NNode = vem.NNode; 
NElemBC = vem.NElemBC  ;
NElem = vem.NElem;
 NeumanBC = zeros(NElemBC,5); 
vem.NodeBC = zeros(NNode,2);
vem.FBC = zeros(NNode,2);
vem.ElemFBC = zeros(NElem,2);
k = 0;
kn = 0;
 for ielebc = 1: NElemBC
       nele = vem.ElementBCJ(ielebc,1);   
       NNodeElem = size(vem.Element{nele},2);  
       side = vem.ElementBCJ(ielebc,2); 
       nodel1 = side;  
       nodel2 = nodel1 + 1;
       if nodel2 > NNodeElem, nodel2 =1; end
        nodebc1 = vem.Element{nele}(nodel1); 
        nodebc2 = vem.Element{nele}(nodel2); 
        edofbc = [nodebc1; nodebc2] ; 
%  Nodal boundary condition Inlet u u(y), v =0
   if  vem.Node(nodebc1,1) <= 0.01 &&  vem.Node(nodebc2,1) <= 0.01  && ... 
       vem.Node(nodebc1,2) <= B  &&  vem.Node(nodebc2,2) <= B
            k = k + 1; 
            kn = kn+1;
            vem.NodeBC(edofbc,1) =1;
            vem.NodeBC(edofbc,2) =1;
            y1 = vem.Node(nodebc1,2);
            yy = y1 ; 
            u = (4/H^2)*Umax1*(yy)*(H-yy)  ; %hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'ro');
            vem.FBC(nodebc1,1) = u;
            y2 = vem.Node(nodebc2,2);
            yy = y2;
            u = (4/H^2)*Umax1*(yy)*(H-yy)  ; %hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'bo');
            vem.FBC(nodebc2,1) = u;
            vem.FBC(nodebc2,2) =0.0; 
      %  Inlet
            NeumanBC(kn,1) = nele;
            NeumanBC(kn,2) = side;
   else
       %  Nodal Boundary condition Wall u = v = 0 
        if   vem.Node(nodebc1,1) >= 0.01 &&  vem.Node(nodebc2,1) >= 0.01  && ... 
             vem.Node(nodebc1,1) <= 3.01 &&  vem.Node(nodebc2,1) <= 3.01  && ...
             vem.Node(nodebc1,2) >= 0.79 &&  vem.Node(nodebc2,2) >= 0.79  || ...
             vem.Node(nodebc1,1) >= A-l2 &&  vem.Node(nodebc2,1) >= A-l2  && ... 
             vem.Node(nodebc1,1) <= A    &&  vem.Node(nodebc2,1) <= A     && ...
             vem.Node(nodebc1,2) >= 0.79 &&  vem.Node(nodebc2,2) >= 0.79  || ...
             vem.Node(nodebc1,1) >= 3.59 &&  vem.Node(nodebc2,1) >= 3.59  && ...
             vem.Node(nodebc1,1) <= 8.41 &&  vem.Node(nodebc2,1) <= 8.41  && ...
             vem.Node(nodebc1,2) >= 0.79 &&  vem.Node(nodebc2,2) >= 0.79  && ...
             vem.Node(nodebc1,2) <= 0.96 &&  vem.Node(nodebc2,2) <= 0.96  || ...
             vem.Node(nodebc1,1) >= 5.49 &&  vem.Node(nodebc2,1) >= 5.49  && ...
             vem.Node(nodebc1,1) <= 6.51 &&  vem.Node(nodebc2,1) <= 6.51  && ...
             vem.Node(nodebc1,2) <= H &&  vem.Node(nodebc2,2) <= H
         
            
         
         hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'mo');
         hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'mo');
         
             k = k + 1;
             vem.NodeBC(edofbc,1) =1;
             vem.NodeBC(edofbc,2) =1;
             vem.FBC(edofbc,1) =0.0;
             vem.FBC(edofbc,2) =0.0;
        else
            %  Nodal Boundary condition Wall| y =0 ==> u = v = 0 
            if  vem.Node(nodebc1,2) <= 0.01 &&  vem.Node(nodebc2,2) <= 0.01
                     
                hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'go');
                hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'go');
                         
                 k = k + 1;
                 vem.NodeBC(edofbc,1) =1;
                 vem.NodeBC(edofbc,2) =1;
                 vem.FBC(edofbc,1) =0.0;
                 vem.FBC(edofbc,2) =0.0;
            else
                    %  Nodal Boundary condition Wall| y = B ==> u = v = 0 
               if  vem.Node(nodebc1,2) >= 4.09 &&  vem.Node(nodebc2,2) >= 4.09 
%                         
                   hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'ko');
                   hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'ko');
                        
                   k = k + 1;
                   vem.NodeBC(edofbc,1) =1;
                   vem.NodeBC(edofbc,2) =1;
                   vem.FBC(edofbc,1) =0.0;
                   vem.FBC(edofbc,2) =0.0;
               else      
                       %  Pressure out x=A, y<H , p = 0.0
                  if vem.Node(nodebc1,1) >= 11.9999 &&  vem.Node(nodebc2,1) >= 11.9999  
                     
                       hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'ro');
                         hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'ro');
                       kn = kn+1;
                       vem.ElemFBC(ielebc,1) =0.0;   
                       NeumanBC(kn,1) = nele; hold on;plot(vem.Node(vem.Element{nele},1),vem.Node(vem.Element{nele},2),'bo');
                       NeumanBC(kn,2) = side;
                      NeumanBC(kn,3) = 1;
                  end
               end
           end
       end
   end
 end
vem.NElemNBC =kn;
vem.NNodeBC = k;      
vem.NeumanBC = NeumanBC(1:kn,:);
end
  