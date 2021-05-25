function [vem]= GraftArteryBoundCond(vem,fluid)
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
 %  Nodal Boundary condition Wall| x= 0, y<1/6    u = u(y) ,v =0
if  vem.Node(nodebc1,1) <= 0.01 &&  vem.Node(nodebc2,1) <= 0.01  && ... 
    vem.Node(nodebc1,2) <= 1/6  &&  vem.Node(nodebc2,2) <= 1/6               
    k = k + 1; 
    kn = kn+1;
    vem.NodeBC(edofbc,1) =1;
    vem.NodeBC(edofbc,2) =1;
    yy = vem.Node(nodebc1,2);
    u = 24*Umax2*(yy)*(1-6*yy);
    vem.FBC(nodebc1,1) = u;
    yy = vem.Node(nodebc2,2);
    u = 24*Umax2*(yy)*(1-6*yy);
    vem.FBC(nodebc2,1) = u;
    vem.FBC(nodebc2,2) =0.0; 
    %  intlet
    NeumanBC(kn,1) = nele;
    NeumanBC(kn,2) = side;
    hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'ro');
    hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'ro');
else
  %  Nodal Boundary condition Wall| x= 0, 1/2<y<2/3    u = u(y) ,v =0
if  vem.Node(nodebc1,1) <= 0.01 &&  vem.Node(nodebc2,1) <= 0.01  && ... 
    vem.Node(nodebc1,2) >= 1/2  &&  vem.Node(nodebc2,2) >= 1/2  && ...
    vem.Node(nodebc1,2) <= 2/3  &&  vem.Node(nodebc2,2) <= 2/3               
    k = k + 1; 
    kn = kn+1;
    vem.NodeBC(edofbc,1) =1;
    vem.NodeBC(edofbc,2) =1;
    yy = vem.Node(nodebc1,2);
    u = -24*Umax1*(1-2*yy)*(2-3*yy);
    vem.FBC(nodebc1,1) = u;
    yy = vem.Node(nodebc2,2);
    u = -24*Umax1*(1-2*yy)*(2-3*yy);
    vem.FBC(nodebc2,1) = u;
    vem.FBC(nodebc2,2) =0.0; 
    %  intlet
    NeumanBC(kn,1) = nele;
    NeumanBC(kn,2) = side;
    hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'bo');
    hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'bo');
else
  %  Nodal Boundary condition Wall u = v = 0 
if  vem.Node(nodebc1,1) <= 0.209 &&  vem.Node(nodebc2,1) <= 0.209  && ... 
    vem.Node(nodebc1,2) >= 0.159 &&  vem.Node(nodebc2,2) >= 0.159  && ...
    vem.Node(nodebc1,2) <= 0.509 &&  vem.Node(nodebc2,2) <= 0.509  || ...
    vem.Node(nodebc1,1) <= 0.209 &&  vem.Node(nodebc2,1) <= 0.209  && ... 
    vem.Node(nodebc1,2) >= 0.659 &&  vem.Node(nodebc2,2) >= 0.659  || ...
    vem.Node(nodebc1,1) >= 1.199 &&  vem.Node(nodebc2,1) >= 1.199  && ...
    vem.Node(nodebc1,2) >= 0.159 &&  vem.Node(nodebc2,2) >= 0.159        
    k = k + 1;
    vem.NodeBC(edofbc,1) =1;
    vem.NodeBC(edofbc,2) =1;
    vem.FBC(edofbc,1) =0.0;
    vem.FBC(edofbc,2) =0.0;
    hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'mo');
    hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'mo');
else    
 %  Nodal Boundary condition Wall| y =0 ==> u = v = 0 
if  vem.Node(nodebc1,2) <= 0.01 &&  vem.Node(nodebc2,2) <= 0.01
    k = k + 1;
    vem.NodeBC(edofbc,1) =1;
    vem.NodeBC(edofbc,2) =1;
    vem.FBC(edofbc,1) =0.0;
    vem.FBC(edofbc,2) =0.0;
    hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'ko');
    hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'ko');
else
 %  Nodal Boundary condition Wall| y = 5/6 ==> u = v = 0 
if  vem.Node(nodebc1,2) >= 0.829 &&  vem.Node(nodebc2,2) >= 0.829
    k = k + 1;
    vem.NodeBC(edofbc,1) =1;
    vem.NodeBC(edofbc,2) =1;
    vem.FBC(edofbc,1) =0.0;
    vem.FBC(edofbc,2) =0.0;
    hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'go');
    hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'go');
else
  %  Nodal Boundary condition Wall| x = 1, outlet pressure
  
if  vem.Node(nodebc1,1) >= 1.39 &&  vem.Node(nodebc2,1) >= 1.39  && ...
    vem.Node(nodebc1,2) <= 1/6  &&  vem.Node(nodebc2,2) <= 1/6             
    kn = kn+1;
    vem.ElemFBC(ielebc,1) =0.0;   
    NeumanBC(kn,1) = nele;
    NeumanBC(kn,2) = side;
    NeumanBC(kn,3) = 1;
    hold on;plot(vem.Node(nodebc1,1),vem.Node(nodebc1,2),'yo');
    hold on;plot(vem.Node(nodebc2,1),vem.Node(nodebc2,2),'yo');
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
  