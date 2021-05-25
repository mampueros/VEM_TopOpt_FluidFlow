%% Non newtonian Topology Optimization VEM
clear;clc;close all;
% 
fprintf('\n\n=========================================================================');
fprintf('\n    *** Topology Optimization for Non-Newtonian flow using VEM  *** ');
fprintf('\n=========================================================================\n');

%% Domains:
addpath(genpath('./PolyMesher'),genpath('./Domain'));
n = input('Enter a number, (0) By-Pass, (1) Arterial Bifurcation , (2) Graft stenosed Artery: \n');
switch n
    case 0
        disp('Graft design for a blocked artery')
        Domain = @Bypass;
        BdBox = Domain('BdBox');
        Nelements = 200;
        [Node,Element,~,~,P] = PolyMesher(Domain,Nelements,100);
        I = find( P(:,2)>0.8);
    case 1
        disp('Arterial bypass design for a stenosed artery')
        Domain = @ArterialBifurcation;
        BdBox = Domain('BdBox');
        Nelements =2000;
        [Node,Element,~,~,P] = PolyMesher(Domain,Nelements,30);
        I = find(P(:,1)>=2/10 & P(:,1)<=12/10); 
    case 2
        disp('Graft design for a stenosed artery') 
        Domain = @GraftArtery;
        BdBox = Domain('BdBox');
        Nelements = 2000;
        [Node,Element,~,~,P] = PolyMesher(Domain,Nelements,100);
        I = find(P(:,1)>= 2/10 & P(:,2)>=1/6); 
    otherwise
        disp('other value')
end
rmpath('./PolyMesher','./Domain');
% verifying elements selected to optimize...      
hold on;      
for i = 1:length(I)
    vert = Node(Element{I(i)},:);
    plot(vert(:,1),vert(:,2),'b*');
end
%  -------------------------------------------------- CREATE 'fluid' STRUCT
fluid = struct( ...
  'visc_nn',1,...   % Newtonian (0)/ Non-Newtonian fluid (1)
  'ro',1.058, ...   % Density
  'mi',0.0345, ... % Viscosity 
  'Umax_1',1.5, ... % max velocity inlet 01
  'Umax_2',1.5 ...  % max velocity inlet 02 (not in this case)
  ) ;
% ---------------------------------------------------- CREATE 'vem' STRUCT
%% construct vem struct vem 
vem = struct( ...
  'NNode',size(Node,1), ...     % Number of nodes
  'NElem',size(Element,1), ...  % Number of elements
  'Node',Node, ...              % [NNode x 2] array of nodes
  'ElemOpt',I,...               % Elements to optimize
  'Element',{Element} ...       % [NElement x Var] cell array of elements
   ); 
%% forming vem projections
   vem.ShapeFnc = cell(vem.NElem,1); 
   for el=1:vem.NElem
       ElementNode=vem.Element{el};
       vx=vem.Node(ElementNode,1);vy=vem.Node(ElementNode,2);
       [wQ,xQ,N_p,S,S_grad,P_grad] = ShapeFncVEM(vx,vy);
       vem.ShapeFnc{el}.W=wQ; 
       vem.ShapeFnc{el}.X=xQ;
       vem.ShapeFnc{el}.N_p=N_p;
       vem.ShapeFnc{el}.S=S;
       vem.ShapeFnc{el}.S_grad=S_grad;
       vem.ShapeFnc{el}.P_grad=P_grad;
   end 
%%
Nmax = max(cellfun(@length,vem.Element));
[vem] = Triangulation(vem,Nmax);
%% Boundary Condition Domain
%% BCs:
addpath(genpath('./Domain'));
if n == 0
    [vem] = ByPassBoundCond(vem,fluid,BdBox);
    disp('Volume fraction (VolFrac):');volfrac = 0.33; disp(volfrac); 
    disp('penalty factor (q): '); q = 0.1; disp(q);
    
elseif n==1
    
    [vem]= ArterialBifurcationBoundCond(vem,fluid);
    disp('Volume fraction (VolFrac):'); volfrac = 0.24; disp(volfrac); 
    disp('penalty factor (q): '); q = 0.4 ;disp(q);  % q =0.4
    
elseif n==2
    [vem]= GraftArteryBoundCond(vem,fluid);
     disp('Volume fraction (VolFrac):'); volfrac = 0.10; disp(volfrac); 
     disp('penalty factor (q): '); q = 0.1 ;disp(q);
end
rmpath('./Domain');
%%
[vem] = DofDriveJ(vem);
%%
vem.ElemNotOpt = setdiff(1:vem.NElem,vem.ElemOpt)';
% verifying elements not optimize...    
hold on;      
for i = 1:length(vem.ElemNotOpt)
  vert = Node(Element{vem.ElemNotOpt(i)},:);
  plot(vert(:,1),vert(:,2),'r*');
end

%% ------------------------------------------ TOPOLOGY  OPTIMIZATION 
itermax0 = 50  ;   % number of iterations
nq = length(vem.ElemOpt); 
vt = 0; v = zeros(nq,1);
for i=1:nq
  v(i,:) = vem.ShapeFnc{vem.ElemOpt(i)}.W;
  vt = vt + v(i,1);
end
kmax = 2500; 
kmin = 0.00; 
vem.xp =   0.99*ones(nq,1);
vem.kappa =zeros(vem.NElem,1);
vem.dkappa =zeros(vem.NElem,1);
xold = vem.xp;
%% DEFINE VARIABLES AND PARAMETERS FOR MMA OPTIMIZATION ALGORITHM 
a0 = 1;
a = zeros(2,1);
c = 100*ones(2,1); 
d = zeros(2,1);
xmin = zeros(nq,1);
xmax = ones(nq,1);
xolder = xold;
low = zeros(nq,1);
upp = ones(nq,1);
% Two constraints
d2fdx2 = zeros(nq,1);
d2gdx2 = zeros(2,nq);
g = zeros(2,1);
dgdx= zeros(2,nq);
%% DESIGN LOOP FOR THE ACTUAL TOPOLOGY OPTIMIZATION
data = zeros(itermax0,2);
tic;
figure;
title('Optimal Structure');
for iter = 1:itermax0 
    
  for i =1:nq
     vem.kappa(vem.ElemOpt(i),1) = kmax +(kmin - kmax)*vem.xp(i,1)*(1+q)/(vem.xp(i,1) +q);
     vem.dkappa(vem.ElemOpt(i),1) = (kmin - kmax)*q*(1+q)/((vem.xp(i,1) +q)^2);
  end 
% ****  solve Navier Stokes Darcy problem  ****
disp('Newton Rapshon Method: ');
disp('Computing ...')
[C,vem] = VEM_NStokes(vem, fluid);
% ****  energy sensitivity analysis *****
[cE,dE] = SensitivityElemOpt(vem,fluid,C);
% ***** process of optimization  ****
f =    cE  ; data(iter,1) = f;
dfdx =  dE   ;
volf_iter = vem.xp'*v/vt ; data(iter,2) = volf_iter;
%  constraint function: volume
g(1) = volf_iter -volfrac;
dgdx(1,:) = v'/vt;
[xnew,yp,z,lambda,ksi,eta,mu,zeta,s,low,upp] = mmasub(2,length(vem.xp),iter, ...
vem.xp,xmin,xmax,xold,xolder,f,dfdx,d2fdx2,g,dgdx,d2gdx2,low,upp,a0,a,c,d);
xolder = xold; 
xold = vem.xp; 
gamma = xnew;
vem.xp = gamma; 
change = max(max(abs(vem.xp-xold)));
disp([' It.: ' sprintf('%4i',iter) ' Obj.: ' sprintf('%10.4f',f) ...
      ' Vol.: ' sprintf('%10.4f', volf_iter ) ...
      ' ch.: ' sprintf('%10.6f',change )])  
% PLOT DENSITIES  
vem.xvarp = ones(vem.NElem,1); vem.xvarp(vem.ElemOpt)=vem.xp;
[FigHandle,FigData] = InitialPlot(vem,vem.xvarp);
set(FigHandle,'FaceColor','flat','CData',vem.xvarp(FigData)); drawnow;
% TEST CONVERGENCE
     if iter >= itermax0 || change < 0.0001
     break
     end
end
temp_opt = toc;
figure;
title('Velocity Field');
PlotVelIntoTop(vem,C)
figure;
title('Pressure Field');
PlotPressIntoTop(vem,C)
figure;
PlotGamma2IntoTop(vem)
%% -----------------------------PLOT OBJECTIVE & CONSTRAINT FUNCTION 
figure(6);
set(6,'Name','Objective and Constraint function','Position',[766   476   560   420]);
[ax,h1,h2]=plotyy([1:1:length(data)],data(:,1)',[1:1:length(data)],data(:,2)');
set(get(ax(1),'YLabel'),'String','Objective Function','FontSize',12);
set(get(ax(2),'YLabel'),'String','Constraint Function','FontSize',12);
set(h1,'LineWidth',1)
set(h2,'LineWidth',1,'LineStyle','--')
xlabel('Iterations','FontSize',12);
 
 