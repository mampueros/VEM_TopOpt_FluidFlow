function [] = PlotGamma2IntoTop(vem)
%% Pressure
hold on;
for iele = 1:vem.NElem
   NNodeEle = size(vem.Element{iele},2);
   Ze = zeros(NNodeEle,1);  Ge=Ze;
   iq = vem.Element{iele}(:) ;
   Xe  = vem.Node(iq,1);
   Ye  = vem.Node(iq,2); 
   ge = vem.Gamma2(iele,1);
   Ge(:,1) = ge;
   if vem.xvarp(iele) > 0.05
       fill3(Xe,Ye,Ze,Ge,'EdgeColor','none');
   else
       fill3(Xe,Ye,Ge,'k'); 
   end
end

colormap(jet); 
axis equal; axis tight; axis off;pause(1e-9);
colorbar('EastOutside');
hold off;

end

 