function [] = PlotPressIntoTop(vem,C)
%% Pressure
hold on;
for iele = 1:vem.NElem
   NNodeEle = size(vem.Element{iele},2);
   Ze = zeros(NNodeEle,1);   Pe=Ze;
   iq = vem.Element{iele}(:) ;
   Xe  = vem.Node(iq,1);
   Ye  = vem.Node(iq,2); 
   p = C(vem.ASSMtrx(2*NNodeEle+1,iele));
   Pe(:,1) = p;
   if vem.xvarp(iele) > 0.05
       fill3(Xe,Ye,Ze,Pe,'EdgeColor','none');
   else
       fill3(Xe,Ye,Pe,'k'); 
   end
end

colormap(jet); 
axis equal; axis tight; axis off;pause(1e-9);
colorbar('EastOutside');
hold off;

end

 