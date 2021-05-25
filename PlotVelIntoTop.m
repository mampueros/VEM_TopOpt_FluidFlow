function [] = PlotVelIntoTop(vem,C)
hold on;
for iele = 1:vem.NElem
   NNodeEle = size(vem.Element{iele},2);
   Ze = zeros(NNodeEle,1);   Ve=Ze;
   iq = vem.Element{iele}(:) ;
   vx = vem.ASSMtrx(1:NNodeEle,iele);
   vy = vem.ASSMtrx(NNodeEle+1:2*NNodeEle,iele);
   VXG  = C(vx);
   VYG  = C(vy); 
    Xe  = vem.Node(iq,1);
    Ye  = vem.Node(iq,2);  
   for ilnode = 1:NNodeEle
      Ve(ilnode) = sqrt(VXG(ilnode)^2+VYG(ilnode)^2);
   end
   if vem.xvarp(iele) > 0.05
      fill3(Xe,Ye,Ze,Ve,'EdgeColor','none'); 
   else
      fill3(Xe,Ye,Ve,'k'); 
   end
end
colormap(jet); 
axis equal; axis tight;   axis on; pause(1e-9);
colorbar('EastOutside');
hold off;
axis('off');
end