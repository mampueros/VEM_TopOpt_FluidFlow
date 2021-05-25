function [data] = Triangulation (data,Nmax)
node_num =  data.NNode; 
poly_num  = data.NElem;
pvol = zeros(node_num,1);
pneighb = zeros(node_num,4);
neighb = zeros(poly_num,Nmax);
neighb(:,:) = -1;
data.Neighb = data.Element;
% **** Neighboring volumes around a node  ****
for nel = 1: poly_num
     nl = size(data.Element{nel},2) ;
      for j = 1:nl
           ip = data.Element{nel}(j);
           if ip > 0
              pvol(ip) = pvol(ip) + 1;
              pneighb(ip, pvol(ip)) = nel;
           end 
      end
      if nl == 3
         neighb(nel,nl+1) = 0;   
      end
end

% **** Neighboring volumes around a volume   ****
  k = 0;
 for nel = 1: poly_num
        nl = size(data.Element{nel},2) ;
        data.Neighb{nel}(1:nl) = -1;
        for j = 1:nl            
            a1 = data.Element{nel}(j);
              jj = j+1;
              if jj > nl
                 jj = 1;
              end
              a2 = data.Element{nel}(jj);
              na1 = pvol(a1);
              na2 = pvol(a2);
              for k1 = 1:na1
                  viz1 = pneighb(a1,k1);
                  if viz1 ~= nel
                      for k2 = 1: na2
                         viz2 = pneighb(a2,k2);
                         if viz1 == viz2
                         %  neighb(nel,j) = viz1;    
                            data.Neighb{nel}(j) = viz1;
                         end
                      end
                  end
              end
              if data.Neighb{nel}(j) == -1, k = k+1; end
        end   
 end
 data.NElemBC = k ;
 %  **** Boundary elements  ****
  data.ElementBCJ = zeros(data.NElemBC,2);
  k = 0;
 for nel = 1: poly_num
        nl = size(data.Element{nel},2) ;
        for j = 1:nl
            if data.Neighb{nel}(j) == -1
                k = k + 1;
                data.ElementBCJ(k,1) = nel;
                data.ElementBCJ(k,2) = j;
            end
        end
 end
                

            