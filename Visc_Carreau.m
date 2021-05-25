function [neta,d_neta,vem] =visc_carreau(dudx,dudy,dvdx,dvdy,vem) 
% ****   non-Newtonian model (Carreau–Yasuda)   *****
% ****   Parameters: 

visc_0 =  0.56;    % poise
visc_oo = 0.0345;  % poise

a = 2.0;
lamb = 3.313;
n =  0.3568;
ex = (n - 1.0)/a;


gamma =  2*dudx*dudx + 2*dvdy*dvdy +(dudy + dvdx)*(dudy + dvdx);      
gamma2 = sqrt(gamma)  ; vem.gamma2 = gamma2;

neta = visc_oo + (visc_0 - visc_oo)*( 1 + (lamb*gamma2)^a)^(ex); 

d_neta = (n-1)*lamb*(visc_0 - visc_oo)*(( 1 + (lamb*gamma2)^a)^(ex-1))*(lamb*gamma2)^(a-1);

end
