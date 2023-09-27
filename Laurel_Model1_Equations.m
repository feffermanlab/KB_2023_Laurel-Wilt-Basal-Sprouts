function dx = Laurel_Model1_Equations(t,x,par)
dx=[0;0;0;0;0;0;0;0;0;0];
r_m=par(1);
r_l=par(2);
gamma=par(3);
g_p=par(4);
g_m=par(5);
g_b=par(6);
mu_d=par(7);
mu_p=par(8);
mu_m=par(9);
mu_l=par(10);
mu_b=par(11);
beta_p=par(12);
beta_m=par(13);
beta_l=par(14);
rho_m=par(15);
rho_l=par(16);
sigma=par(17);
v_p=par(18);
v_m=par(19);
v_l=par(20);
v_b=par(21);
v_0=par(22);
delta_p=par(23);
delta_m=par(24);
delta_l=par(25);
delta_b=par(26);
omega=par(27);


dx(1)=r_m*x(3)+r_l*x(4)-(gamma+mu_d)*x(1) ; %dD/dt SEEDS
dx(2)=gamma*x(1)-(g_p+mu_p)*x(2)-beta_p*x(10)*x(2); %dP_s/dt susceptible saplings
dx(3)=g_p*x(2)+g_b*x(5)-(g_m+mu_m)*x(3)-beta_m*x(10)*x(3); %dM_s/dt susceptible medium trees
dx(4)=g_m*x(3)-mu_l*x(4)-beta_l*x(10)*x(4);  % dL_s/dt susceptible large trees
dx(5)=rho_m/v_m*x(7)+rho_l/v_l*x(8)-(sigma+g_b+mu_b)*x(5); %dB_s/dt susceptible basal sprouts
dx(6)=v_p*beta_p*x(10)*x(2)-delta_p*x(6)*x(10); %dP_i/dt vol infected saplings
dx(7)=v_m*beta_m*x(10)*x(3)-delta_m*x(7)*x(10); % dM_i/dt vol infected medium trees
dx(8)=v_l*beta_l*x(10)*x(4)-delta_l*x(8)*x(10); %dM_i/dt vol infected large trees
dx(9)=v_b*sigma*x(5)-delta_b*x(9)*x(10); %dB_i/dt vol infected basal sprouts
dx(10)=omega*x(10)*(1-v_0*x(10)/(.0000000001+x(6)+x(7)+x(8)+x(9))); %dA/dt Adult female beetle population 

end

