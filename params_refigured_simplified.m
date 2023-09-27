%%params refigured 
%%simplified

%the point of this script it to round the numbers used to make it easier to
%report


par=zeros(1,27);
%parameter script
%this is updated after getting better estimates in a large parameter sweep
daysperyear=365; %number of totaldays in a "year"--might change this to ignore winter
annualsurvivorship_sapling=.97; %annual percent survival for sapling sassafras
annualsurvivorship_medium=.98; %annual percent survival for medium sassafras
annualsurvivorship_large=.99; %annual percent survival for large sassafras

 % par=[r_m; r_l; gama; g_p; g_m; g_b; mu_d; mu_p; mu_m; mu_l; mu_b; beta_p; beta_m; beta_l; rho_m; rho_l; sig; v_p; v_m; v_l; v_b; v_0; delta_p; delta_m; delta_l; delta_b; omega];


par(2)=1/3/daysperyear; %r_l
%was 0.31, changeed to .3333

par(1)=0.9*par(2); %r_m
%was 0.28, changed to 90% of par 2

par(3)=.1/daysperyear; %gamma %between 6 and 16 percent germination 
par(4)=1/(12.5*daysperyear); %g_p%(2.5-0)/.2
par(5)=1/(30*daysperyear); %g_m %(10-2.5)/.25
par(6)=1/(12.5*daysperyear);%g_b %assume same as sapplings
par(7)=1/(3*daysperyear); %mu_d% sass seeds last 6 years in seed bank. so life span is 6*totaldays. rate is 1 over that
%par(8)=1-(annualsurvivorship_sapling)^(1/daysperyear); %mu_p
%par(9)=1-(annualsurvivorship_large)^(1/daysperyear);%mu_m
%par(10)=1-(annualsurvivorship_large)^(1/daysperyear);%mu_l


%updated to havce 1% death each year, consistent with whats known
par(10)=0.01/daysperyear;
 %par(10)=3.46544914251165e-05; %updated 
 par(9)=2*par(10);
 par(8)=3*par(10);

 %%this feels very wrong below. can I change it??
par(11)=3.4247e-06 ;% approx .0013/365; %mu_b


%par(12)=1.00000000000000e-07; %beta_p
%par(13)=4.61428571428571e-05;%beta_m
par(14)=0.025/daysperyear
%par(14)=4.85714285714286e-05;%beta_l
par(12)=par(14)/100;
par(13)=par(14)/4;

%changed 1/3/22
par(15)=1.5/daysperyear;
par(16)=1.5/daysperyear;
par(17)=1/daysperyear;

%par(15)=0.0188888888888889;  %rho_m
%par(16)=0.0216666666666667;%rho_l
%par(17)=0.00215443469003188;%sigma

par(18)=1;%assume basal and sappling are same size %v_p
par(19)=2*par(18);%assume mediuem is twice as big as sappling %v_m
par(20)=4*par(18);%asume large is twice as big as medium. %v_l
par(21)=par(18); %assume basal and sappling are same size %v_p
par(22)=1/3; %vo-updated
%^^changed, was 0.3

%par(23)=4.5000e-05; %delta_p all updated
%par(24)=4.5000e-05; %delta_m
%par(25)=4.5000e-05; %delta_l
%par(26)=4.5000e-05; %delta_b

par(23)=1/60/daysperyear; %delta_p all updated
par(24)=par(23); %delta_m
par(25)=par(23); %delta_l
par(26)=par(23); %delta_b

%par(27)=0.0225; %omega updated
par(27)=8/daysperyear; %omega updated


clear annualsurvivorship_large
clear annualsurvivorship_medium
clear annualsurvivorship_sapling



