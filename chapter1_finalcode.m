%% Final Chapter 1 Analysis and Plots
%This is a copy paste from other files to get the most consise, clear code
%for all the analysis for chapter 1.

%This includes a short term sweep (and analsis of short term effects)
%long term effects in two spots of the parameter regieme 
%and a sweep of the whole parameter space to do stability analysis

% this requires stats and machine learning toolbox


%% Short term sweep

%% Define Parameters
%I want to see what affets rho, sig, and sig/rho have on the short term
%effects of the run
tic;
numberofcombos_short=251; %251 bc also defualt parameters.
 %%%%%%%% %%%%%%%% %%%%%%%% %%%%%%%% %%%%%%%% %%%%%%%%
% Set up combos using base parameters
% %We will vary r_l, gamma, g_p, mu_D and beta_L
params_refigured_simplified ;
 basicparams=par';
 clear par 

 %Define low and high range for each parameter 
 %column 3 is dummy variable
varyparam_info=zeros(27,3);
% vary r, gamma, g,mu and beta  50%-200% of default)
%r_data %vary from 50%-200% of basic par value
    varyparam_info(2,:)=[basicparams(2)/2,basicparams(2)*2,1];
    varyparam_info(3,:)=[basicparams(3)/2,basicparams(3)*2,1];
    varyparam_info(4,:)=[basicparams(4)/2,basicparams(4)*2,1];
    varyparam_info(7,:)=[basicparams(7)/2,basicparams(7)*2,1];
    varyparam_info(14,:)=[basicparams(14)/2,basicparams(14)*2,1];

%now do the varying with latin hyper cube sampling   
	variedindex=find(varyparam_info(:,3)>0);
    numvariedpars= length(variedindex);
	varyparam_info_condensed=varyparam_info( any(varyparam_info(:,3),2), : );
%do hypercube and transform using min max
	zeroonevalues=lhsdesign(numberofcombos_short-1,numvariedpars)';
    transformedlsh_short=varyparam_info_condensed(:,2).*zeroonevalues+  varyparam_info_condensed(:,1).*(1-zeroonevalues);
    clear  daysperyear
%related parameter co vary 
	parcombo_short=zeros(27,numberofcombos_short-1); %initalize matrix to hold all
	parcombo_short(variedindex,:)=transformedlsh_short; %save varied things in right spots
    %betas are related
    parcombo_short(13,:)=.25*parcombo_short(14,:);
    parcombo_short(12,:)=.01*parcombo_short(14,:); 
    %r_m=.9*r_l (based on orignal basic params
    parcombo_short(1,:)=.9*parcombo_short(2,:);
    %g_m is .4 g_m based on orig basic params, g_b=g_p
        parcombo_short(5,:)=.4*parcombo_short(4,:);
%add all other parameters which dont vary
    parcombo_short(~any(parcombo_short,2),:)=repmat(basicparams(~any(parcombo_short,2),1),1,numberofcombos_short-1);
    parcombo_short(6,:)=parcombo_short(4,:);%%
    
%Put each combo in its own spot in s structure (why? idk, but i need it)
 for eachcombo=1:numberofcombos_short-1
combodata(eachcombo).combo=parcombo_short(:,eachcombo);
 end
 %%%%%%%% %%%%%%%% %%%%%%%% %%%%%%%%
 %Set up the values of rho and sigma that I'll use for each par combo: 
 % for wihout BS and then *4 with BS

allrho=[0,.75,1.5,3]/365;
allsig=[0,0.5,1,2]/365;
allrhosigcombo=[];
for rho=1:length(allrho)
    if rho>1
        for sig=1:length(allsig)
        allrhosigcombo=[allrhosigcombo;allrho(rho), allrho(rho), allsig(sig)];
        end
    else
          allrhosigcombo=[allrhosigcombo;allrho(rho), allrho(rho), 0];

    end
end
    numberofpoints=length(allrhosigcombo);
    allrhosigcombo=allrhosigcombo';
 
%determine how many rho and sigs to look at
%    numberofpoints=20; %55; made smaller when just analyizing trends, not making a graph
   outputs_short(1).description="varyrhosig";


 
%% Run ODE over short term
%this runs the ODE over a few years for all parameter combos, for all rho
%and all sigma.
    
%save combos in outputs short
for i=1:4
    outputs_short(i).combospecific=combodata;
end

clear parcombo
 %%%%%%%% %%%%%%%% %%%%%%%% %%%%%%%% %%%%%%%% %%%%%%%%
%Run and save ode45 for all the runs
%define itits
      syms t x y z
    years_short=7; tmax=365*years_short; tspan_short=[0 tmax];
    init.D_0=140;   init.P_s0=80;   init.M_s0=60;   init.L_s0=20;   init.B_s0=0;
    init.P_i0=0;    init.M_i0=0;    init.L_i0=1;    init.B_i0=0;    init.A_0=1;
    inits2=[init.D_0; init.P_s0; init.M_s0; init.L_s0; init.B_s0; init.P_i0; init.M_i0; init.L_i0; init.B_i0; init.A_0];
    options = odeset('MaxStep',365*1e-2);

     type=1;

        for combonum=1:numberofcombos_short
            if combonum==numberofcombos_short
             combotemp=basicparams;

            else
            combotemp=outputs_short(type).combospecific(combonum).combo;
            end
         for range=1:numberofpoints %range is the range of values of rho and sigma 
        %define short term par combo
        combotemptemp=combotemp;
        %define rho and sigma
    
           combotemptemp(15:17,1)=allrhosigcombo(:,range);

        %run ODE
            [T,Y]=ode45(@Laurel_Model1_Equations, tspan_short, inits2, options,combotemptemp);
            outputs_short(type).combospecific(combonum).allTs(range).Ts=T;
            outputs_short(type).combospecific(combonum).allYss(range).Ys=Y;
        %save
        end
        end
     

clear z Y T  combodata combonum comnotemp  eachcombo outputs type transformed lsh t range x y tmax i
toc;
%% All Short term Analysis 

%analysis saved in outputs short

%for each population, save max and time to max
tic;
type=1;
for combonum=1:numberofcombos_short
for range=1:numberofpoints
for state=10
[a,b]=max(outputs_short(type).combospecific(combonum).allYss(range).Ys(:,state));
outputs_short(type).combospecific(combonum).vecpeakVAL_all(range)=a;
outputs_short(type).combospecific(combonum).vecpeakTIME_all(range)=outputs_short(type).combospecific(combonum).allTs(range).Ts(b,1)/365;



if range==1
    outputs_short(type).combospecific(combonum).vecpeakVAL_noBS(range)=a;

outputs_short(type).combospecific(combonum).vecpeakTIME_noBS(range)=outputs_short(type).combospecific(combonum).allTs(range).Ts(b,1)/365;
end

    outputs_short(type).combospecific(combonum).vecpeakTIME_diff=max(abs(outputs_short(type).combospecific(combonum).vecpeakTIME_all-outputs_short(type).combospecific(combonum).vecpeakTIME_all(1)));
    outputs_short(type).combospecific(combonum).vecpeakVAL_diff=max(abs(outputs_short(type).combospecific(combonum).vecpeakVAL_all-outputs_short(type).combospecific(combonum).vecpeakVAL_all(1)));

end
end
end


%time to 80% mortality 
type=1;
bad=[];
for combonum=1:numberofcombos_short
if min(outputs_short(type).combospecific(combonum).allYss(1).Ys(:,4))<4
    for range=1:numberofpoints

        if range==1
        outputs_short(type).combospecific(combonum).Ls90time_noBS(1)=outputs_short(type).combospecific(combonum).allTs(range).Ts(find(outputs_short(type).combospecific(combonum).allYss(range).Ys(:,4)<4,1),1)/365;
        end
        
       outputs_short(type).combospecific(combonum).Ls90time_all(1,range)=outputs_short(type).combospecific(combonum).allTs(range).Ts(find(outputs_short(type).combospecific(combonum).allYss(range).Ys(:,4)<    4.0008,1),1)/365;

    end
    
    outputs_short(type).combospecific(combonum).Ls90time_diff=max(abs(outputs_short(type).combospecific(combonum).Ls90time_all-outputs_short(type).combospecific(combonum).Ls90time_all(1)));

    
    
else
    bad=[bad,combonum];
end
end
toc;
display('NO BASAL sprouts: default parameters, time to 80 ')
outputs_short(type).combospecific(numberofcombos_short).Ls90time_noBS(1,1)
display('NO BASAL sprouts: default parameters, time of beetle peak ')
outputs_short(type).combospecific(numberofcombos_short).vecpeakTIME_noBS(1)

display('NO BASAL sprouts: varied parameters, time to 80 range')
[min([outputs_short(type).combospecific.Ls90time_noBS]),max([outputs_short(type).combospecific.Ls90time_noBS])]
display('NO BASAL sprouts: varied parameters, time of beetle peak ')
[min([outputs_short(type).combospecific.vecpeakTIME_noBS]),max([outputs_short(type).combospecific.vecpeakTIME_noBS])]



display('WITH BASAL sprouts: varied parameters, change to time to 80 range')
max([outputs_short(type).combospecific.Ls90time_diff])



%Proof that he big differences occur when beta is small:
%     find([outputs_short(type).combospecific.Ls90time_diff]>0.5*max([outputs_short(type).combospecific.Ls90time_diff]));
%     max(zeroonevalues(5,ans))
% 
% 
%     find([outputs_short(type).combospecific.Ls90time_diff]>0.05);
%     max(zeroonevalues(5,ans))


display('WITH BASAL sprouts: varied parameters, change to time to beetle peak')

max([outputs_short(type).combospecific.vecpeakTIME_diff])


display('WITH BASAL sprouts: varied parameters, change to number to beetle peak')

max([outputs_short(type).combospecific.vecpeakVAL_diff])

display('these are the values discussed in Results paragraph 1')

% %Proof that the size of the vector peak increses with BOTH increasing rho
% %and increasing sigma:
%     %pull out just the values of the peaks for all combos, all 13 rhosig
%     test=reshape([outputs_short(type).combospecific.vecpeakVAL_all],13,251)';
%     %increasing sigma
%     sum(diff(test(:,[2:5]),1,2)>0);
%     sum(diff(test(:,[6:9]),1,2)>0);
%     sum(diff(test(:,[10:13]),1,2)>0);
%     %increasing rho
%     sum(diff(test(:,[1,2,6,10]),1,2)>0);
%     sum(diff(test(:,[1,3,7,11]),1,2)>0);
%     sum(diff(test(:,[1,4,8,12]),1,2)>0);
%     sum(diff(test(:,[1,5,9,13]),1,2)>0);
%     %all are 251,thus increasing with increasing rho and sigma
% 














%%
%THIS IS A PAPER PLOT (Figure 2)
%this uses the DEFAULT parameter values (except when rho or sig are zero)
% I eant to make one plot that shows that rho and sig do not affect the
% short term dynamics.
close all
%run the model for 3 years for default parameters, with rho=0 sig=0,
%rho=pos,sig=0, rho=pos sig=pos
     years_short=3; tmax=365*years_short; tspan_short=[0 tmax]
         inits2=[init.D_0; init.P_s0; init.M_s0; init.L_s0; init.B_s0; init.P_i0; init.M_i0; init.L_i0; init.B_i0; init.A_0];

rhosigonoff=[ 0,0,0;1,1,0;1,1,1]'.*[1.5,1.5,1]'/365'
    plotstyle=["-k", "--k", ":k"];
    basicparams_rhosig=basicparams;
for i=1:3
    basicparams_rhosig(15:17)=rhosigonoff(:,i);

[T,Y]=ode45(@Laurel_Model1_Equations, tspan_short, inits2, options,basicparams_rhosig);
f=figure(5)
subplot(1,2,1)
hold on
plot(T/365, Y(:,4), plotstyle(i), 'LineWidth',2)
subplot(1,2,2)
hold on
plot(T/365, Y(:,10), plotstyle(i), 'LineWidth',2)

end

 tree=["Susceptible  Sapling", "Susceptible Medium Tree", "Susceptible Large Tree", "Beetle"];
s=3
    subplot(1,2,s-2)
    legend("\rho = 0, \sigma = 0", "\rho > 0, \sigma = 0", "\rho > 0, \sigma > 0",'FontSize', 16)

 xlabel({'Time (Years)', '(a)'}, 'FontSize', 18)
 xlim([0,3])
 ylabel(strcat(tree(s)," Population"),'FontSize', 18)
 if s==4
     ylim([0 200])
 end
s=4
     subplot(1,2,s-2)
    legend("\rho = 0, \sigma = 0", "\rho > 0, \sigma = 0", "\rho > 0, \sigma > 0",'FontSize', 16)

 xlabel({'Time (Years)', '(b)'}, 'FontSize', 18)
 xlim([0,3])
 ylabel(strcat(tree(s)," Population"),'FontSize', 18)
 if s==4
     ylim([0 200])
 end

%% other ch 1 graph: Ms* vs sig/rho and vs rho.(long term)
%Figures 3 and 4 in paper
close all
linestyle_plot1=["--k", "-k"];
linestyle_plot2=["-k", "--k", "-.k", "k:"];


for run=1:2
part2pars=basicparams;      
%update gamma r_Z, g_Z and mu_D
if mod(run,2)==0
part2pars([1,2,3])=2*part2pars([1,2,3]);
end
if mod(run,2)==1
part2pars([1,2,3])=.5*part2pars([1,2,3]);
end


% if run>2
%     if mod(run,2)==0
% part2pars(4:6)=2*part2pars(4:6);end
% if mod(run,2)==1
% part2pars([4,5,6])=.5*part2pars([4,5,6]);
% end
%  end      
%     
tic;
numpoints=60;%make smaller for test runs
medeq=zeros(1,numpoints);
%allrho=logspace(-5,-2,numpoints); %rho from zero to about 3 per year %or 1.75 for about 6 years?
allrho=linspace(0,3/365,numpoints); %rho from zero to about 3 per year %or 1.75 for about 6 years?



part2pars(17)=0;
%run the ode and save the end value
      syms t x y z
years=8000; tmax=365*years; tspan=[0 tmax];
init.D_0=140;   init.P_s0=80;   init.M_s0=60;   init.L_s0=20;   init.B_s0=0;
init.P_i0=0;    init.M_i0=0;    init.L_i0=1;    init.B_i0=0;    init.A_0=1;
inits2=[init.D_0; init.P_s0; init.M_s0; init.L_s0; init.B_s0; init.P_i0; init.M_i0; init.L_i0; init.B_i0; init.A_0];
options=[];
medeq_rho=zeros(1,numpoints);

for i=1:numpoints
    part2pars(15:16)=[allrho(i), allrho(i)];
    [T,Y]=ode45(@Laurel_Model1_Equations, tspan, inits2, options,part2pars);
    
medeq_rho(i)=Y(length(Y),3)';
%figure(100+i)
%plot(T,Y)
end

long1=figure(100)
hold on
plot(allrho*365, medeq_rho, linestyle_plot1(run),'LineWidth',2)
toc;


        tic;

%now do the rho sig thing
slope=[0,logspace(-5,1,numpoints-1)];
%part2pars([1,2,3,4,5,6])=2*part2pars([1,2,3,4,5,6]);
%part2pars(7)=.5*part2pars(7);
%part2pars(12:14)=.25*part2pars(12:14);
allrho=[0.01,.75, 1.5, 3]/365;
medeq_rhosig=zeros(length(allrho),numpoints);

    for r=2:length(allrho)
   % for r=1:length(allrho) %use this one to also display near zero rho
    part2pars(15:16)=[allrho(r), allrho(r)];      

    for s=1:numpoints
    part2pars(17)=allrho(r)*slope(s);    
    
    
        [T,Y]=ode45(@Laurel_Model1_Equations, tspan, inits2, options,part2pars);

    medeq_rhosig(r,s)=Y(length(Y),3)';

    end
  long2=  figure(101)
subplot(1,2,run)

semilogx(slope, medeq_rhosig(r,:), linestyle_plot2(r),'LineWidth',2)
  hold on  
    end


    toc;



end


figure(100)

% title('Effect of Basal Sprout Production with no Secondary Infection','FontSize', 14)

 xlabel('Basal Sprout Production Rate (\rho)', 'FontSize', 16)
 ylabel('Medium Tree Population at Equilibrium','FontSize', 16)
% legend("High Seed Production and Germination","Poor Host Growth")
 
 figure(101)
 for run=1:2
 subplot(1,2,run)
  xlabel({'Basal Sprout Infection to Production Ratio (\sigma/\rho)',['(',char('a'+(run-1)),')']}, 'FontSize', 16)
 ylabel('Medium Tree Population at Equilibrium','FontSize', 16)
 xlim([0,10])

 %legend('Near Zero \rho','Small \rho', 'Moderate \rho', 'Large \rho','FontSize', 12)
%use above legend when changing line 484
 legend('Small \rho', 'Moderate \rho', 'Large \rho','FontSize', 12)
%  title('Effect of Basal Sprout Infection Rate','FontSize', 14)
 end
 



 %export
exportgraphics(long1,'long1.png','Resolution',300)
exportgraphics(long2,'long2.png','Resolution',300)


 %% Endemic eq value at default parameters (in appendix)
 
       syms t x y z
years=8000; tmax=365*years; tspan=[0 tmax];
init.D_0=140;   init.P_s0=80;   init.M_s0=60;   init.L_s0=20;   init.B_s0=0;
init.P_i0=0;    init.M_i0=0;    init.L_i0=1;    init.B_i0=0;    init.A_0=1;
inits2=[init.D_0; init.P_s0; init.M_s0; init.L_s0; init.B_s0; init.P_i0; init.M_i0; init.L_i0; init.B_i0; init.A_0];
options=[];



    [T,Y]=ode45(@Laurel_Model1_Equations, tspan, inits2, options,basicparams);
    
Y(length(Y), :)
 
%% stability of system:
%% 4: Initiate parameter sweeps
%number of parameter combos to consider
    vary.numberofcombos=20;

    
%ode.initate base params and ranges
params_refigured_simplified ;   basicparams=par';
 clear par
vary.param_info=zeros(27,3);
%1:2- pars.r_l. (pars.r_m related)
    vary.param_info(2,:)=[basicparams(2)/2,basicparams(2)*2,0];%vary from 50%-200% of basic par value. 
%3- pars.gamma_data=
    varyparam_info(3,:)=[basicparams(3)/2,basicparams(3)*2,1];
%4:6- pars.g_p (note: gm and gb depdend on pars.g_p later on)
    vary.param_info(4,:)=[basicparams(4)/2,basicparams(4)*2,0];%vary from 50%-200% of basic par value
%7- pars.mu_d=
    varyparam_info(7,:)=[basicparams(7)/2,basicparams(7)*2,1];
%8:10 pars.mu_l (other mup and mum related below)
    vary.param_info(10,:)=[basicparams(10)/2, basicparams(10)*2,0];
%11 pars.mu_b 
    vary.param_info(11,:)=[basicparams(11)/2, basicparams(11)*2,0];
%12:14- betal (note: betal and beta m depend on beta l later on)  
    vary.param_info(14,:)=[basicparams(14)/2,basicparams(14)*2,0];
% 15:16- pars.rho_m (note pars.rho_l equal
    vary.param_info(15,:)=[0.00001,0.035,0];
%17-pars.sigma
    vary.param_info(17,:)=[0.0001,0.03,0];
%18:22 related to volumes of trees. Do not vary.
%23:26-pars.delta_p (other deltas equal)
    vary.param_info(23,:)=[basicparams(23)/2, basicparams(23)*2,0];
%23-pars.omega
    vary.param_info(27,:)=[basicparams(27)/2, basicparams(27)*2,0];
%vary.type is user defined in decisions. based on this, we vary
       %2r, 3pars.gamma, 4g, 7mud, 10mut, 11mub, 14beta, 15rho, 17pars.sigma, 23detlta, 27
       %pars.omega
    vary.param_info([2,3,4,7,10,11,14,15,17,23,27],3)=1;

vary.variedindex=find(vary.param_info(:,3)>0);
vary.numvariedpars= length(vary.variedindex);
  vary.param_info_condensed=vary.param_info( any(vary.param_info(:,3),2), : );
  %do hypercube
  vary.zeroonevalues=lhsdesign(vary.numberofcombos,vary.numvariedpars)';
  
  %transform using min max
  
 transformedlsh=vary.param_info_condensed(:,2).*vary.zeroonevalues+  vary.param_info_condensed(:,1).*(1-vary.zeroonevalues);

 clear  daysperyear
  %add values diretly related to others so they covary
  parcombo=zeros(27,vary.numberofcombos);
  parcombo(vary.variedindex,:)=transformedlsh;

%alter related quantities for small
%rs are related. pars.r_m=.9*pars.r_l (based on orignal basic params
    parcombo(1,:)=.9*parcombo(2,:);   
%gs are related   %gm=.4 pars.g_p
    parcombo(5,:)=.4*parcombo(4,:);         parcombo(6,:)=parcombo(4,:);
%betas are related
    parcombo(13,:)=.25*parcombo(14,:);      parcombo(12,:)=.01*parcombo(14,:); 
%rhos are realted 
    parcombo(16,:)=parcombo(15,:); 
     vary.type="some" %added this bc of error below
if vary.type=="all"
 %mus are related %pars.mu_m=2pars.mu_l       $pars.mu_p=3 pars.mu_l
    parcombo(9,:)=2*parcombo(10,:);         parcombo(8,:)=3*parcombo(10,:);
%deltas are related 
    parcombo(24,:)=parcombo(23,:);          parcombo(25,:)=parcombo(23,:);
    parcombo(26,:)=parcombo(23,:); 
end
    %add basic params to any parameter not varied
parcombo(~any(parcombo,2),:)=repmat(basicparams(~any(parcombo,2),1),1,vary.numberofcombos);

clear basicparams transformedlsh

%%

 
 %%
%functions
function value=evaluateatvariedparams(expression,parvars, parcombo)

 %parvars=[r_m; r_l; gamma; g_p; g_m; g_b; mu_d; mu_p; mu_m; mu_l; mu_b; beta_p; beta_m; beta_l; rho_m; rho_l; sigma; v_p; v_m; v_l; v_b; v_0; delta_p; delta_m; delta_l; delta_b; omega];

for k=1:length(parvars)
    expression=subs(expression, parvars(k), parcombo(k));
end
value=expression;
end

function tempoutput=determinestability(tempeigs)
realpart=real(tempeigs);
             imagpart=imag(tempeigs);
             negrealpart=realpart(realpart<0);
if min(abs(realpart))<1e-5
    tempoutput="zeroeigs";
else
             if length(negrealpart)==length(realpart) %if real parts are all neg
                if isempty(imagpart(imagpart~=0))==0 %if complex
                    tempoutput="stab_cycle"; %stab cycle
                else
                    tempoutput="stab_node"; %   stab node';
                end
             else
                 %  some real parts pos
                     if isempty(imagpart(imagpart~=0))==0
                         tempoutput="unstab_cycle"; %'unstab cycle';
                     else
                        tempoutput="unstab_node";%'unstab node';
                    end
             end
end
end
