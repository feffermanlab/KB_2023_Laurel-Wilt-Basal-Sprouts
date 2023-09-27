%condensed stability analysis
%this copy and pastes the shortest and best code for equilibirum existence
%and stability analysis code into one place. It goes nicely with the
%latex document which summarizes it (supplemental material)
clear
clc

%% 0: decisions


%% 1: universal set up

clear diff
syms D Ps Ms Ls Bs Pi Mi Li Bi A  
symbolicnames.states=["D", "Ps", "Ms", "Ls", "Bs", "Pi", "Mi", "Li", "Bi", "A"];
statestemp=[D Ps Ms Ls Bs Pi Mi Li Bi A  ];
for s=1:10
    states.(symbolicnames.states(s))=statestemp(s);
end
clear D Ps Ms Ls Bs Pi Mi Li Bi A statestemp

syms  r_m r_l gamma g_p g_m g_b mu_d mu_p mu_m mu_l mu_b 
syms beta_p beta_m beta_l rho_m rho_l sigma v_p v_m v_l v_b v_0 delta_p delta_m delta_l delta_b omega

symbolicnames.pars=[ "r_m", "r_l", "gamma", "g_p", "g_m", "g_b", "mu_d", "mu_p", "mu_m", "mu_l", "mu_b",  "beta_p", "beta_m", "beta_l", "rho_m", "rho_l", "sigma", "v_p", "v_m", "v_l", "v_b", "v_0", "delta_p", "delta_m", "delta_l", "delta_b", "omega"];
parstemp=[r_m; r_l; gamma; g_p; g_m; g_b; mu_d; mu_p; mu_m; mu_l; mu_b; beta_p; beta_m; beta_l; rho_m; rho_l; sigma; v_p; v_m; v_l; v_b; v_0; delta_p; delta_m; delta_l; delta_b; omega];
for p=1:27
   pars.(symbolicnames.pars(p))=parstemp(p) ;
end
clear  r_m r_l gamma g_p g_m g_b mu_d mu_p mu_m mu_l mu_b 
clear beta_p beta_m beta_l rho_m rho_l sigma v_p v_m v_l v_b v_0 delta_p delta_m delta_l delta_b omega
clear parstemp
%for k=1:length(parvars)
 %   assume(parvars(k), 'positive');
%end

syms Mistar Pistar
stars.Mi=Mistar; stars.Pi=Pistar;
clear Mistar Pistar

eqns.dD=pars.r_m*states.Ms+pars.r_l*states.Ls-(pars.gamma+pars.mu_d)*states.D ; %dstates.D/dt SEEDS
eqns.dPs=pars.gamma*states.D-(pars.g_p+pars.mu_p)*states.Ps-pars.beta_p*states.A*states.Ps; %dP_s/dt susceptible saplings
eqns.dMs=pars.g_p*states.Ps+pars.g_b*states.Bs-(pars.g_m+pars.mu_m)*states.Ms-pars.beta_m*states.A*states.Ms; %dM_s/dt susceptible medium trees
eqns.dLs=pars.g_m*states.Ms-pars.mu_l*states.Ls-pars.beta_l*states.A*states.Ls;  % dL_s/dt susceptible large trees
eqns.dBs=pars.rho_m/pars.v_m*states.Mi+pars.rho_l/pars.v_l*states.Li-(pars.sigma+pars.g_b+pars.mu_b)*states.Bs; %dB_s/dt susceptible basal sprouts
eqns.dPi=pars.v_p*pars.beta_p*states.A*states.Ps-pars.delta_p*states.Pi*states.A; %dP_i/dt vol infected saplings
eqns.dMi=pars.v_m*pars.beta_m*states.A*states.Ms-pars.delta_m*states.Mi*states.A; % dM_i/dt vol infected medium trees
eqns.dLi=pars.v_l*pars.beta_l*states.A*states.Ls-pars.delta_l*states.Li*states.A; %dM_i/dt vol infected large trees
eqns.dBi=pars.v_b*pars.sigma*states.Bs-pars.delta_b*states.Bi*states.A; %dB_i/dt vol infected basal sprouts
eqns.dA=pars.omega*states.A*(1-pars.v_0*states.A/(states.Pi+states.Mi+states.Li+states.Bi)); %dA/dt Adult female beetle population 

symbolicnames.eqns=["dD", "dPs", "dMs", "dLs", "dBs", "dPi", "dMi", "dLi", "dBi", "dA"];

clear dD dPs dMs dLs dBs dPi dMi dLi dBi dA

syms J
for i=1:10
    for j=1:10
J(i,j)=diff(eqns.(symbolicnames.eqns(i)),states.(symbolicnames.states(j)));
    end
end
clear   s q p i j
%now, go either to section 2 or 3
origeqns=eqns;
%% 2: Vector free equilibria first. TLDR: two exist. both unstable.
%don't run 2 if you want to run 3
%start with A bc that gives eqn for A
disp('Set dA=0. Then')
stars.A=solve(eqns.dA==0,states.A)

%case 1: A=0
%infected classes
%Case 1: Astar=0. Sub A=0 into infected P M and L. ')
eqns.dPi=subs(eqns.dPi,states.A,stars.A(1)); eqns.dMi=subs(eqns.dMi,states.A,stars.A(1)); eqns.dLi=subs(eqns.dLi,states.A,stars.A(1));
%'That implies [dPi,dMi,dLi]    infected trees rates are zero. so we dont know Mi Li or Pi yet')
%Now consider dBi=0. that implies')
eqns.dBi=subs(eqns.dBi,states.A,stars.A(1)); stars.Bs=solve(eqns.dBi==0,states.Bs)
%Thus we have no sus basal sprouts. Now, we set dB=0 and plug in Bsstar=0.
%That implies that Li is a NEGATIVE mulitple of Mi. 
eqns.dBs=subs(eqns.dBs,states.Bs,stars.Bs);    eqns.dBs=subs(eqns.dBs,states.Mi,stars.Mi); stars.Li=solve(eqns.dBs==0, states.Li);
%^This is not bio-feasible. So it must be that Mi=Li=0'
%'Now sus classes: We plug Astar=0, Bsstrar=0 into all')
eqns.dPs=subs(eqns.dPs,states.A,stars.A(1));   eqns.dMs=subs(eqns.dMs,states.A,stars.A(1));
eqns.dMs=subs(eqns.dMs,states.Bs,stars.Bs);    eqns.dLs=subs(eqns.dLs,states.A,stars.A(1));
%We set dLs=0 and solve for Ls in terms of Ms and plug that into all cases of Ls')
stars.Ls=solve(eqns.dLs==0, states.Ls); eqns.dD=subs(eqns.dD, states.Ls, stars.Ls);
%now set dMs=0 and sovle for Ps in terms of Ms and plug that into instances of Ps')
stars.Ps=solve(eqns.dMs==0, states.Ps); 
eqns.dPs=subs(eqns.dPs, states.Ps, stars.Ps);
eqns.dMs=subs(eqns.dD, states.Ps, stars.Ps);
%now we set dD=0 and sove for D in terms of Ms')
stars.D=solve(eqns.dD==0,states.D);
%'now we plug that into dP ')
eqns.dPs=subs(eqns.dPs, states.D, stars.D);
vectorfreestuff.equalzero=factor(eqns.dPs);
%either Ms=0 which makes everything zeroor this parameter exression is zero. boundary case')
vectorfreestuff.gamma_boundary=solve(vectorfreestuff.equalzero(3)==0, pars.gamma); %factor numerator and demom to get nice expression for this
%note: This parameter equality is bio feasible. Thus we have a disease free eq
%Technically, Pi is open. We said Li=Mi=0. and Bi is open as well.
%summary: two vector free equilbirums. One is total extinction (with some
%leftover unitilized material) and one is disease free (under boundary
%condition for pars.gamma)

%gamma boundary needs to be positive for any growth to happen (because if
%gamma boundary is neg, the disease free eq doesnt exist, so extinction is
%only option for population without
[~,d]=numden(vectorfreestuff.gamma_boundary)
d=-d;

d=subs(d, pars.mu_p, 3*pars.mu_l);
d=subs(d, pars.mu_m, 2*pars.mu_l);
d=subs(d, pars.g_m, 0.4*pars.g_p);
d=subs(d, pars.r_m, 0.9*pars.r_l);

vectorfreestuff.rboundary=solve(d==0, pars.r_l);




%STABILITY of vector free. Start with jacobain
%%%MUst first do intro part
%sub A=0 into jacobiam
vectorfreestuff.Jvectorfree=J; syms x
vectorfreestuff.Jvectorfree=subs(vectorfreestuff.Jvectorfree,states.A, 0);
vectorfreestuff.charpolyvectorfree=charpoly(vectorfreestuff.Jvectorfree,x);
%note: this doesn't have any state variables in it after only subbing A
factor(vectorfreestuff.charpolyvectorfree);
%%%%note: the characterstic polynomial has a factor of x^4. but also it factors,
%%%%and (pars.omega-x) is a factor. since pars.omega is positive, we have a postive
%%%%root, and an unstable system. for BOTH equilibiurm.

vectorfreestuff.Jdiseasefree_eigs=eig(vectorfreestuff.Jvectorfree);
%this is useless^ 

vectorfreestuff.charpolydiseasefree=charpoly(vectorfreestuff.Jvectorfree,x)
factor(vectorfreestuff.charpolydiseasefree)
%%%%the characterstic polynomial has a factor of x^4. but also it factors,
%%%%and (pars.omega-x) is a factor. since pars.omega is positive, we have a postive
%%%%root, and an unstable system

display('this supports the factoring on page 15 (appendix)')
%% 3: Endemic equilibrium 

%Before running this, rerun the universial set up (section 1)

%this section gives the equilibrium expresions and the boundary equations 

%for k=1:length(parstemp)
%    assume(parstemp(k), 'positive');
%end

for i=1:10
   stars.(symbolicnames.states(i))=states.(symbolicnames.states(i)) ;
end


%start with A bc that gives eqn for Set dA=0. take A nonzero
Astar_all=solve(eqns.dA==0,states.A);   stars.A=Astar_all(2);
clear Astar_all
%We divide out A in spaces where we can factor out A 
%and we plug that expression for A into each equaiton with an A 
eqns.dPi=simplify(eqns.dPi/states.A);  
eqns.dMi=simplify(eqns.dMi/states.A);
eqns.dLi=simplify(eqns.dLi/states.A);  

for q=1:10
eqns.(symbolicnames.eqns(q))=subs(eqns.(symbolicnames.eqns(q)), states.A, stars.A);
end
%Then we solve for infected L M P in terms of sus using dinfecteds. sub
%into all stars and eqns
for q=6:8
stars.(symbolicnames.states(q))=solve(eqns.(symbolicnames.eqns(q))==0, states.(symbolicnames.states(q)));
for qq=1:10
 eqns.(symbolicnames.eqns(qq))=subs(eqns.(symbolicnames.eqns(qq)), states.(symbolicnames.states(q)), stars.(symbolicnames.states(q)));
stars.(symbolicnames.states(qq))=subs(stars.(symbolicnames.states(qq)), states.(symbolicnames.states(q)), stars.(symbolicnames.states(q)));
end
end



%we now have ceqna 1-5 and 9 in terms of vars 1-5 and 9
%lets keep solve dBi==0 for Bs because 1, its the shortest, and 2, because
%Bs only shows up in the Bs and M equations, but Ps Ms and LS show up a lot
%more

stars.Bs=solve(eqns.dBi==0, states.Bs);

for q=1:10
     eqns.(symbolicnames.eqns(q))=subs(eqns.(symbolicnames.eqns(q)), states.Bs, stars.Bs);   
 stars.(symbolicnames.states(q))=subs(stars.(symbolicnames.states(q)), states.Bs, stars.Bs);
end

%solve dLs==0 for Ps bc Ps show up the least times and its the
%smallest soln to solve for Ps
stars.Ps=solve(eqns.dLs==0, states.Ps);
for q=1:10
    eqns.(symbolicnames.eqns(q))=subs(eqns.(symbolicnames.eqns(q)), states.Ps, stars.Ps);   
 stars.(symbolicnames.states(q))=subs(stars.(symbolicnames.states(q)), states.Ps, stars.Ps);
end

%okay now we have eqns 1, 2,3 and 5 with vars D Ls Ms and Bi
%s solve for Bi using dBi. it
%comes with some conitions which are probably true
stars.Bi=solve(eqns.dBs==0, states.Bi);
for q=1:10
    eqns.(symbolicnames.eqns(q))=subs(eqns.(symbolicnames.eqns(q)), states.Bi, stars.Bi);   
 stars.(symbolicnames.states(q))=subs(stars.(symbolicnames.states(q)), states.Bi, stars.Bi);
end
%down to eqns 1-3 with vars D Ms and Ls. lets see whats easy.

% solve dD  for Ms.
stars.Ms=solve(eqns.dD==0, states.Ms);
for q=1:10
    eqns.(symbolicnames.eqns(q))=subs(eqns.(symbolicnames.eqns(q)), states.Ms, stars.Ms);   
 stars.(symbolicnames.states(q))=subs(stars.(symbolicnames.states(q)), states.Ms, stars.Ms);
end

%remaining: eqns 2 and 3 and vars Ls and D. solve for ls in dMs. We cannot
%sovle either for any state variable. But, we can factor each and retreive
%the factor which must be equal to zero. (note: I call these relation1 and
%relation2).
ceqns3_factor=factor(eqns.dMs); ceqns2_factor=factor(eqns.dPs);
relation1=ceqns3_factor(2); relation2=ceqns2_factor(1);
clear ceqns3_factor ceqns2_factor
%note: I know which root to choose by inspection
%for both, the following must be non n zero
    %D*pars.g_m*pars.gamma + D*pars.g_m*pars.mu_d - Ls*pars.g_m*pars.r_l - Ls*pars.mu_l*pars.r_m
    %ie Ls not equal to D*pars.g_m*pars.gamma + D*pars.g_m*pars.mu_d)/(pars.g_m*pars.r_l + pars.mu_l*pars.r_m)

% Set up jacobian.
Jendemic=J;
Jendemic(6:8,10)=zeros(3,1);
Jendemic=subs(Jendemic, states.Bi+states.Pi+states.Mi+states.Li,states.A*pars.v_0);
%plug all the algebraic expressions for stars into J endemic so that it
%only contain Ls and Ps following the backwards sovleorder
solveorder=[10;8;7;6;5;2;9;3];
for k=1:length(solveorder)
   var=solveorder(k);      
   Jendemic=subs(Jendemic, states.(symbolicnames.states(var)), stars.(symbolicnames.states(var)));
end
%%now the equilibirum stars, the changed DEs, and jacobian are only in terms of parametrs and D and Ls 
% eigenvalues cannot be solved rn so we have to sub in numbers and go from
% there. 
clear solveorder q qq i j ceqns var k 
%% 4: Initiate parameter sweeps
%step 1: choose number of parameter combos to consider
    vary.numberofcombos=1000 %10000 for asmpyt;
%step 2: chose a subset of parameters to vary.
    vary.type="all";
    %vary.type="small";
    %vary.type="all_nopars.gamma"
%ode.initate base params (for when unvaried)
params_refigured_simplified ;   basicparams=par';
 clear par
 %sweepkind="asymp"
 sweepkind="standard";
vary.scaleoveralllow=1.0000e-02; %changed 3/16 orig 1.0000e-04
vary.scaleoverallhigh=100;
vary.scalelowDISEASE=vary.scaleoveralllow;
vary.scalehighDISEASE=vary.scaleoverallhigh;
if sweepkind=="asymp"
  
    vary.scalelowTREEGOOD=vary.scaleoveralllow;
    vary.scalehighTREEGOOD=vary.scaleoverallhigh/10;
    vary.scalelowTREEBAD=vary.scaleoveralllow*10;
    vary.scalehighTREEBAD=vary.scaleoverallhigh;
else
    vary.scalelowTREEGOOD=vary.scaleoveralllow;
    vary.scalehighTREEGOOD=vary.scaleoverallhigh;
    vary.scalelowTREEBAD=vary.scaleoveralllow;
    vary.scalehighTREEBAD=vary.scaleoverallhigh;
end   
vary.parnums=[2,3,4,7,10,11,14,15,17,23,27];
vary.parnumsTREEGOOD=[2,3,4];
vary.parnumsTREEBAD=[7,10,11];

vary.parnumsDISEASE=[14,15,17,23,27];
vary.param_info(vary.parnumsTREEBAD,:)=basicparams(vary.parnumsTREEBAD)*[vary.scalelowTREEBAD,vary.scalehighTREEBAD];
vary.param_info(vary.parnumsTREEGOOD,:)=basicparams(vary.parnumsTREEGOOD)*[vary.scalelowTREEGOOD,vary.scalehighTREEGOOD];
vary.param_info(vary.parnumsDISEASE,:)=basicparams(vary.parnumsDISEASE)*[vary.scalelowDISEASE,vary.scalehighDISEASE];


  %%%%expand the range!!!!!!!
%vary.param_info([1:14,18:27],1:2)= vary.param_info([1:14,18:27],1:2).*[.02,50];



vary.numvariedpars= length(vary.parnums);
  vary.param_info_condensed=vary.param_info( vary.parnums, : );
  %%%%%%%%%%%%%


  
  %do hypercube
  vary.zeroonevalues=lhsdesign(vary.numberofcombos,vary.numvariedpars)';
  
  %transform using min max
  
 transformedlsh=vary.param_info_condensed(:,2).*vary.zeroonevalues+  vary.param_info_condensed(:,1).*(1-vary.zeroonevalues);

 clear  daysperyear
  %add values diretly related to others so they covary
  allparcombo=zeros(27,vary.numberofcombos);
  allparcombo(vary.parnums,:)=transformedlsh;

%alter related quantities for small
%rs are related. pars.r_m=.9*pars.r_l (based on orignal basic params
    allparcombo(1,:)=.9*allparcombo(2,:);   
%gs are related   %gm=.4 pars.g_p
    allparcombo(5,:)=.4*allparcombo(4,:);         allparcombo(6,:)=allparcombo(4,:);
%betas are related
    allparcombo(13,:)=.25*allparcombo(14,:);      allparcombo(12,:)=.01*allparcombo(14,:); 
%rhos are realted 
    allparcombo(16,:)=allparcombo(15,:);  
if vary.type=="all"
 %mus are related %pars.mu_m=2pars.mu_l       $pars.mu_p=3 pars.mu_l
    allparcombo(9,:)=2*allparcombo(10,:);         allparcombo(8,:)=3*allparcombo(10,:);
%deltas are related 
    allparcombo(24,:)=allparcombo(23,:);          allparcombo(25,:)=allparcombo(23,:);
    allparcombo(26,:)=allparcombo(23,:); 
end
    %add basic params to any parameter not varied
allparcombo(~any(allparcombo,2),:)=repmat(basicparams(~any(allparcombo,2),1),1,vary.numberofcombos);

allparcombo=[allparcombo,basicparams]; %add default paramerers to end
clear  transformedlsh

%% 4b: Filter Out par combos which lead to extinciton w/o diseases
error1=[];
error2=[];
tic;
parcombo_filtered=allparcombo;
for eachcombo=1:vary.numberofcombos

%evaluate rboundary
rboundary_eval=vectorfreestuff.rboundary;


%r broundary is good evaluate gamma boundary
   for eachpar=[4,10]
     rboundary_eval=subs(rboundary_eval, pars.(symbolicnames.pars(eachpar)), allparcombo(eachpar,eachcombo));
   end
rboundary_eval=double(rboundary_eval);

filterdata(eachcombo).rdiff=allparcombo(2,eachcombo)-rboundary_eval;


if allparcombo(2,eachcombo)<rboundary_eval
  error1=[error1,eachcombo]; %#ok<AGROW>
  
    
    
else
    gammaboundary_eval=vectorfreestuff.gamma_boundary;
     for eachpar=1:10
     gammaboundary_eval=subs(gammaboundary_eval, pars.(symbolicnames.pars(eachpar)), allparcombo(eachpar,eachcombo));
     end
     gammaboundary_eval=double(gammaboundary_eval);

     
     filterdata(eachcombo).gdiff=allparcombo(3,eachcombo)-gammaboundary_eval;
     if allparcombo(3,eachcombo)<gammaboundary_eval
      error2=[error2,eachcombo]; %#ok<AGROW>
     end
     
end

end

length(union(error1,error2))

parcombo_filtered(:,union(error1,error2))=[];
toc;

close all
plotnum=0;
 figure(1)

for parnum=[2,3,4,7,10,11]
   plotnum=plotnum+1; 

 subplot(3,2,plotnum)
 nbins=10;
 bounds=linspace(vary.param_info(parnum,1),vary.param_info(parnum,2),nbins+1);
 h=histogram(parcombo_filtered(parnum,:),bounds);
title(symbolicnames.pars(parnum))
 ylim([0,ceil(vary.numberofcombos/nbins)])
xlim(vary.param_info(parnum,1:2))



E = size(parcombo_filtered,2)/nbins*ones(nbins,1); % expected value (equal for uniform dist)

[h,p,stats] = chi2gof(parcombo_filtered(parnum,:),'Expected',E,'Edges',bounds);
allp(plotnum)=p;
   

uniformdist=makedist('Uniform',vary.param_info(parnum,1),vary.param_info(parnum,2));
kstest(parcombo_filtered(parnum,:),'CDF', uniformdist);

end
allp;
size(parcombo_filtered,2)




plotnum=0;
 figure(2)

for parnum=[14,15,17,23,27]
   plotnum=plotnum+1; 

 subplot(3,2,plotnum)
 nbins=10;
 bounds=linspace(vary.param_info(parnum,1),vary.param_info(parnum,2),nbins+1);
 h=histogram(parcombo_filtered(parnum,:),bounds);
title(symbolicnames.pars(parnum))
 ylim([0,ceil(vary.numberofcombos/nbins)])
xlim(vary.param_info(parnum,1:2))



E = size(parcombo_filtered,2)/nbins*ones(nbins,1); % expected value (equal for uniform dist)

[h,p,stats] = chi2gof(parcombo_filtered(parnum,:),'Expected',E,'Edges',bounds);
allp(plotnum)=p;
   

uniformdist=makedist('Uniform',vary.param_info(parnum,1),vary.param_info(parnum,2));
kstest(parcombo_filtered(parnum,:),'CDF', uniformdist);

end

%% dont need this except for exploraiton

rboundary_eval=vectorfreestuff.rboundary;
for eachcombo=1:vary.numberofcombos

%r broundary is good evaluate gamma boundary
   for eachpar=[4,10]
     rboundary_eval=subs(rboundary_eval, pars.(symbolicnames.pars(eachpar)), allparcombo(eachpar,eachcombo));
   end
rboundary_eval=double(rboundary_eval);
res(eachcombo).rboundarydiff=allparcombo(2,eachcombo)-rboundary_eval;
res(eachcombo).rboundaryrat=allparcombo(2,eachcombo)/rboundary_eval;

  gammaboundary_eval=vectorfreestuff.gamma_boundary;
     for eachpar=1:10
     gammaboundary_eval=subs(gammaboundary_eval, pars.(symbolicnames.pars(eachpar)), allparcombo(eachpar,eachcombo));
     end
     gammaboundary_eval=double(gammaboundary_eval);
     res(eachcombo).gammaboundarydiff=allparcombo(3,eachcombo)-gammaboundary_eval;
     res(eachcombo).gammaboundaryrat=allparcombo(3,eachcombo)/gammaboundary_eval;

end

%% 5: Stability analysis for the sweep


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show that endemic equilibirum always exists. this invovles
%showing that relation has a postive root and that the solution to that
%positve root also makes the other equation zero. S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this takes a long time.

%requires parallel computing

tic;
parcombo=parcombo_filtered;
vary.numberofcombos=size(parcombo,2);


%run ode for N centuries years to get an idea of what Ls and D might be
    ode.years=100; ode.tmax=365*ode.years; ode.tspan=[0 ode.tmax];
ode.init.D_0=140;   ode.init.P_s0=80;   ode.init.M_s0=60;   ode.init.L_s0=20;   ode.init.B_s0=0;
ode.init.P_i0=0;    ode.init.M_i0=0;    ode.init.L_i0=1;    ode.init.B_i0=0;    ode.init.A_0=1;
ode.init.all=[ode.init.D_0; ode.init.P_s0; ode.init.M_s0; ode.init.L_s0; ode.init.B_s0; ode.init.P_i0; ode.init.M_i0; ode.init.L_i0; ode.init.B_i0; ode.init.A_0];
%ode.options=odeset('NonNegative',1:10);
ode.options=[];

%Y_end=zeros(10, vary.numberofcombos);
clear runs eachcombo
century=zeros(1,vary.numberofcombos);
fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,vary.numberofcombos) '\n\n']);
parfor eachcombo=1:vary.numberofcombos
  %  parfor eachcombo=80:80

  
  if sweepkind=="asymp"
[Y_end(:,eachcombo),century(1,eachcombo)]=getYend(ode.tspan, ode.init.all, ode.options, parcombo(:,eachcombo)', 0.00000001); %0.00001  0.00000001
  end
  
    
  if sweepkind=="standard"
[Y_end(:,eachcombo),century(1,eachcombo)]=getYend(ode.tspan, ode.init.all, ode.options, parcombo(:,eachcombo)', 0.00001); %0.00001  0.00000001
  end

fprintf('\b|\n');

end
toc;
%% idk what this is (9/23)
inithold=ode.init.all;
inithold(1)=300
options2 = odeset('AbsTol',1e-20,'RelTol',1e-15)
%eachcombo=338
eachcombo=233

[hold]=getYend(ode.tspan, ode.init.all, options2, parcombo(:,eachcombo)',0.000000000000000001) %0.00001  0.00000001 0.00000001


 %%   %5a Run the model to equilibiurm for all parcombo
     %think this is trash, updated above.

parcombo=parcombo_filtered;
vary.numberofcombos=size(parcombo,2);

tic;

%run ode for N centuries years to get an idea of what Ls and D might be
    ode.years=100; ode.tmax=365*ode.years; ode.tspan=[0 ode.tmax];
ode.init.D_0=140;   ode.init.P_s0=80;   ode.init.M_s0=60;   ode.init.L_s0=20;   ode.init.B_s0=0;
ode.init.P_i0=0;    ode.init.M_i0=0;    ode.init.L_i0=1;    ode.init.B_i0=0;    ode.init.A_0=1;
ode.init.all=[ode.init.D_0; ode.init.P_s0; ode.init.M_s0; ode.init.L_s0; ode.init.B_s0; ode.init.P_i0; ode.init.M_i0; ode.init.L_i0; ode.init.B_i0; ode.init.A_0];
ode.options=[];
%Y_end=zeros(10, vary.numberofcombos);
clear runs eachcombo
runs(vary.numberofcombos)=struct();
tic;
parfor eachcombo=1:vary.numberofcombos
    runs(eachcombo).Y_end=zeros(10,1);
    %  if mod(eachcombo,5)==0
    %    eachcombo
    %   end
   %run for one centry
 [runs(eachcombo).T,runs(eachcombo).Y]=ode45(@Laurel_Model1_Equations, ode.tspan, ode.init.all, ode.options, parcombo(:,eachcombo)');      
 runs(eachcombo).Y_end=runs(eachcombo).Y( size(runs(eachcombo).Y,1),:);
 runs(eachcombo).century=1;
 runs(eachcombo).done="no";
% %if seed and large tree populations havent stabilized, run for another
% %century
  while runs(eachcombo).done=="no"
      [~,runs(eachcombo).fiftyyearidx]=min(abs(runs(eachcombo).T-365*50));
  if max(max(runs(eachcombo).Y(runs(eachcombo).fiftyyearidx:length(runs(eachcombo).Y),[1,4])-runs(eachcombo).Y_end([1,4],1)'))>0.01
     [runs(eachcombo).T,runs(eachcombo).Y]=ode45(@Laurel_Model1_Equations, ode.tspan, runs(eachcombo).Y_end(:,eachcombo), ode.options, parcombo(:,eachcombo)');      
     runs(eachcombo).Y_end=runs(eachcombo).Y( size(runs(eachcombo).Y,1),:);
   runs(eachcombo).century=runs(eachcombo).century+1;
      if runs(eachcombo).century>100
         runs(eachcombo).done="quit";
%       
      end
      else
      runs(eachcombo).done="yes";
  end  
  end

end
toc;

%threshold of 0.010 takes about 0.12 sec per combo. 
%threshold of 0.005 takes about 0.55 sec per combo. 

%% 5b: Evaluate
%clear findintersection
tic;
nowbad=zeros(vary.numberofcombos,1);
%Decision 1: How fine should the grid be resolved?
    decisions.numgrids=[6,16,62];
    decisions.zoomfactor=[0.08,0.02,0.005];



  % eachcombo=185;
%initiate the sturcure res (formally findintersection (But that got way to long to type)  
clear res
res(vary.numberofcombos)=struct();

fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,vary.numberofcombos) '\n\n']);


parfor eachcombo=1:vary.numberofcombos
%parfor eachcombo=80:80


 %   if ismember(eachcombo,bad)
     trynum=1;
    thiscombonumgrids=decisions.numgrids(1);
     
       
while trynum<4
    
     thiscombonumgrids=decisions.numgrids(trynum);

%Look at Ls close to  the equilibium given by the ode solver
trynum=3
res(eachcombo).Lsspace=linspace(Y_end(4, eachcombo)*(1-decisions.zoomfactor(trynum)),Y_end(4, eachcombo)*(1+decisions.zoomfactor(trynum)),thiscombonumgrids);

%we now have the D roots for each relation.NOTE: I'm assuming we always get 3 roots. This should be true
%because both relation1 and relation2 have D^3 as highest power of D. Can't get 4 roots 

%solve for the roots of both relations at the defined Ls for the combo
[res(eachcombo).Dsolve1,res(eachcombo).Dsolve2]=evaluateRelations(res(eachcombo).Lsspace, parcombo(:,eachcombo),relation1, relation2, pars,states,symbolicnames);


%keep track of when solve switched order so we can check for issues later
    if sum(abs(sum(sign(res(eachcombo).Dsolve2)))==size(res(eachcombo).Dsolve2,1))==1
          nowbad(eachcombo)=1;
    end


[res(eachcombo).intersectLs, res(eachcombo).intersectD,res(eachcombo).switchidx,res(eachcombo).branchidx]=getIntersectionPoint(res(eachcombo).Dsolve1,res(eachcombo).Dsolve2,  res(eachcombo).Lsspace,parcombo(:,eachcombo),relation1, relation2, pars,states,symbolicnames);
%
if res(eachcombo).switchidx< 0
    trynum=5;
else

 %calculate distance of intersection(s) from Y_end       
      res(eachcombo).Lsdist=res(eachcombo).intersectLs-Y_end(4,eachcombo);
      res(eachcombo).Ddist=res(eachcombo).intersectD-Y_end(1,eachcombo);
      res(eachcombo).euclidian=sqrt(( res(eachcombo).Lsdist).^2+( res(eachcombo).Ddist).^2);
      if min(res(eachcombo).euclidian)<0.01
          res(eachcombo).totaltries=trynum;
          trynum=4
      else
          trynum=trynum+1;
      end
end  
end %end try
    if length(res(eachcombo).intersectLs)>1
      [~,res(eachcombo).bestintersection]=min(res(eachcombo).euclidian);
      %save long list elsewhere
      res(eachcombo).intersectLsall=res(eachcombo).intersectLs;
      res(eachcombo).intersectDall=res(eachcombo).intersectD;
    %  res(eachcombo).switchidxall=res(eachcombo).switchidx;
      res(eachcombo).branchidxall=res(eachcombo).branchidx;
      res(eachcombo).euclidianall=res(eachcombo).euclidian;

      res(eachcombo).intersectLs=res(eachcombo).intersectLsall(res(eachcombo).bestintersection);
      res(eachcombo).intersectD=res(eachcombo).intersectDall(res(eachcombo).bestintersection);
     % res(eachcombo).switchidx= res(eachcombo).switchidxall(res(eachcombo).bestintersection);
      res(eachcombo).branchidx=res(eachcombo).branchidxall(:,res(eachcombo).bestintersection);  
      res(eachcombo).euclidian=res(eachcombo).euclidianall(res(eachcombo).bestintersection);
    end
    
  %  end  
   fprintf('\b|\n');
end % end for each combo

toc;
clear relation1_eval relation2_eval relation1_evalL relation2_evalL Y T maxL 
clear L1 L2 j  eachcombo eachpar idx %Dsolve1 Dsolve2 Lsspace

%% 5c: Check for a common issue I couldn't figure out how to fix in general

%sometimes, the function solve switches the order in which we get roots. if
%that happens around the intersection point, then we dont get the right
%value. We inspect the "bad" ones and look for that issue. When we find it,
%we recalculate  the intersection point

%for eachcombo=1:vary.numberofcombos
%if ismember(eachcombo,nowbad)
 eachcombo=338
 % eachcombo=233

spec1=["kx-","k*-", "k+-"] 
spec2=["rx-","r*-", "r+-"] 

      figure(100+eachcombo)
hold off
for branch=1:3
  %view the intersection plot

plot(res(eachcombo).Lsspace, res(eachcombo).Dsolve1(:,branch),convertStringsToChars(spec1(branch)), 'LineWidth',2)
   hold on
plot(res(eachcombo).Lsspace, res(eachcombo).Dsolve2(:,branch),spec2(branch))

   plot(Y_end(4,eachcombo), Y_end(1,eachcombo), '*b', 'MarkerSize', 10)
  plot(res(eachcombo).intersectLs, res(eachcombo).intersectD, '*g', 'MarkerSize', 10)

end
%end
%end

%% 5d Fix the common issue
badbad=[21,283,697] %these are the "bad" ones that need to be fixed
eachcombo=badbad(3)% FILL THIS IN
Lsspace=res(eachcombo).Lsspace;


%relation 1 and 2 come from cequns3 and ceqns 2=0. We want to find values
%of D and Ls that satisfy each relation simultaneiously 
    relation1_eval=relation1;    relation2_eval=relation2;
%plug in parameters into relations
    for eachpar=1:27
    relation1_eval=subs(relation1_eval, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
    relation2_eval=subs(relation2_eval, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
    end  
clear   Dsolve1 Dsolve2

for j=1:length(Lsspace)
%Sub Ls into things 
    relation1_evalL(j)=subs(relation1_eval, states.Ls, Lsspace(j));
    relation2_evalL(j)=subs(relation2_eval, states.Ls, Lsspace(j));
end %end each entry of Lsspace


 %evaluate rel1 and rel2 at Ds at 1/3 and 1/3 of switch idx, take mean
          Lsspacezoom=[Lsspace(eachidx:eachidx+1)*[1/3; 2/3], Lsspace(eachidx:eachidx+1)*[2/3; 1/3]];
          
            for j=1:2
            %Sub Ls into things 
                relation1_evalLzoom(j)=subs(relation1_eval, states.Ls, Lsspacezoom(j));
                relation2_evalLzoom(j)=subs(relation2_eval, states.Ls, Lsspacezoom(j));
            %Compute the value of D which makes each relation zero when Ls is from Lsspace
                Dsolve1zoom(j,:)=double(solve(relation1_evalLzoom(j)==0, states.D));
                Dsolve2zoom(j,:)=double(solve(relation2_evalLzoom(j)==0, states.D));
            end %end each entry of Lsspace
%       

    
      %21 okay
      % 283
            res(eachcombo).intersectD = mean(Dsolve2zoom(:,3));
            
     %697
            res(eachcombo).intersectD =mean([Dsolve2zoom(2,2),Dsolve2zoom(1,3)])
            
%% zoom in more (when euclidian distance is big). do this manually 
% use find([res.euclidian]>0.001) to determine which values need to be
% individually inpsected
% we did this for 94 111         383  375 334 291 85 215  167 162 312 find([res.euclidian]>0.001)

%366 59 291 520 521 499 233 380 252 497 562 563 48 415 435 54 427 52 485
%158 341 48 134 162 174 182 413 415 417 457 521 677 735 521 366

        decisions.zoomfactor=[0.08,0.02,0.0005]; %0.00005 for 111 and 0.0005 otherwise
        thiscombonumgrids=462;
eachcombo=366  % 520 %define each combo
%Look at Ls close to  the equilibium given by the ode solver
trynum=3;
res(eachcombo).Lsspace=linspace(Y_end(4, eachcombo)*(1-decisions.zoomfactor(trynum)),Y_end(4, eachcombo)*(1+decisions.zoomfactor(trynum)),thiscombonumgrids);

%we now have the D roots for each relation.NOTE: I'm assuming we always get 3 roots. This should be true
%because both relation1 and relation2 have D^3 as highest power of D. Can't get 4 roots 

%solve for the roots of both relations at the defined Ls for the combo
[res(eachcombo).Dsolve1,res(eachcombo).Dsolve2]=evaluateRelations(res(eachcombo).Lsspace, parcombo(:,eachcombo),relation1, relation2, pars,states,symbolicnames);


%keep track of when solve switched order so we can check for issues later
    if sum(abs(sum(sign(res(eachcombo).Dsolve2)))==size(res(eachcombo).Dsolve2,1))==1
          nowbad(eachcombo)=1;
    end


[res(eachcombo).intersectLs, res(eachcombo).intersectD,res(eachcombo).switchidx,res(eachcombo).branchidx]=getIntersectionPoint(res(eachcombo).Dsolve1,res(eachcombo).Dsolve2,  res(eachcombo).Lsspace,parcombo(:,eachcombo),relation1, relation2, pars,states,symbolicnames);
%
if res(eachcombo).switchidx< 0
    trynum=5;
else

 %calculate distance of intersection(s) from Y_end       
      res(eachcombo).Lsdist=res(eachcombo).intersectLs-Y_end(4,eachcombo);
      res(eachcombo).Ddist=res(eachcombo).intersectD-Y_end(1,eachcombo);
      res(eachcombo).euclidian=sqrt(( res(eachcombo).Lsdist).^2+( res(eachcombo).Ddist).^2);
      if min(res(eachcombo).euclidian)<0.01
          res(eachcombo).totaltries=trynum;
          trynum=4
      else
          trynum=trynum+1;
      end
end 

holdeuclidian= [res.euclidian]; 
    [a,b]=max(holdeuclidian);
    
    holdeuclidian(b)=0;
     [c,d]=max(holdeuclidian)
       holdeuclidian(d)=0;
     [e,f]=max(holdeuclidian);
  
   res(eachcombo).euclidian        
      
%%
%find where i didn't get an intersection 
noint=[];
for i=1:vary.numberofcombos
    
  if  isempty(res(i).switchidx)
      noint=[noint,i];
  end
end



%%
%evalutate relation1 and realtion 2 around Y_end(80) and subtract them. show neg on one side and pos on other



%%


%% 5d:Executive Summary

tic;
if length([res.switchidx])==vary.numberofcombos
    disp('Intrsection Report: intersection found for all combos');
else 
    disp('Intrsection Report: intersection NOT found in at least one spot');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Find max euclidian distance from (Dinf,Linf) and (D*,L*) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  parfor eachcombo=1:vary.numberofcombos
      
   %   Lsspace=linspace(Y_end(4, eachcombo)*.5,Y_end(4, eachcombo)*1.5,decisions.numgrids);
%findintersection(eachcombo).Lsspace=Lsspace;

    %  res(eachcombo).Lsdist=res(eachcombo).intersectLs-Y_end(4,eachcombo);
      
     % res(eachcombo).Ddist=res(eachcombo).intersectD-Y_end(1,eachcombo);
     % res(eachcombo).euclidian=sqrt((res(eachcombo).Lsdist)^2+(res(eachcombo).Ddist)^2);
      
      
      res(eachcombo).Dpercentdiff=abs(res(eachcombo).intersectD-Y_end(1,eachcombo))/((res(eachcombo).intersectD+Y_end(1,eachcombo))/2);
      res(eachcombo).lspercentdiff=abs(res(eachcombo).intersectLs-Y_end(4,eachcombo))/((res(eachcombo).intersectLs+Y_end(4,eachcombo))/2);

      
  end
  
  disp('Closeness Report: the max euclidian distance is')
  max([res.euclidian])
  
holdeuclidian= [res.euclidian]; 
    [a,b]=max(holdeuclidian)
    
    holdeuclidian(b)=0;
     [c,d]=max(holdeuclidian)
       holdeuclidian(d)=0;
     [e,f]=max(holdeuclidian)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Evaluate relations at the computed intersection point %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor eachcombo=1:vary.numberofcombos
 
    relation1_eval=relation1;    relation2_eval=relation2;
%plug in parameters into relations
    for eachpar=1:27
    relation1_eval=subs(relation1_eval, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
    relation2_eval=subs(relation2_eval, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
    end  
    
    relation1_eval=subs(relation1_eval, states.D, Y_end(1,eachcombo));
    relation1_eval=subs(relation1_eval, states.Ls, Y_end(4,eachcombo));
    relation2_eval=subs(relation2_eval, states.D, Y_end(1,eachcombo));
    relation2_eval=subs(relation2_eval, states.Ls, Y_end(4,eachcombo));
    
    res(eachcombo).rel1atYend=double(relation1_eval);
    res(eachcombo).re21atYend=double(relation2_eval);
    
end
disp('Relation Report: The relationsevaluated at the par is at most') 
max( max([res.rel1atYend]), max([res.re21atYend]))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Eigenvalues: Show that the eigenvalues of the endemic system are always
%negative. (Steal code from compare BS and no BS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   %sub parcombos and the values of D and Ls into J endemic
%for eachcombo=1:vary.numberofcombos
clear eigeninfo
eigeninfo(vary.numberofcombos)=struct();
parfor eachcombo=1:vary.numberofcombos
    tempJ_par=Jendemic;
    for p=1:27
    tempJ_par=subs(tempJ_par, pars.(symbolicnames.pars(p)), parcombo(p, eachcombo));
    end
    
    tempJ_inf=subs(tempJ_par, states.D, Y_end(1,eachcombo));
    tempJ_inf=subs(tempJ_inf, states.Ls, Y_end(4,eachcombo));
    eigeninfo(eachcombo).jacobian_inf=double(tempJ_inf);
    eigeninfo(eachcombo).eigenvalues_inf=  eig(eigeninfo(eachcombo).jacobian_inf);
    eigeninfo(eachcombo).neg_inf=all(real(eigeninfo(eachcombo).eigenvalues_inf)<0);
    
    tempJ_int=subs(tempJ_par, states.D, res(eachcombo).intersectD);
    tempJ_int=subs(tempJ_int, states.Ls,res(eachcombo).intersectLs);
    eigeninfo(eachcombo).jacobian_int=double(tempJ_int);
    eigeninfo(eachcombo).eigenvalues_int=  eig(eigeninfo(eachcombo).jacobian_int);
    eigeninfo(eachcombo).neg_int=all(real(eigeninfo(eachcombo).eigenvalues_int)<0);
  
end

if sum([eigeninfo.neg_inf])==vary.numberofcombos
    display('Stability Report: stable throughout all region (inf)')
else 
    display('Stability Report: unstable somewhere (inf)')

end

if sum([eigeninfo.neg_int])==vary.numberofcombos
    display('Stability Report: stable throughout all region (int)')
else 
    display('Stability Report: unstable somewhere (int)')

end
toc;





%% this has to do with the MVT proof. not fully flushed out but not needed
if objective==2
    
if decisions.computecoef=="yes"
%Now I want to demonstrate that relation 1 and 2 are both postive for Ls
%small enough and for some D. This is necessary for IVT

%we need the coefficents of a and b. this takes a while to get.
    powersofD=children(collect(relation1, states.D));
    D2powersofL=children(collect(powersofD{2},states.Ls));
    D1powersofL=children(collect(powersofD{3},states.Ls));
    D0powersofL=children(collect(powersofD{4},states.Ls));

    %relation 1 terms (b coefs)
        coef1.D3L0=simplify(powersofD{1,1}); %D^3    short    positive
        coef1.D2L1=simplify(D2powersofL{1,1});%D^2 L long    unknown
        coef1.D2L0=simplify(D2powersofL{1,2});%D^2 long   negative
        coef1.D1L2=simplify(D1powersofL{1,1});%D Ls^2 long unknown
        coef1.D1L1=simplify(D1powersofL{1,2});%D LS short   positive
        coef1.D0L3=simplify(D0powersofL{1,1});% Ls^3 long uninown
        coef1.D0L2=simplify(D0powersofL{1,2} );% LS^2short  negative

    powersofD=children(collect(relation2, states.D));
    D3powersofL=children(collect(powersofD{1},states.Ls))
    D2powersofL=children(collect(powersofD{2},states.Ls))
    D1powersofL=children(collect(powersofD{3},states.Ls))
    D0powersofL=children(collect(powersofD{4},states.Ls))
    %relation 2 terms (a coefs)
        coef2.D3L1=simplify(D3powersofL{1,1}); %D^3  Ls short   positive
        coef2.D3L0=simplify(D3powersofL{1,2}); %D^3     short   negative
        coef2.D2L2=simplify(D2powersofL{1,1});%D^2 L^2  short   positive
        coef2.D2L1=simplify(D2powersofL{1,2});%D^2 L    short   unknown
        %^sign depends on 3*pars.beta_p*pars.g_m*pars.r_l - pars.beta_l*pars.g_p*pars.r_m - pars.beta_l*pars.mu_p*pars.r_m + 3*pars.beta_p*pars.mu_l*pars.r_m
        coef2.D1L3=simplify(D1powersofL{1,1});%D Ls^3   long    unknown
        coef2.D1L2=simplify(D1powersofL{1,2});%D LS^2   short   unknown
        %^sign depends on 2*pars.beta_l*pars.g_p*pars.r_m - 3*pars.beta_p*pars.g_m*pars.r_l + 2*pars.beta_l*pars.mu_p*pars.r_m - 3*pars.beta_p*pars.mu_l*pars.r_m
        coef2.D0L4=simplify(D0powersofL{1,1});% Ls^4    short    unknown
        coef2.D0L3=simplify(D0powersofL{1,2} );% LS32   short   unknown
        %^* sign depends on pars.beta_p*pars.g_m*pars.r_l - pars.beta_l*pars.g_p*pars.r_m - pars.beta_l*pars.mu_p*pars.r_m + pars.beta_p*pars.mu_l*pars.r_m
clear powersofD D2powersofL D1powersofL D0powersofL D3powersofL

%Now we want to find alpha which is the min of alpha 1 and alpha 2
%plug in parameters for a_04 and b03 
for eachcombo=1:5
       coef1.eval.D0L3=subs(coef1.D0L3, states.Ls, 1);
       coef1.eval.D0L2=subs(coef1.D0L2, states.Ls, 1);
       coef2.eval.D0L4=subs(coef2.D0L4, states.Ls,1);
       coef2.eval.D0L3=subs(coef2.D0L3, states.Ls,1);
  for eachpar=1:27
      coef1.eval.D0L3=subs(coef1.eval.D0L3, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
      coef2.eval.D0L4=subs(coef2.eval.D0L4, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
      coef1.eval.D0L2=subs(coef1.eval.D0L2, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
      coef2.eval.D0L3=subs(coef2.eval.D0L3, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));

  end  
holdsign(1,eachcombo)=sign(coef1.eval.D0L3);
holdsign(2,eachcombo)=sign(coef2.eval.D0L4);
holdsign(3,eachcombo)=sign(coef1.eval.D0L2);
holdsign(4,eachcombo)=sign(coef2.eval.D0L3);

alpha(1,eachcombo)=abs(double(coef1.eval.D0L2)/double(coef1.eval.D0L3));
alpha(2,eachcombo)=abs(double(coef2.eval.D0L3/coef2.eval.D0L4));

end

%let Ls vary from 0 to alpha and evaluate each reation at a D value
for eachcombo=1:75
    relation1_eval=relation1;
    relation2_eval=relation2;
   % Lssmall=linspace(tol, 2*min(alpha(:,eachcombo)),30);
   Lssmall=linspace(tol, .2,25);
    for eachpar=1:27
      relation1_eval=subs(relation1_eval, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
      relation2_eval=subs(relation2_eval, pars.(symbolicnames.pars(eachpar)), parcombo(eachpar,eachcombo));
    
    end  
    %sub in D value 
     relation1_eval=subs(relation1_eval, states.D, .2);
     relation2_eval=subs(relation2_eval, states.D, .2);
     %sub in Lssmall
     for eachLs=1:20
     relationsave1(eachLs, eachcombo)=subs(relation1_eval, states.Ls, Lssmall(1,eachLs));
     relationsave2(eachLs, eachcombo)=subs(relation2_eval, states.Ls, Lssmall(1,eachLs));
     end
    
end
relationsave1=double(relationsave1)
relationsave2=double(relationsave2)

end





end%end objective

toc
  %%

function [Y_end,century]=getYend(tspan, init, opt,pars,threshold)
[T,Y]=ode45(@Laurel_Model1_Equations, tspan, init, opt, pars);      
Y_end=Y( size(Y,1),:);
century=1;
done="no";
%if seed and large tree populations havent stabilized, run for another
%century
while done=="no"
    [~,fiftyyearidx]=min(abs(T-365*50));
if max(max(Y(fiftyyearidx:length(Y),[1,4])-Y_end([1,4])))>threshold
   [T,Y]=ode45(@Laurel_Model1_Equations, tspan, Y_end, opt, pars);      
   Y_end=Y( size(Y,1),:);
 century=century+1;
    if century>100
        done="quit";
    
    end
    else
    done="yes";
end  
end
end


function [Dsolve1,Dsolve2]=evaluateRelations(Lsspace, parvalues,relation1_eval, relation2_eval, pars,states,symbolicnames)
%this function takes in a set of Ls values and a parameter combo, and it
%evaluates both relations at them and sovles for D. It returns the set of D

%plug in parameters into relations
    for eachpar=1:27
    relation1_eval=subs(relation1_eval, pars.(symbolicnames.pars(eachpar)), parvalues(eachpar));
    relation2_eval=subs(relation2_eval, pars.(symbolicnames.pars(eachpar)), parvalues(eachpar));
    end  
relation1_evalL=sym(zeros(length(Lsspace),1));
relation2_evalL=sym(zeros(length(Lsspace),1));
Dsolve1=zeros(length(Lsspace),3);
Dsolve2=zeros(length(Lsspace),3);

for j=1:length(Lsspace)
%Sub Ls into things 
    relation1_evalL(j)=subs(relation1_eval, states.Ls, Lsspace(j));
    relation2_evalL(j)=subs(relation2_eval, states.Ls, Lsspace(j));
%Compute the value of D which makes each relation zero when Ls is from Lsspace
    Dsolve1(j,:)=double(solve(relation1_evalL(j)==0, states.D));
    Dsolve2(j,:)=double(solve(relation2_evalL(j)==0, states.D));
end %end each entry of Lsspace

end



function [intersectLs, intersectD,switchidx,branchidx]=getIntersectionPoint(Dsolve1,Dsolve2,  Lsspace,parvalue,relation1, relation2, pars,states,symbolicnames)
    %it needs Dsolve1 and Dsovle2 in addition to everuthing that
    %evaluteRelations needs (since it calls evaluateRelations)
branchA_opt=sum(Dsolve1>0)>0 & sum(real(Dsolve1)~=0);
branchB_opt=sum(Dsolve2>0)>0 & sum(real(Dsolve2)~=0);
%d%%decide what to dow ith this
     intersectLs=[];  
     intersectD=[];
     switchidx=[];
     finalswitchidx=[];
     branchidx=[];
for brancha=1:3
    if branchA_opt(brancha)
    
    for branchb=1:3
        if branchB_opt(branchb)
            diffD=Dsolve1(:,brancha)-Dsolve2(:,branchb);

  if ~or(all(diffD>0), all(diffD<0)) %if there's a sign change, it's a candidate for intersection
      %find the location of the sign change: multiply each entry of diff by the next entry of diff, and find when
      %the product is negative 
        switchidx=find(real(diffD(1:length(diffD)-1)).*real(diffD(2:length(diffD)))<0);
      
      %Remove switch idx where both sides are imaginary
      switchidx= union(intersect(switchidx,find(diffD==real(diffD))),intersect(switchidx+1,find(diffD==real(diffD)))-1);
  
      if ~isempty(switchidx)
      for i=1:length(switchidx)
          eachidx=switchidx(i);
      %evaluate rel1 and rel2 at Ds at 1/3 and 1/3 of switch idx, take mean
          Lsspacezoom=[Lsspace(eachidx:eachidx+1)*[1/3; 2/3], Lsspace(eachidx:eachidx+1)*[2/3; 1/3]];
[~,Dsolve2zoom]=evaluateRelations(Lsspacezoom, parvalue,relation1, relation2, pars,states,symbolicnames);
  
%replaced with evaluateRelations
        %%%%%%%%%%%%%%%%
           % for j=1:2
            %Sub Ls into things 
              %  relation1_evalLzoom(j)=subs(relation1_eval, states.Ls, Lsspacezoom(j));
             %   relation2_evalLzoom(j)=subs(relation2_eval, states.Ls, Lsspacezoom(j));
            %Compute the value of D which makes each relation zero when Ls is from Lsspace
             %   Dsolve1zoom(j,:)=double(solve(relation1_evalLzoom(j)==0, states.D));
              %  Dsolve2zoom(j,:)=double(solve(relation2_evalLzoom(j)==0, states.D));
           % end %end each entry of Lsspace
%       %%%%%%%%%%%%%%%%%%%
 
      intersectLs=[intersectLs,mean(Lsspacezoom)]; %#ok<AGROW>
      intersectD =[intersectD, mean(Dsolve2zoom(:,branchb))]; %#ok<AGROW>
      finalswitchidx=[finalswitchidx,eachidx]; %#ok<AGROW>
      branchidx=[branchidx,[brancha;branchb]]; %#ok<AGROW>
      end
     
      end  
      
  end %end has a sign change
        end
    end %end branchb
    end
end %end brancha
if isempty(finalswitchidx)
    finalswitchidx=-1;
end

end %end function



