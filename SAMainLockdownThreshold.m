% Dynamics of SARS-CoV-2 transmission dynamics in the population where the
% vaccination rollout started and the compliance with physical distancing
% measures decreases as the vaccination coverage grows. The government
% implements a lockdown if the prevalence exceeds a threshold.
% This script performs the sensitivity analysis of the systems dynamics
% with respect to threshold value at which the lockdown is initiated and
% removed. Original variant circulating. Two vaccination uptake rates are
% considered.

% prepare state space
clc;
close all;
clear variables;
format long;
%format of the legend
formatSpec = '%.2e';

%prepare control parameter settings

meshsize=30;

pcarr=linspace(50,1000,meshsize)/1.7e7;

figc=1;

%Set parameters

%contact rate of non-compliant before the lockdown at the start of the
%epidemic
chat=14.9;
%alpha and gamma are fixed here to calculate c and r1 and initial data
%1/alpha duration of the latent period
alpha=1/4;
%1/gamma duration of the infectious period
gamma=1/7;
R0=2.5;
%calculate epsilon
epsilon=R0*gamma/chat;

%compliance waning rate when the vaccination coverage is zero
mu0=1/30;
delta=4e-5;

%fraction of population that has been vaccinated when the compliance state
%on average lasts 7 days
fracarr=1/3;
mu1=3./(fracarr*5.1e8);

minupsilon=5e-4;%5.9e-4;
maxupsilon=6e-3;%0.0053;% vaccination rate in Israel
upsilonarr=[minupsilon,maxupsilon];

%set vaccine efficacy
omega=0.6;

k1=1;
k2=1;

%Set initial conditions
%set up initial data: necessary for the calculation of initial conditions,
%prevalence of infectious cases
TotalInf=112435;
Incid=TotalInf*gamma;
%seroprevalence
SP=0.08;
%total population
N=1.7e7;
%percentage of compliant people
PerCompl=0.65;

N0=N*(1-PerCompl);
Nc0=N*PerCompl;
TotalRec=SP*N;
TotalE=TotalInf*gamma/alpha;
TotalS=N-TotalInf-TotalRec-TotalE;
%setting up of initial data
S0=(1-PerCompl)*TotalS;
E0=(1-PerCompl)*TotalInf*gamma/alpha;
I0=(1-PerCompl)*TotalInf;
R0=(1-PerCompl)*TotalRec;
Sc0=PerCompl*TotalS;
Ec0=PerCompl*TotalInf*gamma/alpha;
Ic0=PerCompl*TotalInf;
Rc0=PerCompl*TotalRec;

infect0=E0+I0+R0+Ec0+Ic0+Rc0;
popOut=1e5;
infect0pop=popOut*infect0/N;

V0=0;
SV0=0;
EV0=0;
IV0=0;
RV0=0;

TV0=V0+SV0+EV0+IV0+RV0;
init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0,0,0,0];

%calculation of the contact rate of non-compliant individuals c and the
%reduction factor of compliant individuals contact rate when compared to
%non-compliant
%define arrays of reduction factors and contact rates
r1num=60;
r1arr=linspace(0,1,r1num);
%define the array of r1
num=2e4;
carr=linspace(0,15,num);

%for each c in carr calculate r1
%set up equation Re-1.1=0

Res=nan(r1num,2);
counter=1;
for r1=r1arr
    beta=carr*epsilon;
    Re1=beta*S0./(gamma*(N0+r1*Nc0))+beta*r1*Sc0.*(mu0*(alpha+gamma+mu0)+alpha*gamma*r1)./(gamma*(alpha+mu0)*(gamma+mu0).*(N0+r1*Nc0))-1.1;
    nulleqn=Re1(1:num-1).*Re1(2:num);
    ind=find(nulleqn<0);
    if numel(ind)==1
        Res(counter,1)=carr(ind);
        Res(counter,2)=r1;
        counter=counter+1;
    elseif numel(ind)>0
        error('RcContact: more than one root');
    end
end

%set up the array of initial average contacts:
%c*ProportionNonCompl+c*r1*ProportionCompl
temparr=Res(:,1)*(1-PerCompl)+Res(:,2).*Res(:,1)*PerCompl;
%find index where the average contact rate is closest to 5
temparr0=temparr-5;
eqn=temparr0(2:r1num).*temparr0(1:(r1num-1));
ind=find(eqn<0,1);
c=Res(ind,1);
r1=Res(ind,2);
% vaccinated increase their contact rate to the levels of pre-pandemic
r2=13.47/c;
%contact rate in lockdown
clockD=3;

%integration options
Atol=1e-5;
opts = odeset('RelTol',1e-4,'AbsTol',Atol);
%integrating time
T=200;

subfigc=1;

%calculate the no-vaccination scenario
upsilon0=0;
beta=c*epsilon;
parsNoLock=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
[t0,y0]=ode23s(@(t,y)COVIDVaccineRHS2(t,y,parsNoLock),[0,T], init,opts);

cum0=(y0(:,2)+y0(:,3)+y0(:,4)+y0(:,6)+y0(:,7)+y0(:,8)+y0(:,11)+y0(:,12)+y0(:,13)+y0(:,15)+y0(:,16)+y0(:,17))-infect0;

%find indices of 3 and 6 months
ind03=find(t0>3*30,1);
ind06=find(t0>6*30,1);
%cumulative infected
Cum03=100*cum0(ind03)/N;
Cum06=100*cum0(ind06)/N;

for upsilon=upsilonarr
    %allocate result arrays
    Cum3=zeros(1,meshsize);
    Cum6=zeros(1,meshsize);
    MaxPrev=zeros(1,meshsize);
    for i1=1:meshsize
        pc=pcarr(i1);
                        
        pars=[clockD,c,epsilon,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega,pc];
        [t,y]=ode23s(@(t,y)COVIDVaccIntervenRHS(t,y,pars),[0,T], init,opts);
        pop=sum(y(:,1:17),2)-y(:,14);
        
        cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;
        prev=y(:,3)+y(:,7)+y(:,12)+y(:,16);
        
        %find indices of 3 and 6 months
        ind3=find(t>3*30,1);
        ind6=find(t>6*30,1);
        %cumulative infected
        cum3=100*cum(ind3)/N;
        cum6=100*cum(ind6)/N;
        Cum3(1,i1)=cum3;
        Cum6(1,i1)=cum6;
        MaxPrev(1,i1)=max(prev);
    end
    f1=figure(1);
    subplot(1,2,subfigc);
    
    h1(1)=plot(100*pcarr,Cum3,'LineWidth',4);hold on;
    h1(2)=plot(100*pcarr,Cum6,'-.','LineWidth',4);hold on;
    xlim([min(100*pcarr),max(100*pcarr)]);
    ylim([0.6,1.8]);
    xlabel({'Infectious individuals threshold';'for the commencement of the lockdown (\%)'},'interpreter','latex');

    % output of the colorbar only for subplot 2,2,4
    if subfigc==2
        legend(h1,'3 months since the vaccination rollout','6 months since the vaccination rollout','Location','southeastoutside','NumColumns',2);
    end
    
    if upsilon==upsilonarr(1) 
        title('Slow vaccination','interpreter','latex');
    else
        title('Fast vaccination','interpreter','latex');
    end
    set(gca,'FontSize',20);
    
    f2=figure(2);
    subplot(1,2,subfigc);
    
    h2(1)=plot(100*pcarr,100*(Cum3-Cum03)./Cum03,'LineWidth',4);hold on;
    h2(2)=plot(100*pcarr,100*(Cum6-Cum06)./Cum06,'-.','LineWidth',4);hold on;
    xlim([min(100*pcarr),max(100*pcarr)]);
    %ylim([0.6,1.8]);
    xlabel({'Infectious individuals threshold';'for the commencement of the lockdown (\%)'},'interpreter','latex');

    % output of the colorbar only for subplot 2,2,4
    if subfigc==2
        legend(h1,'3 months since the vaccination rollout','6 months since the vaccination rollout','Location','southeastoutside','NumColumns',2);
    end
    
    if upsilon==upsilonarr(1) 
        title('Slow vaccination','interpreter','latex');
    else
        title('Fast vaccination','interpreter','latex');
    end
    set(gca,'FontSize',20);
    
 
    subfigc=subfigc+1; 
end

han1=axes(f1,'visible','off');

han1.Title.Visible='on';
han1.XLabel.Visible='on';
han1.YLabel.Visible='on';

yl1=ylabel(han1,{'Relative difference in the cumulative';'number of new infections (\%)'},'interpreter','latex');
set(yl1,'position',[-0.07 0.500000476837158 0])
set(yl1,'FontSize',25);

%annotations
figure(1);annotation('textbox', [0.08, 0.99, 0, 0], 'string', 'a','FontWeight','bold','FontSize',30)
figure(1);annotation('textbox', [0.53, 0.99, 0, 0], 'string', 'b','FontWeight','bold','FontSize',30)

han2=axes(f2,'visible','off');

han2.Title.Visible='on';
han2.XLabel.Visible='on';
han2.YLabel.Visible='on';

yl2=ylabel(han2,{'Relative difference in the cumulative';'number of new infections (\%)'},'interpreter','latex');
set(yl2,'position',[-0.07 0.500000476837158 0])
set(yl2,'FontSize',25);

%annotations
figure(2);annotation('textbox', [0.08, 0.99, 0, 0], 'string', 'c','FontWeight','bold','FontSize',30)
figure(2);annotation('textbox', [0.53, 0.99, 0, 0], 'string', 'd','FontWeight','bold','FontSize',30)