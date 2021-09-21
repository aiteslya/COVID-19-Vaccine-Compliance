% This is the script that outputs figures showing the sensitivity of outputs
% to the variations of initial data with respect to the initial percentage
% of exposed individuals, the rest of the epidemiological compartments
% remain with the same size, save for susceptible which gets adjusted such
% that the total population size adds up to the population size of the
% Netherlands
% keeping the rest of the parameters fixed
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 with Re(0)=1.1 at the start of the
% vaccination rollout
% the outputs are the cumulative number of infected, peak of prevalence and
% timing of the prevalence peak
% time units are in days
%prepare settings
clc;
clear variables;
close all;
format long;

%set the parameters

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

%set up initial data: necessary for the calculation of initial conditions,
%but will be varied in the main loop
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
%compliance waning rate when the vaccination coverage is zero
mu0=1/30;

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

figure(20);
plot(Res(:,1).*Res(:,2),Res(:,1),'o-','LineWidth',2);
xlabel({'Contact rate of compliant individuals';'(individuals/day)'},'interpreter','latex');
ylabel({'Contact rate of non-compliant';' individuals (individuals/day)'},'interpreter','latex');
xlim([0.9*min(Res(:,1).*Res(:,2)),1.01*max(Res(:,1).*Res(:,2))]);
set(gca,'FontSize',20);
%set up the array of initial average contacts:
%c*ProportionNonCompl+c*r1*ProportionCompl
temparr=Res(:,1)*(1-PerCompl)+Res(:,2).*Res(:,1)*PerCompl;
%find index where the average contact rate is closest to 5
temparr0=temparr-5;
eqn=temparr0(2:r1num).*temparr0(1:(r1num-1));
ind=find(eqn<0,1);
c=Res(ind,1);
r1=Res(ind,2);

%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=200;

k1=1;
k2=1;

%format of the legend
formatSpec = '%.2e';

% set up arrays
% vaccinated increase their contact rate to the levels of pre-pandemic
r2=13.47/c;
%fraction of population that has been vaccinated when the compliance state
%on average lasts 7 days
fracarr=1/3;
mu1=3./(fracarr*5.1e8);
% set up main arrays for the heat maps
meshsize=25;
minupsilon=5e-4;%5.9e-4;
maxupsilon=6e-3;%0.0053;% vaccination rate in Israel
upsilonarr=[minupsilon,maxupsilon];
% output the figure with minimum and maximum vaccination curves

figc=1;
counter=1;
tarr=linspace(0,200,100);
% numpoints=100;
% tdiscr=linspace(0,T,numpoints);
omega=0.6;
%define alpha and gamma arrays

delta=4e-5;

% original variant
beta=epsilon*c;

%outer loop: loops through the vaccination uptake rates
%figure counter
subfigc=1;
% define compliance proportions
ExpProp=linspace(0.001,0.01,meshsize);

%allocate result arrays
Cum03=zeros(1,meshsize);
Cum06=zeros(1,meshsize);
upsilon=0;
for i1=1:meshsize 
    % set the currently infectious individuals
    % calculate initial conditions
    %Incid=TotalInf*gamma;
    %seroprevalence
    %SP=0.08;
    %total population
    %N=1.7e7;

    %N0=N*(1-PerCompl);
    %Nc0=N*PerCompl;
    %TotalRec=SP*N;
    TotalE=ExpProp(i1)*N;
    TotalS=N-TotalInf-TotalRec-TotalE;
    %setting up of initial data
    S0=(1-PerCompl)*TotalS;
    %E0=(1-PerCompl)*TotalInf*gamma/alpha;
    E0=(1-PerCompl)*TotalE;
    %R0=(1-PerCompl)*TotalRec;
    Sc0=PerCompl*TotalS;
    %Ec0=PerCompl*TotalInf*gamma/alpha;
    Ec0=PerCompl*TotalE;
    %Rc0=PerCompl*TotalRec;

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
    %set parameters
    pars=[beta,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega];
    [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
    % collect and process the output
    % cumulative new infections from the start of the vaccination
    % rollout
    cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;
    %find indices of 3 and 6 months
    ind3=find(t>3*30,1);
    ind6=find(t>6*30,1);
    %cumulative infected
    cum3=100*cum(ind3)/N;
    cum6=100*cum(ind6)/N;
    Cum03(1,i1)=cum3;
    Cum06(1,i1)=cum6;
end
    
for upsilon=upsilonarr
    
    %allocate result arrays
    Cum3=zeros(1,meshsize);
    Cum6=zeros(1,meshsize);
    for i1=1:meshsize 
        % set the currently infectious individuals
        % calculate initial conditions
        %Incid=TotalInf*gamma;
        %seroprevalence
        %SP=0.08;
        %total population
        %N=1.7e7;
        
        %N0=N*(1-PerCompl);
        %Nc0=N*PerCompl;
        %TotalRec=SP*N;
        TotalE=ExpProp(i1)*N;
        TotalS=N-TotalInf-TotalRec-TotalE;
        %setting up of initial data
        S0=(1-PerCompl)*TotalS;
        %E0=(1-PerCompl)*TotalInf*gamma/alpha;
        E0=(1-PerCompl)*TotalE;
        %R0=(1-PerCompl)*TotalRec;
        Sc0=PerCompl*TotalS;
        %Ec0=PerCompl*TotalInf*gamma/alpha;
        Ec0=PerCompl*TotalE;
        %Rc0=PerCompl*TotalRec;

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
        %set parameters
        pars=[beta,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega];
        [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
        % collect and process the output
        % cumulative new infections from the start of the vaccination
        % rollout
        cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;
        %find indices of 3 and 6 months
        ind3=find(t>3*30,1);
        ind6=find(t>6*30,1);
        %cumulative infected
        cum3=100*cum(ind3)/N;
        cum6=100*cum(ind6)/N;
        Cum3(1,i1)=cum3;
        Cum6(1,i1)=cum6;
    end

    f1=figure(1);
    subplot(1,2,subfigc);
    
    h1(1)=plot(100*ExpProp,Cum3,'LineWidth',4);hold on;
    h1(2)=plot(100*ExpProp,Cum6,'-.','LineWidth',4);hold on;
    xlim([min(100*ExpProp),max(100*ExpProp)]);
    ylim([1,4]);
    xlabel({'Percentage of exposed individuals';'at the start of the vaccination rollout (\%)'},'interpreter','latex');

    if upsilon==upsilonarr(1) 
        title('Slow vaccination','interpreter','latex');
    else
        title('Fast vaccination','interpreter','latex');
    end
    set(gca,'FontSize',20);

    f2=figure(2);
    subplot(1,2,subfigc);
    
    h2(1)=plot(100*ExpProp,100*(Cum3-Cum03)./Cum03,'LineWidth',4);hold on;
    h2(2)=plot(100*ExpProp,100*(Cum6-Cum06)./Cum06,'-.','LineWidth',4);hold on;
    xlim([min(100*ExpProp),max(100*ExpProp)]);
    %ylim([1,4]);
    xlabel({'Percentage of exposed individuals';'at the start of the vaccination rollout (\%)'},'interpreter','latex');

    % output of the colorbar only for subplot 2,2,4
    if subfigc==2
        legend(h2,'3 months since the vaccination rollout','6 months since the vaccination rollout','Location','southeastoutside','NumColumns',2);
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

yl1=ylabel(han1,{'Cumulative number of new infections (\%)'},'interpreter','latex');
set(yl1,'position',[-0.07 0.500000476837158 0])
set(yl1,'FontSize',20);

%annotations
figure(1);annotation('textbox', [0.08, 0.99, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.53, 0.99, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)

han2=axes(f2,'visible','off');

han2.Title.Visible='on';
han2.XLabel.Visible='on';
han2.YLabel.Visible='on';

yl2=ylabel(han2,{'Relative difference in the cumulative';'number of new infections (\%)'},'interpreter','latex');
set(yl2,'position',[-0.07 0.500000476837158 0])
set(yl2,'FontSize',20);

%annotations
figure(2);annotation('textbox', [0.08, 0.99, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.53, 0.99, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)
