%This is the script that plots heat maps for varying variant varying upsilon,
%omega and mu1 and r2, keeping the rest of the parameters fixed
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 with Re(0)=1.1 at the start of the
% vaccination rollout
% the outputs are the cumulative
% percentage of infected relative to no vaccination scenario adjusted to
% not count the infected and recovered at time 0
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
%1/gamma duration of infectious period
gamma=1/7;
alpha=1/4;
R0=2.5;
%calculate epsilon
epsilon=R0*gamma/chat;

%set up initial data
% factor of detected/total
X=1;
%prevalence of infectious cases
TotalInf=112435;%37706;
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
E0=(1-PerCompl)*TotalInf*gamma/alpha;%7541;
I0=(1-PerCompl)*TotalInf;% 13197;
R0=(1-PerCompl)*TotalRec;
Sc0=PerCompl*TotalS;
Ec0=PerCompl*TotalInf*gamma/alpha;%14005;
Ic0=PerCompl*TotalInf;%24509;%PerCompl*TotalInf/2;
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
%compliance
mu0=1/30;

%define array of contact rates
r1num=60;
r1arr=linspace(0,1,r1num);
%define the array of r1
num=2e4;
carr=linspace(0,15,num);

%for each c in carr calculate r1
%set up equation Re=1.1

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

%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time
T=200;

delta=4e-5;
k1=1;
k2=1;

%format of the legend
formatSpec = '%.2e';

% set up arrays
% vaccinated can increase their contact rate as high as non-compliant
% during the lockdown or even go up to the levels of pre-pandemic
r2arr=[chat/c,1];
%fraction of population that has been vaccinated when the compliance state
%on average lasts 7 days
fracarr=[1/3,2/3];
mu1arr=3./(fracarr*5.1e8);
% set up main arrays for the heat maps
meshsize=25;
minupsilon=5.9e-4;
maxupsilon=0.0053;% vaccination rate in Israel
upsilonarr=linspace(minupsilon,maxupsilon,meshsize);%[5.9e-4,4.9e-3];
% output the figure with minimum and maximum vaccination curves
%define colors for the vaccination rates

l=5;
red = [0.01, 0.4, 0.76];%[1 0 0];
pink = [0.45,0.76,0.98];
colors_p = [linspace(pink(1),red(1),l)', linspace(pink(2),red(2),l)', linspace(pink(3),red(3),l)'];
col_p=colors_p([1,5],:);
yearcol=[0.17,0.17,0.17];
novacc=[0.9,0.13,0.13];

figc=1;
counter=1;
tarr=linspace(0,200,100);

beta=epsilon*c;
r2=1;
numpoints=100;
tdiscr=linspace(0,T,numpoints);

% output data points for Belarusia and Israel
%read the CSV files and load them into array
TblBEL = readtable('DataVaccine/Belarus.csv');
TblISR = readtable('DataVaccine/Israel.csv');

% convert dates to serial numbers and set the first date in the date set as
% zero
dates=TblBEL.date;
DATEBEL=datenum(dates)-datenum(dates(1))+1;

dates=TblISR.date;
DATEISR=datenum(dates)-datenum(dates(1))+1;

TBLBEL=[DATEBEL';(TblBEL.people_vaccinated)';(TblBEL.people_fully_vaccinated)';]';
TBLBEL(:,1)=TBLBEL(:,1)-1;

TBLISR=[DATEISR';(TblISR.people_vaccinated)';(TblISR.people_fully_vaccinated)';]';
TBLISR=[0 0 0 ;TBLISR];

TBLISR(isnan(TBLISR(:,3)),3)=0;

%define the population size
%Nsize=[Belarusia,Israel];
Nsize=[9449323,9053000];

figure(1);
subplot(1,2,1);
% h2(5)=plot(TBLBEL(:,1),100*TBLBEL(:,2)/Nsize(1),'*','MarkerSize',3,'color',col_p(1,:));hold on;
% h2(6)=plot(TBLISR(:,1),100*TBLISR(:,2)/Nsize(2),'*','MarkerSize',3,'color',col_p(2,:));hold on;
h2(3)=plot(TBLBEL(:,1),100*TBLBEL(:,2)/Nsize(1),'*','MarkerSize',3,'color',col_p(1,:));hold on;
h2(4)=plot(TBLISR(:,1),100*TBLISR(:,2)/Nsize(2),'*','MarkerSize',3,'color',col_p(2,:));hold on;
for upsilon=[minupsilon maxupsilon]
   figure(figc); 
   subplot(1,2,1)
   h1(counter)=plot(tarr,100*(1-exp(-upsilon*tarr)),'color',col_p(counter,:),'LineWidth',4); hold on;
   mucounter=1;
   for mu1=mu1arr(2)
       init1=[(1-PerCompl)*N,PerCompl*N,0];
       pars1=[delta,mu0,mu1,upsilon,Incid];
       [t1,y1]=ode45(@(t1,y1)VaccComplReduced(t1,y1,pars1),[0,T], init1,opts);
       switch mu1
           case mu1arr(1)
               subplot(1,2,2);h2(2*(counter-1)+mucounter)=plot(t1,100*y1(:,2)/N,'LineWidth',4,'color',col_p(counter,:),'MarkerEdgeColor',col_p(counter,:));hold on;
           otherwise
               %subplot(1,2,2);h2(2*(counter-1)+mucounter)=plot(t1,100*y1(:,2)/N,'-.','MarkerSize',8,'LineWidth',4,'color',col_p(counter,:),'MarkerEdgeColor',col_p(counter,:));hold on;
               subplot(1,2,2);h2(counter)=plot(t1,100*y1(:,2)/N,'-.','MarkerSize',8,'LineWidth',4,'color',col_p(counter,:),'MarkerEdgeColor',col_p(counter,:));hold on;
        end
        mucounter=mucounter+1;
   end
   counter=counter+1;
end

subplot(1,2,1)
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
xlim([0,200]);
ylim([0,100]);

xlabel('Time (days)','interpreter','latex');
ylabel('Vaccination coverage (\%)','interpreter','latex');
%legend(h1,'Slow vaccination','Fast vaccination','Location','southoutside','NumColumns',2);
set(gca,'FontSize',25);

subplot(1,2,2)
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
xlim([0,200]);
ylim([0,100]);

xlabel('Time (days)','interpreter','latex');
ylabel('Compliant population (\%)','interpreter','latex');
%legend(h2,'Slow vaccination, slow compliance waning','Slow vaccination, fast compliance waning','Fast vaccination, slow compliance waning','Fast vaccination, fast compliance waning', 'Belarus data', 'Israel data','Location','southoutside','NumColumns',3);
legend(h2,'Slow vaccination','Fast vaccination','Belarus data', 'Israel data','Location','southoutside','NumColumns',2);
set(gca,'FontSize',25);

figure(figc);annotation('textbox', [0.01, 0.999, 0, 0], 'string', 'a','FontWeight','bold','FontSize',35)
figure(figc);annotation('textbox', [0.51, 0.999, 0, 0], 'string', 'b','FontWeight','bold','FontSize',35)

