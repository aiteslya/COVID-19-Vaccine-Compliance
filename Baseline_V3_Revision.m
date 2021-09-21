%This is the script that plots time series for the original variant,
% a B.1.1.7-like variant with 1.5 infectiousness of the original variant,
% and delta-like variant with 2.25 infectiousness of the original variant,
% with and without the vaccination with high compliance rise and
%decay rates (intervention vaccinated have the same contact rate as
%non-compliant individuals
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 and for the new 3.75
% such that Re(0)=1.1 for the old strain and Re(0)=1.65 for the new strain
% the cumulative output is adjusted for the initial number of infections
%prepare settings
clc;
clear variables;
close all;
format long;

%set the parameters

%define colors
l=5;
red = [0.01, 0.4, 0.76];
pink = [0.45,0.76,0.98];
medblue=(red+pink)/2;
novacc=[0.9,0.13,0.13];
yearcol=[0.17,0.17,0.17];
colors_p = [linspace(pink(1),red(1),l)', linspace(pink(2),red(2),l)', linspace(pink(3),red(3),l)'];
col_p=colors_p([1,5],:);

%contact rate of non-compliant before the lockdown at the start of the
%epidemic
chat=14.9;
%1/gamma duration of infectious period
gamma=1/7;
alpha=1/4;
R0arr=[2.5, 2.5*1.5 4.92];
%calculate epsilon
epsilonarr=R0arr*gamma/chat;

%set up initial data for the original variant and B.1.1.7-like variants
%seroprevalence
SP=0.08;
%total population
N=1.7e7;
TotalInf=112435;
Incid=TotalInf*gamma;
mu0=1/30;
delta=4e-5;
%percentage of compliant people
PerCompl=0.65;
N0=N*(1-PerCompl);
Nc0=N*PerCompl;

exposed=Incid/alpha;
TotalRec=SP*N;
TotalS=N-TotalInf-TotalRec-exposed;
%setting up of initial data
S0=(1-PerCompl)*TotalS;
E0=(1-PerCompl)*exposed;
I0=(1-PerCompl)*TotalInf;
R0=(1-PerCompl)*TotalRec;
Sc0=PerCompl*TotalS;
Ec0=PerCompl*exposed;
Ic0=PerCompl*TotalInf;
Rc0=PerCompl*TotalRec;
V0=0;
SV0=0;
EV0=0;
IV0=0;
RV0=0;
TV0=SV0+EV0+IV0+RV0;
%compliance

infect0=E0+I0+R0+Ec0+Ic0+Rc0;
popOut=1e5;
infect0pop=popOut*infect0/N;

V0=0;
SV0=0;
EV0=0;
IV0=0;
RV0=0;

TV0=SV0+EV0+IV0+RV0;
init=[S0,E0,I0,R0,Sc0,Ec0,Ic0,Rc0,V0,SV0,EV0,IV0,RV0,TV0,0,0,0];

%define array of contact rates
r1num=60;
r1arr=linspace(0,1,r1num);
%define the array of r1
num=2e4;
carr=linspace(0,15,num);

%for each c in carr calculate r1
%set up equation Re=1

Res=nan(r1num,2);
counter=1;
for r1=r1arr
    beta=carr*epsilonarr(1);
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
%average contact rate
ACR0=5;

ACR=Res(:,1)*(1-PerCompl)+Res(:,2).*Res(:,1)*PerCompl-ACR0;
eqn=ACR(1:r1num-1).*ACR(2:r1num);
ind=find(eqn<0,1);
c=Res(ind,1);
r1=Res(ind,2);
%integration options
Atol=1e-11;
opts = odeset('RelTol',1e-12,'AbsTol',Atol);
%integrating time: slightly more than 6 months
T=200;

numout=2;
minupsilon=5.9e-4;
maxupsilon=0.0053;%israel vaccination rate
upsilonarr=[0,minupsilon,maxupsilon];
k1=1;
k2=1;
omega=0.6;
%format of the legend
formatSpec = '%.2e';

r2=13.47/c;
frac=1/3;
mu1=3./(frac*5.1e8);
mu1hat=mu1*N;

%legend for the vaccinated plot
counter=1;
numpoints=30;
tdiscr=linspace(0,T,numpoints);

XpopOut=N/popOut;
%the following the counter for the type of variant
i1=1;
%container for vaccination cumulative infected at mark of 3 months and 6
%months
CumContainer=zeros(4,2);
for epsilon=epsilonarr(1:2)
    beta=epsilon*c;
    i2=1;
    for upsilon0=upsilonarr
        pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
        [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
                
        infectious=(y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16))/XpopOut;
        cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))/XpopOut-infect0pop;
        %new
        cumvacc=(y(:,15)+y(:,16)+y(:,17))/XpopOut;
        compl=(y(:,5)+y(:,6)+y(:,7)+y(:,8))/XpopOut;
        
        if upsilon0==0 %plot 
            %%%
            if i1==1
                figure(6);
                subplot(1,3,1);
                plot(t,infectious,'-.','color',novacc,'LineWidth',4);hold on;
                
            else
                figure(6);
                subplot(1,3,2);
                plot(t,infectious,'color',novacc,'LineWidth',4);hold on;
                
            end
            %collect baseline cumulative infections at mark of 1 year and 2
            %years
            %find index of the first year
            ind1year=find(t>=90,1);
            baseCum(i1,1)=cum(ind1year);
            %find index of the second year
            ind2year=find(t>=180,1);
            baseCum(i1,2)=cum(ind2year);
            disp('Baseline cumulative infections have been collected');
        elseif upsilon0==upsilonarr(2)
            if i1==1
                figure(6);
                subplot(1,3,1);
                plot(t,infectious,'LineWidth',4,'color',col_p(1,:));hold on;
                                
                figure(8);% v1, slow
                subplot(1,3,1);
                h8(2)=area(t,100*cum/popOut,'facecolor',novacc,'edgecolor','none');hold on;
                h8(3)=area(t,100*cumvacc/popOut,'facecolor',medblue,'edgecolor','none');hold on;
                xlim([0,T]);
                %ylim([0,40]);
            else
                figure(6);
                subplot(1,3,2);
                plot(t,infectious,'LineWidth',4,'color',col_p(1,:));hold on;
                                
                figure(8);% v2, slow
                subplot(1,3,2);
                h9(2)=area(t,100*cum/popOut,'facecolor',novacc,'edgecolor','none');hold on;
                h9(3)=area(t,100*cumvacc/popOut,'facecolor',medblue,'edgecolor','none');hold on;
                xlim([0,T]);
                %ylim([0,40]);
            end
        else
            if i1==1
                figure(6);
                subplot(1,3,1);
                plot(t,infectious,'LineWidth',4,'color',col_p(2,:));hold on;
                                
                figure(9);% v1, fast
                subplot(1,3,1);
                h8(4)=area(t,100*cum/popOut,'facecolor',novacc,'edgecolor','none');hold on;
                h8(5)=area(t,100*cumvacc/popOut,'facecolor',medblue,'edgecolor','none');hold on;
                xlim([0,T]);
               % ylim([0,40]);
            else
                figure(6);
                subplot(1,3,2);
                plot(t,infectious,'LineWidth',4,'color',col_p(2,:));hold on;
                  
                figure(9);% v2, fast
                subplot(1,3,2);
                h9(4)=area(t,100*cum/popOut,'facecolor',novacc,'edgecolor','none');hold on;
                h9(5)=area(t,100*cumvacc/popOut,'facecolor',medblue,'edgecolor','none');hold on;
                xlim([0,T]);
                %ylim([0,40]);
            end
        end
        if upsilon0>0
            ind1year=find(t>=90,1);
            CumContainer(2*i1-2+i2,1)=cum(ind1year);
            ind2year=find(t>=180,1);
            CumContainer(2*i1-2+i2,2)=cum(ind2year);
            i2=i2+1;
        end
        
    end
    i1=i1+1;
end
%figure settings

figure(6);
subplot(1,3,1);
xlim([0, T]);
ylim([0,2300]);

xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
ylabel('Infected individuals (1/100,000)','interpreter','latex');
title('Original variant');
legend('Baseline, no vaccination','Slow vaccination','Fast vaccination','location','southoutside','NumColumns',3);
set(gca,'FontSize',25);

%figure settings
subplot(1,3,2);
xlim([0, T]);
ylim([0,2300]);
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
title('Alpha-like variant');
set(gca,'FontSize',25);

%create bar charts for old and new variants cumulative infections
CumContainer(1,1)=100*(CumContainer(1,1)-baseCum(1,1))/baseCum(1,1);
CumContainer(1,2)=100*(CumContainer(1,2)-baseCum(1,2))/baseCum(1,2);
CumContainer(2,1)=100*(CumContainer(2,1)-baseCum(1,1))/baseCum(1,1);
CumContainer(2,2)=100*(CumContainer(2,2)-baseCum(1,2))/baseCum(1,2);

temp=CumContainer(1,2);
CumContainer(1,2)=CumContainer(2,1);
CumContainer(2,1)=temp;

CumContainer(3,1)=100*(CumContainer(3,1)-baseCum(2,1))/baseCum(2,1);
CumContainer(3,2)=100*(CumContainer(3,2)-baseCum(2,2))/baseCum(2,2);
CumContainer(4,1)=100*(CumContainer(4,1)-baseCum(2,1))/baseCum(2,1);
CumContainer(4,2)=100*(CumContainer(4,2)-baseCum(2,2))/baseCum(2,2);

temp=CumContainer(3,2);
CumContainer(3,2)=CumContainer(4,1);
CumContainer(4,1)=temp;

figure(5);
subplot(1,3,1);
p1=bar(CumContainer(1:2,:));
ylim([0,150])
set(p1(1),'FaceColor',col_p(1,:));%[116 120 128]/255);
set(p1(2),'FaceColor',col_p(2,:));%[24 74 69]/255);

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'3 months' '6 months'})
title('Original variant');
ylabel('Excess infections (\%)','interpreter','latex');
%ylim([-5,70]);
legend('Slow vaccination','Fast vaccination','Location','southoutside','NumColumns',2,'interpreter','latex');
set(gca,'FontSize',25);

subplot(1,3,2);
p2=bar(CumContainer(3:4,:));
ylim([0,150])
set(p2(1),'FaceColor',col_p(1,:));
set(p2(2),'FaceColor',col_p(2,:));
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'3 months' '6 months'})
title('Alpha-like variant');
%ylim([-5,70]);
set(gca,'FontSize',25);

%output settings for figures 8 and 9
f8=figure(8);
set(gca,'FontSize',25);
subplot(1,3,1);
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
ylim([0,30]);
%ylim([0,6000]);
xlabel('Time (days)','interpreter','latex');
%ylabel({'Cumulative number of infected';'individuals (1/100,000)'},'interpreter','latex');
ylabel('Attack rate ($$\%$$)','interpreter','latex');
title('Original variant','FontWeight','normal');
%legend(h8,'Baseline, no vaccination    ','Total, with vaccination     ','Vaccinated, with vaccination','Total, with vaccination     ','Vaccinated, with vaccination','NumColumns',2,'location','southoutside');
legend(h8(4:5),'Non-vaccinated','Vaccinated','NumColumns',2,'location','southoutside');
set(gca,'FontSize',25);
subplot(1,3,2);
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
ylim([0,30]);
%ylim([0,100]);
xlabel('Time (days)','interpreter','latex');
title('Alpha-like variant','FontWeight','normal');


% han2.XLabel.Visible='on';
% han2.YLabel.Visible='on';

figure(9);
subplot(1,3,1);
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
ylim([0,30]);
title('Original variant','FontWeight','normal')
%ylim([0,2.3e4]);
xlabel('Time (days)','interpreter','latex');
%ylabel({'Cumulative number of infected';'individuals (1/100,000)'},'interpreter','latex');
ylabel('Attack rate ($$\%$$)','interpreter','latex');
set(gca,'FontSize',25);
subplot(1,3,2);
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
ylim([0,30]);
xlabel('Time (days)','interpreter','latex');
title('Alpha-like variant','FontWeight','normal');
legend(h9(4:5),'Non-vaccinated','Vaccinated','NumColumns',2,'location','southoutside');
set(gca,'FontSize',25);

%output the figures for delta-like variant
R0delta=4.92;
epsilonDelta=R0delta*gamma/chat;
beta=epsilonDelta*c;
i2=1;
%initial conditions
%find the time when delta became 50% of all infections 23 June 2021 delta
%started to be the dominant strain
%find the compliant proportion on June 29
% PerCompl=0.63;
% %find the infectious proporion on 23 of June
% TotalInf=175466;
% %calculate incidence - these also include vaccinated
% Incid=TotalInf*gamma;
% %calculate exposed;
% exposed=Incid/alpha;
% %calculate recovered
% SP=0.2;
% recovered=N*SP;
% %find the coverage of vaccinated with one or two doses
% vaccinated=2659985;
% vaccprop=vaccinated/N;
% perNonCompl=1-PerCompl-vaccprop;
% vaccnonimmunprop=vaccprop*(1-omega);
% vaccimune=omega*vaccinated;
% %these also include infected, exposed and recovered!!
% vaccnonimmune=(1-omega)*vaccinated;
% %all susceptibles, including vaccinated
% susceptible=N-TotalInf-exposed-recovered-vaccimune;
% snoncompl=perNonCompl*susceptible/(perNonCompl+PerCompl+vaccnonimmunprop);
% scompl=PerCompl*susceptible/(perNonCompl+PerCompl+vaccnonimmunprop);
% 
% enoncompl=perNonCompl*exposed/(perNonCompl+PerCompl+vaccnonimmunprop);
% ecompl=PerCompl*exposed/(perNonCompl+PerCompl+vaccnonimmunprop);
% evacc=vaccnonimmunprop*exposed/(perNonCompl+PerCompl+vaccnonimmunprop);
% 
% inoncompl=perNonCompl*TotalInf/(perNonCompl+PerCompl+vaccnonimmunprop);
% icompl=PerCompl*TotalInf/(perNonCompl+PerCompl+vaccnonimmunprop);
% ivacc=vaccnonimmunprop*TotalInf/(perNonCompl+PerCompl+vaccnonimmunprop);
% 
% rnoncompl=perNonCompl*recovered/(perNonCompl+PerCompl+vaccnonimmunprop);
% rcompl=PerCompl*recovered/(perNonCompl+PerCompl+vaccnonimmunprop);
% rvacc=vaccnonimmunprop*recovered/(perNonCompl+PerCompl+vaccnonimmunprop);
% 
% svacc=vaccinated-vaccimune-evacc-ivacc-rvacc;

% set initial conditions
%init=[snoncompl,enoncompl,inoncompl,rnoncompl,scompl,ecompl,icompl,rcompl,vaccimune,svacc,evacc,ivacc,rvacc,vaccinated,0,0,0];
%infect0=enoncompl+inoncompl+rnoncompl+ecompl+icompl+rcompl+evacc+ivacc+rvacc;
%infect0pop=infect0/XpopOut;
%set the total
i1=1;
i2=1;
baseCum=zeros(1,2);
CumContainer=zeros(2,2);
r2=13.47/c;
fracarr=2/3;
mu1=3./(fracarr*5.1e8);
for upsilon0=upsilonarr
    pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
    [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);

    prev=(y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16))/XpopOut;
    cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))/XpopOut-infect0pop;
    %new
    cumvacc=(y(:,15)+y(:,16)+y(:,17))/XpopOut;
    compl=(y(:,5)+y(:,6)+y(:,7)+y(:,8))/XpopOut;

    if upsilon0==0 %plot 
        figure(6);
        subplot(1,3,3);
        plot(t,prev,'-.','color',novacc,'LineWidth',4);hold on;
        %collect baseline cumulative infections at mark of 1 year and 2
        %years
        %find index of the first year
        ind1year=find(t>=90,1);
        baseCum(i1,1)=cum(ind1year);
        %find index of the second year
        ind2year=find(t>=180,1);
        baseCum(i1,2)=cum(ind2year);
        figure(10);plot(t,compl,'-.','color',novacc,'LineWidth',4);hold on;
        figure(11);plot(t,cum,'-.','color',novacc,'LineWidth',4);hold on;
    elseif upsilon0==upsilonarr(2)
        figure(6);
        subplot(1,3,3);
        plot(t,prev,'LineWidth',4,'color',col_p(1,:));hold on;
        figure(10);plot(t,compl,'LineWidth',4,'color',col_p(1,:));hold on;
        figure(11);plot(t,cum,'LineWidth',4,'color',col_p(1,:));hold on;
        figure(8);% v1, slow
        subplot(1,3,3);
        h8(2)=area(t,100*cum/popOut,'facecolor',novacc,'edgecolor','none');hold on;
        h8(3)=area(t,100*cumvacc/popOut,'facecolor',medblue,'edgecolor','none');hold on;
        xlim([0,T]);
        
    else
        figure(6);
        subplot(1,3,3);
        plot(t,prev,'LineWidth',4,'color',col_p(2,:));hold on;
        figure(10);plot(t,compl,'LineWidth',4,'color',col_p(2,:));hold on;
        figure(11);plot(t,cum,'LineWidth',4,'color',col_p(2,:));hold on;
        figure(9);% v1, fast
        subplot(1,3,3);
        h8(4)=area(t,100*cum/popOut,'facecolor',novacc,'edgecolor','none');hold on;
        h8(5)=area(t,100*cumvacc/popOut,'facecolor',medblue,'edgecolor','none');hold on;
        xlim([0,T]);
    end
    if upsilon0>0
        ind1year=find(t>=90,1);
        CumContainer(2*i1-2+i2,1)=cum(ind1year);
        ind2year=find(t>=180,1);
        CumContainer(2*i1-2+i2,2)=cum(ind2year);
        i2=i2+1;
    end
        
end

%figure settings

figure(6);
subplot(1,3,3);
xlim([0, T]);
ylim([0,2300]);

xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);

xlabel('Time (days)','interpreter','latex');
title('Delta-like variant');
set(gca,'FontSize',25);

%create bar charts for old and new variants cumulative infections
CumContainer(1,1)=100*(CumContainer(1,1)-baseCum(1,1))/baseCum(1,1);
CumContainer(1,2)=100*(CumContainer(1,2)-baseCum(1,2))/baseCum(1,2);
CumContainer(2,1)=100*(CumContainer(2,1)-baseCum(1,1))/baseCum(1,1);
CumContainer(2,2)=100*(CumContainer(2,2)-baseCum(1,2))/baseCum(1,2);

temp=CumContainer(1,2);
CumContainer(1,2)=CumContainer(2,1);
CumContainer(2,1)=temp;

figure(5);
subplot(1,3,3);
p1=bar(CumContainer(1:2,:));
ylim([0,150])
set(p1(1),'FaceColor',col_p(1,:));%[116 120 128]/255);
set(p1(2),'FaceColor',col_p(2,:));%[24 74 69]/255);

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'3 months' '6 months'})
title('Delta-like variant');

set(gca,'FontSize',25);

%output settings for figures 8 and 9
figure(8);
set(gca,'FontSize',25);
subplot(1,3,3);
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
ylim([0,30]);
xlabel('Time (days)','interpreter','latex');

title('Delta-like variant','FontWeight','normal');
han8=axes(f8,'visible','off');

han8.Title.Visible='on';
tit=title(han8,{'Slow vaccination'},'interpreter','latex');
set(tit,'FontWeight','bold');
set(tit,'FontSize',30);

f9=figure(9);
set(gca,'FontSize',25);
subplot(1,3,3);
xline(90,'-.','color',yearcol,'LineWidth',2);
xline(180,'-.','color',yearcol,'LineWidth',2);
ylim([0,30]);
title('Delta-like variant','FontWeight','normal');
%ylim([0,2.3e4]);
xlabel('Time (days)','interpreter','latex');

han9=axes(f9,'visible','off');

han9.Title.Visible='on';
tit=title(han9,{'Fast vaccination'},'interpreter','latex');
set(tit,'FontWeight','bold');
set(tit,'FontSize',30);

%set(gca,'FontSize',25);

%labeling of the panels
figure(6);annotation('textbox', [0.05, 0.99, 0, 0], 'string', 'a','FontWeight','bold','FontSize',30)
figure(6);annotation('textbox', [0.35, 0.99, 0, 0], 'string', 'b','FontWeight','bold','FontSize',30)
figure(6);annotation('textbox', [0.66, 0.99, 0, 0], 'string', 'c','FontWeight','bold','FontSize',30)

figure(5);annotation('textbox', [0.05, 0.99, 0, 0], 'string', 'd','FontWeight','bold','FontSize',30)
figure(5);annotation('textbox', [0.35, 0.99, 0, 0], 'string', 'e','FontWeight','bold','FontSize',30)
figure(5);annotation('textbox', [0.66, 0.99, 0, 0], 'string', 'f','FontWeight','bold','FontSize',30)

figure(8);annotation('textbox', [0.05, 0.96, 0, 0], 'string', 'a','FontWeight','bold','FontSize',30)
figure(8);annotation('textbox', [0.35, 0.96, 0, 0], 'string', 'b','FontWeight','bold','FontSize',30)
figure(8);annotation('textbox', [0.66, 0.96, 0, 0], 'string', 'c','FontWeight','bold','FontSize',30)

figure(9);annotation('textbox', [0.05, 0.96, 0, 0], 'string', 'd','FontWeight','bold','FontSize',30)
figure(9);annotation('textbox', [0.35, 0.96, 0, 0], 'string', 'e','FontWeight','bold','FontSize',30)
figure(9);annotation('textbox', [0.66, 0.96, 0, 0], 'string', 'f','FontWeight','bold','FontSize',30)