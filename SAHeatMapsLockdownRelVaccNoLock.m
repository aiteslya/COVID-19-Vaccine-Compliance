%This is the script that plots heat maps for the difference in the
% cumulative number of new infections in the vaccine uptake - with different
% threshold for commencement of the lockdown and the strictness of the lockdown
% keeping the rest of the parameters fixed
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 with Re(0)=1.1 at the start of the
% vaccination rollout
% the outputs are the relative difference in the cumulative
% number of new infections relative to the vaccination scenario without lockdown
% adjusted to not count the infected and recovered at time 0
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

figure(20);
plot(Res(:,1).*Res(:,2),Res(:,1),'o-','LineWidth',2);
xlabel({'Contact rate of compliant individuals';'(individuals/day)'},'interpreter','latex');
ylabel({'Contact rate of non-compliant';' individuals (individuals/day)'},'interpreter','latex');
xlim([0.9*min(Res(:,1).*Res(:,2)),1.01*max(Res(:,1).*Res(:,2))]);
set(gca,'FontSize',25);
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
r2=13.47/c;
%fraction of population that has been vaccinated when the compliance state
%on average lasts 7 days
fracarr=1/3;
mu1=3./(fracarr*5.1e8);
%set up lockdown options
cLockarr=[3,5];
thresharr=[50,500]./(1.7e7);
% set up main arrays for the heat maps
meshsize=10;
minupsilon=5e-4;%5.9e-4;
maxupsilon=6e-3;%0.0053;% vaccination rate in Israel
upsilonarr=linspace(minupsilon,maxupsilon,meshsize);%[5.9e-4,4.9e-3];
% output the figure with minimum and maximum vaccination curves
%define colors for the vaccination rates

figc=1;
counter=1;
tarr=linspace(0,200,100);
numpoints=100;
tdiscr=linspace(0,T,numpoints);
omegaarr=linspace(0.55,0.95,meshsize);
[Uarr,Oarr]=meshgrid(upsilonarr,omegaarr);

%%
% original variant
beta=epsilon*c;

%outer loop: circles through cLock: the average contact rate in the
%lockdown
%figure counter
subfigc=1;

% preset xtics and their labels
num_tics=7;
xtpoints=linspace(minupsilon,maxupsilon,num_tics);
xtpoints=xtpoints(2:num_tics-1);
%numeric values: vaccination coverage 3 and 6 months after the vaccination
%rollout

%xtCoverage3=100*(1-exp(-xtpoints*90));

xtCoverage3=10:5:40;
xtpoints3=(1/(90))*log(100./(100-xtCoverage3))/log(exp(1));
for count=1:numel(xtCoverage3)
   xtlabelsStr3{count}=num2str(round(xtCoverage3(count))); 
end

xtCoverage6=20:5:70;%100*(1-exp(-xtpoints*180));
xtpoints6=(1/(180))*log(100./(100-xtCoverage6))/log(exp(1));
for count=1:numel(xtCoverage6)
   xtlabelsStr6{count}=num2str(round(xtCoverage6(count))); 
end

clockD=3;
%the next loop: circles through mu1: sensitivity of the compliance
%waning to vaccination coverage
for pc=thresharr
    %allocate result arrays
    PeakRes=zeros(meshsize,meshsize);
    Cum3=zeros(meshsize,meshsize);
    Cum6=zeros(meshsize,meshsize);
    for i1=1:meshsize
        for i2=1:meshsize 
                       
            upsilon=Uarr(i1,i2);
            omega=Oarr(i1,i2);
            
            % generate no lockdown vaccination statistics statistics
            pars=[beta,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega];
            [t0,y0]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
            % collect outputs
            % cumulatively infected
            cum0=(y0(:,2)+y0(:,3)+y0(:,4)+y0(:,6)+y0(:,7)+y0(:,8)+y0(:,11)+y0(:,12)+y0(:,13)+y0(:,15)+y0(:,16)+y0(:,17))-infect0;

            %find indices of 3,6
            ind03=find(t0>3*30,1);
            ind06=find(t0>6*30,1);
            %cumulative infected
            cum03=cum0(ind03);
            cum06=cum0(ind06);

            %generate lockdown statistics
            pars=[clockD,c,epsilon,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega,pc];
            [t,y]=ode23s(@(t,y)COVIDVaccIntervenRHS(t,y,pars),[0,T], init,opts);
            % collect and process the output
            %prevalence
            prev=y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16);
            cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;
            maxi1i2=max(prev);
            %find indices of 3,6,12,24 months
            ind3=find(t>3*30,1);
            ind6=find(t>6*30,1);
            %cumulative infected
            cum3=cum(ind3);
            cum6=cum(ind6);

            PeakRes(i1,i2)=maxi1i2;
            Cum3(i1,i2)=100*(cum3-cum03)/cum03;
            Cum6(i1,i2)=100*(cum6-cum06)/cum06;
        end
    end

    f1=figure(figc);
    subplot(2,2,2*subfigc-1);
    s=surf(Uarr,100*Oarr,Cum3);
    colormap jet
    caxis([-75,-40]);

    f=subplot(2,2,2*subfigc-1);
    hp4 = get(subplot(2,2,2*subfigc-1),'Position')
    f.Position=[hp4(1) hp4(2) hp4(3),0.2];

    s.EdgeColor = 'none';
    hold on;
    view(2);
    [M,C]=contour(Uarr,100*Oarr,Cum3,[0,0]);
    C.LineColor='m';
    C.LineWidth=4;
    % colorbar will need to be updated, depending on the row of the
    % figure
    
    set(gca,'xtick',xtpoints3)
    set(gca,'xticklabel',xtlabelsStr3)

    %caxis([-65,-40]);
    %xlabel({'Vaccination coverage after six months (\%)'},'interpreter','latex');
    title('Three months','interpreter','latex');
    if pc==thresharr(1) & subfigc==1
        ylabel({'\textbf{Lockfown threshold}';'\textbf{50 individuals}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
    end

    set(gca,'FontSize',15);

    subplot(2,2,2*subfigc);
    s=surf(Uarr,100*Oarr,Cum6);
    colormap jet
    % colorbar will need to be updated, depending on the row of the
    caxis([-75,-40]);
    s.EdgeColor = 'none';
    hold on;
    view(2);
    [M,C]=contour(Uarr,100*Oarr,Cum6,[0,0]);
    C.LineColor='m';
    C.LineWidth=4;

    f=subplot(2,2,2*subfigc);
    hp4 = get(subplot(2,2,2*subfigc),'Position')
    f.Position=[f.Position(1) f.Position(2) f.Position(3),0.2];

    %yline(60,'g','LineWidth',4);
    if subfigc==1
        xlabel({'Vaccination coverage after three months (\%)'},'interpreter','latex');
    else
        xlabel({'Vaccination coverage after six months (\%)'},'interpreter','latex');
    end
    
    if pc==thresharr(2) 
        ylabel({'\textbf{Lockdown threshold}';'\textbf{500 individuals}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
    end
    
    set(gca,'xtick',xtpoints6)
    set(gca,'xticklabel',xtlabelsStr6)
    title('Six months','interpreter','latex');

    set(gca,'FontSize',15);

    colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.04  hp4(2)*0.2+hp4(3)*0.1])

    subfigc=subfigc+1;
end

figc=figc+1;
