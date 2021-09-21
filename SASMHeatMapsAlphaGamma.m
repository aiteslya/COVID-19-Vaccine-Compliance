 %This is the script that plots heat maps for the original variant,
% vaccine efficacy of 60% and fast compliance waning as the vaccination 
% coverage grows and significant foregoing of compliance with variation in 
% the durations of the latent and infectious periods keeping the rest of the parameters fixed
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

%set up initial data
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
%define colors for the vaccination rates

figc=1;
counter=1;
tarr=linspace(0,200,100);
% numpoints=100;
% tdiscr=linspace(0,T,numpoints);
omega=0.6;
%define alpha and gamma arrays
alphaarr=linspace(1/6,1/2,meshsize);
gammaarr=linspace(1/9,1/5,meshsize);
[Aarr,Garr]=meshgrid(alphaarr,gammaarr);

%%
% original variant
beta=epsilon*c;

%outer loop: loops through the vaccination uptake rates
%figure counter
subfigc=1;
for upsilon=upsilonarr
    
    %allocate result arrays
    PeakRes=zeros(meshsize,meshsize);
    PeakTime=zeros(meshsize,meshsize);
    Cum3=zeros(meshsize,meshsize);
    Cum6=zeros(meshsize,meshsize);
    Cum03=zeros(meshsize,meshsize);
    Cum06=zeros(meshsize,meshsize);
    for i1=1:meshsize
        for i2=1:meshsize 
            alpha=Aarr(i1,i2);
            gamma=Garr(i1,i2);
            %collect outputs in the no-vaccination scenario
            upsilontemp=0;
            pars=[beta,r1,r2,delta,mu0,mu1,upsilontemp,alpha,gamma,k1,k2,omega];
            [t0,y0]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
            cum0=(y0(:,2)+y0(:,3)+y0(:,4)+y0(:,6)+y0(:,7)+y0(:,8)+y0(:,11)+y0(:,12)+y0(:,13)+y0(:,15)+y0(:,16)+y0(:,17))-infect0;
            %find indices of 3 and 6 months
            ind03=find(t0>3*30,1);
            ind06=find(t0>6*30,1);
            %cumulative infected
            cum03=100*cum0(ind03)/N;
            cum06=100*cum0(ind06)/N;
            Cum03(i1,i2)=cum03;
            Cum06(i1,i2)=cum06;
            
            pars=[beta,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega];
            [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
            % collect and process the output
            %prevalence
            prev=y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16);
            cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;
            maxi1i2=max(prev);
            indmax=find(prev==maxi1i2,1);
            tmax=t(indmax);
            %find indices of 3 and 6 months
            ind3=find(t>3*30,1);
            ind6=find(t>6*30,1);
            %cumulative infected
            cum3=100*cum(ind3)/N;
            cum6=100*cum(ind6)/N;
            PeakRes(i1,i2)=100*maxi1i2/N;
            PeakTime(i1,i2)=tmax;
            Cum3(i1,i2)=cum3;
            Cum6(i1,i2)=cum6;
        end
    end

    f1=figure(1);
    subplot(2,2,subfigc);
    
    
    s=surf(1./Aarr,1./Garr,Cum3);
    ylabel('Infectious period (days)','interpreter','latex');
    colormap jet

    s.EdgeColor = 'none';
    hold on;
    view(2);
    [M,C]=contour(1./Aarr,1./Garr,Cum3,[0,0]);
    C.LineColor='m';
    C.LineWidth=4;

    caxis([0.5,7])
    
    if upsilon==upsilonarr(1) 
        title('Slow vaccination','interpreter','latex');
    else
        title('Fast vaccination','interpreter','latex');
    end
    set(gca,'FontSize',20);

    f1=figure(1);
    subplot(2,2,2+subfigc);
    s=surf(1./Aarr,1./Garr,Cum6);
    ylabel('Infectious period (days)','interpreter','latex');
    ylim([5,9]);
    colormap jet
    colorbar
    
    xlabel('Exposed period (days)','interpreter','latex');
    
    caxis([0.5,7])
    s.EdgeColor = 'none';
    hold on;
    view(2);
    
    if upsilon==upsilonarr(1) 
        title('Slow vaccination','interpreter','latex');
    else
        title('Fast vaccination','interpreter','latex');
    end
    set(gca,'FontSize',20);
    
    %output relative differences
    
    f2=figure(2);
    subplot(2,2,subfigc);
    
    
    s=surf(1./Aarr,1./Garr,100*(Cum3-Cum03)./Cum03);
    ylabel('Infectious period (days)','interpreter','latex');
    colormap jet

    s.EdgeColor = 'none';
    hold on;
    view(2);
    [M,C]=contour(1./Aarr,1./Garr,100*(Cum3-Cum03)./Cum03,[0,0]);
    C.LineColor='m';
    C.LineWidth=4;

    caxis([2,100])
    
    if upsilon==upsilonarr(1) 
        title('Slow vaccination','interpreter','latex');
    else
        title('Fast vaccination','interpreter','latex');
    end
    set(gca,'FontSize',20);

    f1=figure(2);
    subplot(2,2,2+subfigc);
    s=surf(1./Aarr,1./Garr,100*(Cum6-Cum06)./Cum06);
    ylabel('Infectious period (days)','interpreter','latex');
    ylim([5,9]);
    colormap jet
    if subfigc==2
        colorbar
    end
    
    s.EdgeColor = 'none';
    hold on;
    view(2);
    [M,C]=contour(1./Aarr,1./Garr,100*(Cum6-Cum06)./Cum06,[0,0]);
    C.LineColor='m';
    C.LineWidth=4;
    
    caxis([2,100])
    
    xlabel('Exposed period (days)','interpreter','latex');
    
    %caxis([0.5,7])
    s.EdgeColor = 'none';
    hold on;
    view(2);
    
    if upsilon==upsilonarr(1) 
        title('Slow vaccination','interpreter','latex');
    else
        title('Fast vaccination','interpreter','latex');
    end
    set(gca,'FontSize',20);
        
    subfigc=subfigc+1;
end

tit=suptitle({'Cumulative number of',' new infections at three months (\%)'});
set(tit,'FontSize',20);

%annotations
figure(1);annotation('textbox', [0.08, 0.99, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.53, 0.99, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)

figure(1);annotation('textbox', [0.08, 0.49, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)
figure(1);annotation('textbox', [0.53, 0.49, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)

tit=suptitle({'Relative difference in the number';'of new infections at three months (\%)'});
set(tit,'FontSize',20);

%annotations
figure(2);annotation('textbox', [0.08, 0.99, 0, 0], 'string', 'a','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.53, 0.99, 0, 0], 'string', 'b','FontWeight','bold','FontSize',25)

figure(2);annotation('textbox', [0.08, 0.49, 0, 0], 'string', 'c','FontWeight','bold','FontSize',25)
figure(2);annotation('textbox', [0.53, 0.49, 0, 0], 'string', 'd','FontWeight','bold','FontSize',25)
