 %This is the script that plots heat maps varying upsilon,
%omega and mu1 and r2, keeping the rest of the parameters fixed
% during the lockdown for COVID-19-Compliance-Vaccine model
% R0 for the old virus strain is 2.5 with Re(0)=1.1 at the start of the
% vaccination rollout
% the outputs are the cumulative
% percentage of hospitalized relative to no vaccination scenario adjusted to
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

% compliance acquisition rate
delta=4e-5;
k1=1;
k2=1;

% hospitalization probability 2x3: rows unvaccinated/vaccinates, columns:
% original, alpha, delta
%vaccine efficacy in preventing hospitalization for original, alpha and
%delta variant
vacc_eff_arr=[0.95,0.85,0.85];

hosp_prob_arr(1,1)=0.04;%probability of becoming hospitalized when infected with original strain and non-vaccinated
hosp_prob_arr(1,2)=1.4*hosp_prob_arr(1,1);%0.11;%probability of becoming hospitalized when infected with alpha strain and non-vaccinated, slightly higher than for original
hosp_prob_arr(1,3)=1.85*hosp_prob_arr(1,2);%probability of becoming hospitalized when infected with delta strain and non-vaccinated, 1.85 higher than for alpha
hosp_prob_arr(2,:)=(1-vacc_eff_arr).*(hosp_prob_arr(1,:));


%format of the legend
formatSpec = '%.2e';

% set up arrays
% vaccinated can increase their contact rate as high as non-compliant
% during the lockdown or even go up to the levels of pre-pandemic
r2arr=[13.47/c,1];
%fraction of population that has been vaccinated when the compliance state
%on average lasts 7 days
fracarr=[1/3,2/3];
mu1arr=3./(fracarr*5.1e8);
% set up main arrays for the heat maps
meshsize=25;
minupsilon=5e-4;%5.9e-4;
maxupsilon=6e-3;%0.0053;% vaccination rate in Israel
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
%tarr=linspace(0,200,100);
r2=1;
%numpoints=100;
%tdiscr=linspace(0,T,numpoints);
%discretize the mesh for heat maps
omegaarr=linspace(0.55,0.95,meshsize);
[Uarr,Oarr]=meshgrid(upsilonarr,omegaarr);

%%
% original variant
beta=epsilon*c;

%generate no-vaccination statistics
upsilon=0;

% these parameters are not active in no vaccination scenario (upsilon=0)
mu1=mu1arr(1);
omega=1;

pars=[beta,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega];
[t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
% collect outputs
% cumulatively hospitalized
cum_hosp=hosp_prob_arr(1,1)*(y(:,3)+y(:,4)+y(:,7)+y(:,8)+y(:,12)+y(:,13)-(y(1,3)+y(1,4)+y(1,7)+y(1,8)+y(1,12)+y(1,13)));

%find indices of 3,6
ind3=find(t>3*30,1);
ind6=find(t>6*30,1);
%cumulative hospitalized following 3 and 6 months after the start of
%simulation
cum_hosp03=cum_hosp(ind3);
cum_hosp06=cum_hosp(ind6);

%outer loop: circles through r2: increase of the contact rate of the
%vaccinated
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

for r2=r2arr
    %the next loop: circles through mu1: sensitivity of the compliance
    %waning to vaccination coverage
    for mu1=mu1arr
        %allocate result arrays
        PeakRes=zeros(meshsize,meshsize);
        Cum3=zeros(meshsize,meshsize);
        Cum6=zeros(meshsize,meshsize);
        for i1=1:meshsize
            for i2=1:meshsize 
                upsilon=Uarr(i1,i2);
                omega=Oarr(i1,i2);
                           
                % in the next two lines the column index depends on the variant
                hosp_prob_vacc=hosp_prob_arr(2,1);
                hosp_prob_unvacc=hosp_prob_arr(1,1);
        
                pars=[beta,r1,r2,delta,mu0,mu1,upsilon,alpha,gamma,k1,k2,omega];
                [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
                % collect and process the output
                %prevalence
                prev=y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16);
                cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;

                %cumulative number of new hospitalizations since the start of the vaccination
                cum_hosp=hosp_prob_unvacc*(y(:,3)+y(:,4)+y(:,7)+y(:,8)+y(:,12)+y(:,13)-(y(1,3)+y(1,4)+y(1,7)+y(1,8)+y(1,12)+y(1,13)))+hosp_prob_vacc*(y(:,16)-y(1,16)+y(:,17)-y(1,17));

% cum=N-y(:,1)-y(:,5)-infect0;
% % cumulative hospitalized
% cum_hosp=hosp_prob_arr(1,1)*cum;

                maxi1i2=max(prev);
                %find indices of 3,6,12,24 months
                ind3=find(t>3*30,1);
                ind6=find(t>6*30,1);
                %cumulative infected
                cum3=cum_hosp(ind3);
                cum6=cum_hosp(ind6);
                
                PeakRes(i1,i2)=maxi1i2;
                Cum3(i1,i2)=100*(cum3-cum_hosp03)/cum_hosp03;
                Cum6(i1,i2)=100*(cum6-cum_hosp06)/cum_hosp06;
            end
        end
        
        f1=figure(figc);
        subplot(4,2,2*subfigc-1);
        s=surf(Uarr,100*Oarr,Cum3);
        colormap jet
        
        f=subplot(4,2,2*subfigc-1);
        hp4 = get(subplot(4,2,2*subfigc-1),'Position')
        f.Position=[hp4(1) hp4(2) hp4(3),0.2];
        
        s.EdgeColor = 'none';
        hold on;
        view(2);
        [M,C]=contour(Uarr,100*Oarr,Cum3,[0,0]);
        C.LineColor='m';
        C.LineWidth=4;
        % colorbar will need to be updated, depending on the row of the
        % figure
        
        caxis([-30,20]);
        if subfigc==4
            set(gca,'xtick',xtpoints3)
            set(gca,'xticklabel',xtlabelsStr3)
            xlabel({'Vaccination coverage after three months (\%)'},'interpreter','latex');
        else 
            set(gca,'xtick',xtpoints3)
            set(gca,'xticklabel',[])
        end
        
         if r2==r2arr(1) & mu1==mu1arr(1)
             ylabel({'\textbf{Fast compliance decay,}';'\textbf{high contact rate}';'\textbf{of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         elseif r2==r2arr(1) & mu1==mu1arr(2)
             ylabel({'\textbf{Slow compliance decay,}';'\textbf{high contact rate}';'\textbf{of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         elseif r2==r2arr(2) & mu1==mu1arr(1)
             ylabel({'\textbf{Fast compliance decay,}';'\textbf{low contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         else
             ylabel({'\textbf{Slow compliance decay,}';'\textbf{low contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         end
         if subfigc==1
            title('Three months since vaccination start','interpreter','latex','FontSize',23);
         end
         set(gca,'FontSize',12);
       
        subplot(4,2,2*subfigc);
        s=surf(Uarr,100*Oarr,Cum6);
        colormap jet
        % colorbar will need to be updated, depending on the row of the
        caxis([-30,20]);
        s.EdgeColor = 'none';
        hold on;
        view(2);
        [M,C]=contour(Uarr,100*Oarr,Cum6,[0,0]);
        C.LineColor='m';
        C.LineWidth=4;
        
        f=subplot(4,2,2*subfigc);
        hp4 = get(subplot(4,2,2*subfigc),'Position')
        f.Position=[f.Position(1) f.Position(2) f.Position(3),0.2];
        
        %yline(60,'g','LineWidth',4);
        if subfigc==4
            set(gca,'xtick',xtpoints6)
            set(gca,'xticklabel',xtlabelsStr6)
            xlabel({'Vaccination coverage after six months (\%)'},'interpreter','latex');
        else 
            set(gca,'xtick',xtpoints6)
            set(gca,'xticklabel',[])
        end
        if subfigc==1
            title('Six months since vaccination start','interpreter','latex','FontSize',23);
        end
        set(gca,'FontSize',12);
        
        colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.04  hp4(2)*0.2+hp4(3)*0.1])
       
        subfigc=subfigc+1;
    end
end

figc=figc+1;

%%
% Alpha-like variant
epsilonAlpha=1.5*epsilon;
beta=epsilonAlpha*c;

%generate no-vaccination statistics
upsilon0=0;

mu1=mu1arr(1);
omega=1;
r2=13.47/c;

pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
[t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
% collect outputs
cum_hosp=100*(hosp_prob_arr(1,2)*(y(:,3)+y(:,4)+y(:,7)+y(:,8)+y(:,12)+y(:,13)-(y(1,3)+y(1,4)+y(1,7)+y(1,8)+y(1,12)+y(1,13))))/N;

if abs(cum(1))<1e-9
    cum(1)=0;
else
    error('The cumulative number of newly infected at time 0 is not equal to 0');
end

%find indices of 3,6
ind3=find(t>3*30,1);
ind6=find(t>6*30,1);

%cumulative hospitalized following 3 and 6 months after the start of
%simulation
cum_hosp03=cum_hosp(ind3);
cum_hosp06=cum_hosp(ind6);

%outer loop: circles through r2: increase of the contact rate of the
%vaccinated
%figure counter
subfigc=1;
for r2=r2arr
    %the next loop: circles through mu1: sensitivity of the compliance
    %waning to vaccination coverage
    for mu1=mu1arr
        %allocate result arrays
        PeakRes=zeros(meshsize,meshsize);
        Cum3=zeros(meshsize,meshsize);
        Cum6=zeros(meshsize,meshsize);
        for i1=1:meshsize
            for i2=1:meshsize 
                upsilon0=Uarr(i1,i2);
                omega=Oarr(i1,i2);
                % in the next two lines the column index depends on the variant
                hosp_prob_vacc=hosp_prob_arr(2,2);
                hosp_prob_unvacc=hosp_prob_arr(1,2);

                pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
                [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
                % collect and process the output
                prev=y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16);
                cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;
                cum_hosp=100*(hosp_prob_unvacc*(y(:,3)+y(:,4)+y(:,7)+y(:,8)+y(:,12)+y(:,13)-(y(1,3)+y(1,4)+y(1,7)+y(1,8)+y(1,12)+y(1,13)))+hosp_prob_vacc*(y(:,16)-y(1,16)+y(:,17)-y(1,17)))/N;
                
                maxi1i2=max(prev);
                %find indices of 3,6,12,24 months
                ind3=find(t>3*30,1);
                ind6=find(t>6*30,1);
                %cumulative infected
                cum3=cum_hosp(ind3);
                cum6=cum_hosp(ind6);
                
                PeakRes(i1,i2)=maxi1i2;
                Cum3(i1,i2)=100*(cum3-cum_hosp03)/cum_hosp03;
                Cum6(i1,i2)=100*(cum6-cum_hosp06)/cum_hosp06;
            end
        end
        
        f1=figure(figc);
        subplot(4,2,2*subfigc-1);
        s=surf(Uarr,100*Oarr,Cum3);
        colormap jet
        
        f=subplot(4,2,2*subfigc-1);
        hp4 = get(subplot(4,2,2*subfigc-1),'Position')
        f.Position=[hp4(1) hp4(2) hp4(3),0.2];
        
        s.EdgeColor = 'none';
        hold on;
        view(2);
        [M,C]=contour(Uarr,100*Oarr,Cum3,[0,0]);
        C.LineColor='m';
        C.LineWidth=4;
        % colorbar will need to be updated, depending on the row of the
        % figure
        
        caxis([-40,65]);
        if subfigc==4
            set(gca,'xtick',xtpoints3)
            set(gca,'xticklabel',xtlabelsStr3)
            xlabel({'Vaccination coverage after three months (\%)'},'interpreter','latex');
        else 
            set(gca,'xtick',xtpoints3)
            set(gca,'xticklabel',[])
        end
        
         if r2==r2arr(1) & mu1==mu1arr(1)
             ylabel({'\textbf{Fast compliance decay,}';'\textbf{high contact rate}';'\textbf{of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         elseif r2==r2arr(1) & mu1==mu1arr(2)
             ylabel({'\textbf{Slow compliance decay,}';'\textbf{high contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         elseif r2==r2arr(2) & mu1==mu1arr(1)
             ylabel({'\textbf{Fast compliance decay,}';'\textbf{low contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         else
             ylabel({'\textbf{Slow compliance decay,}';'\textbf{low contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         end
         if subfigc==1
            title('Three months since vaccination start','interpreter','latex','FontSize',23);
         end
         set(gca,'FontSize',12);
       
        subplot(4,2,2*subfigc);
        s=surf(Uarr,100*Oarr,Cum6);
        colormap jet
        % colorbar will need to be updated, depending on the row of the
        caxis([-40,65]);
        s.EdgeColor = 'none';
        hold on;
        view(2);
        [M,C]=contour(Uarr,100*Oarr,Cum6,[0,0]);
        C.LineColor='m';
        C.LineWidth=4;
        
        f=subplot(4,2,2*subfigc);
        hp4 = get(subplot(4,2,2*subfigc),'Position')
        f.Position=[f.Position(1) f.Position(2) f.Position(3),0.2];
        
        %yline(60,'g','LineWidth',4);
        if subfigc==4
            set(gca,'xtick',xtpoints6)
            set(gca,'xticklabel',xtlabelsStr6)
            xlabel({'Vaccination coverage after six months (\%)'},'interpreter','latex');
        else 
            set(gca,'xtick',xtpoints6)
            set(gca,'xticklabel',[])
        end
        if subfigc==1
            title('Six months since vaccination start','interpreter','latex','FontSize',23);
        end
        set(gca,'FontSize',12);
        
        colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.04  hp4(2)*0.2+hp4(3)*0.1])
        
        subfigc=subfigc+1;
    end
end

figc=figc+1;
%%
% delta-like variant

R0delta=4.92;
epsilonDelta=R0delta*gamma/chat;
beta=epsilonDelta*c;

%generate no-vaccination statistics
upsilon0=0;

mu1=mu1arr(1);
omega=1;
r2=13.47/c;

pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
[t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
% collect outputs
cum_hosp=hosp_prob_arr(1,3)*(y(:,3)+y(:,4)+y(:,7)+y(:,8)+y(:,12)+y(:,13)-(y(1,3)+y(1,4)+y(1,7)+y(1,8)+y(1,12)+y(1,13)));

% maximum of prevalence
max0=max(prev);
%find indices of 3,6
ind3=find(t>3*30,1);
ind6=find(t>6*30,1);
%cumulative infected
cum_hosp03=cum_hosp(ind3);
cum_hosp06=cum_hosp(ind6);

XpopOut=N/popOut;

%outer loop: circles through r2: increase of the contact rate of the
%vaccinated
%figure counter
subfigc=1;
for r2=r2arr
    %the next loop: circles through mu1: sensitivity of the compliance
    %waning to vaccination coverage
    for mu1=mu1arr
        %allocate result arrays
        PeakRes=zeros(meshsize,meshsize);
        Cum3=zeros(meshsize,meshsize);
        Cum6=zeros(meshsize,meshsize);
        for i1=1:meshsize
            for i2=1:meshsize 
                upsilon0=Uarr(i1,i2);
                omega=Oarr(i1,i2);
                % in the next two lines the column index depends on the variant
                hosp_prob_vacc=hosp_prob_arr(2,3);
                hosp_prob_unvacc=hosp_prob_arr(1,3);

                pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega];
                [t,y]=ode45(@(t,y)COVIDVaccineRHS2(t,y,pars),[0,T], init,opts);
                % collect and process the output
                prev=y(:,2)+y(:,3)+y(:,6)+y(:,7)+y(:,11)+y(:,12)+y(:,15)+y(:,16);

                cum_hosp=hosp_prob_unvacc*(y(:,3)+y(:,4)+y(:,7)+y(:,8)+y(:,12)+y(:,13)-(y(1,3)+y(1,4)+y(1,7)+y(1,8)+y(1,12)+y(1,13)))+hosp_prob_vacc*(y(:,16)-y(1,16)+y(:,17)-y(1,17));
                cum=(y(:,2)+y(:,3)+y(:,4)+y(:,6)+y(:,7)+y(:,8)+y(:,11)+y(:,12)+y(:,13)+y(:,15)+y(:,16)+y(:,17))-infect0;
                maxi1i2=max(prev);
                %find indices of 3,6,12,24 months
                ind3=find(t>3*30,1);
                ind6=find(t>6*30,1);
                %cumulative infected
                cum3=cum_hosp(ind3);
                cum6=cum_hosp(ind6);
                
                PeakRes(i1,i2)=maxi1i2;
                Cum3(i1,i2)=100*(cum3-cum_hosp03)/cum_hosp03;
                Cum6(i1,i2)=100*(cum6-cum_hosp06)/cum_hosp06;
            end
        end
        
        f1=figure(figc);
        subplot(4,2,2*subfigc-1);
        s=surf(Uarr,100*Oarr,Cum3);
        colormap jet
        
        f=subplot(4,2,2*subfigc-1);
        hp4 = get(subplot(4,2,2*subfigc-1),'Position')
        f.Position=[hp4(1) hp4(2) hp4(3),0.2];
        
        s.EdgeColor = 'none';
        hold on;
        view(2);
        [M,C]=contour(Uarr,100*Oarr,Cum3,[0,0]);
        C.LineColor='m';
        C.LineWidth=4;
        % colorbar will need to be updated, depending on the row of the
        % figure
        
        caxis([-50,90]);
        if subfigc==4
            set(gca,'xtick',xtpoints3)
            set(gca,'xticklabel',xtlabelsStr3)
            xlabel({'Vaccination coverage after three months (\%)'},'interpreter','latex');
        else 
            set(gca,'xtick',xtpoints3)
            set(gca,'xticklabel',[])
        end
        
         if r2==r2arr(1) & mu1==mu1arr(1)
             ylabel({'\textbf{Fast compliance decay,}';'\textbf{high contact rate}';'\textbf{of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         elseif r2==r2arr(1) & mu1==mu1arr(2)
             ylabel({'\textbf{Slow compliance decay,}';'\textbf{high contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         elseif r2==r2arr(2) & mu1==mu1arr(1)
             ylabel({'\textbf{Fast compliance decay,}';'\textbf{low contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         else
             ylabel({'\textbf{Slow compliance decay,}';'\textbf{low contact rate of vaccinated}';' ';'Vaccine efficacy (\%)'},'interpreter','latex');
         end
         if subfigc==1
            title('Three months since vaccination start','interpreter','latex','FontSize',23);
         end
         set(gca,'FontSize',12);
       
        subplot(4,2,2*subfigc);
        s=surf(Uarr,100*Oarr,Cum6);
        colormap jet
        % colorbar will need to be updated, depending on the row of the
        caxis([-50,90]);
        s.EdgeColor = 'none';
        hold on;
        view(2);
        [M,C]=contour(Uarr,100*Oarr,Cum6,[0,0]);
        C.LineColor='m';
        C.LineWidth=4;
        
        f=subplot(4,2,2*subfigc);
        hp4 = get(subplot(4,2,2*subfigc),'Position')
        f.Position=[f.Position(1) f.Position(2) f.Position(3),0.2];
        
        %yline(60,'g','LineWidth',4);
        if subfigc==4
            set(gca,'xtick',xtpoints6)
            set(gca,'xticklabel',xtlabelsStr6)
            xlabel({'Vaccination coverage after six months (\%)'},'interpreter','latex');
        else 
            set(gca,'xtick',xtpoints6)
            set(gca,'xticklabel',[])
        end
        if subfigc==1
            title('Six months since vaccination start','interpreter','latex','FontSize',23);
        end
        set(gca,'FontSize',12);
        
        colorbar('Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.04  hp4(2)*0.2+hp4(3)*0.1])
        
        subfigc=subfigc+1;
    end
end
