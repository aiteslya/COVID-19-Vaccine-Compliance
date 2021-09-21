% this script loads the cleaned up vaccination uptake data for the Belarus, Denmark, Israel and the Netherlands
% from http://ourworldindata.org
% both fully vaccinated individuals and individuals vaccinated with at
% least one dose
% the script approximates the vaccination rates 

%prepare the space
clc;
clear variables;
close all;
format long;

% read tables

%read the CSV files and load them into array
TblBEL = readtable('DataVaccine/Belarus.csv');
TblDEN = readtable('DataVaccine/Denmark.csv');
TblISR = readtable('DataVaccine/Israel.csv');
TblNED = readtable('DataVaccine/Netherlands.csv');

% convert dates to serial numbers and set the first date in the date set as
% zero
dates=TblBEL.date;
DATEBEL=datenum(dates)-datenum(dates(1))+1;

dates=TblDEN.date;
DATEDEN=datenum(dates)-datenum(dates(1))+1;

dates=TblISR.date;
DATEISR=datenum(dates)-datenum(dates(1))+1;

dates=TblNED.date;
DATENED=datenum(dates)-datenum(dates(1))+1;

% we are interested in the date, people vaccinated and
% convert dates into number of days since the first vaccination
% people_fully_vaccinated columns
TBLBEL=[DATEBEL';(TblBEL.people_vaccinated)';(TblBEL.people_fully_vaccinated)';]';
TBLBEL(:,1)=TBLBEL(:,1)-1;
TBLDEN=[DATEDEN';(TblDEN.people_vaccinated)';(TblDEN.people_fully_vaccinated)';]';
TBLDEN=[0 0 0 ;TBLDEN];
TBLISR=[DATEISR';(TblISR.people_vaccinated)';(TblISR.people_fully_vaccinated)';]';
TBLISR=[0 0 0 ;TBLISR];
TBLNED=[DATENED';(TblNED.people_vaccinated)';(TblNED.people_fully_vaccinated)';]';
TBLNED=[0 0 0 ;TBLNED];
% replace NANs with zeros in Denmark and Israel Data

TBLISR(isnan(TBLISR(:,3)),3)=0;
TBLDEN(isnan(TBLDEN(:,3)),3)=0;

%for each country plot total vaccinated and total fully vaccinated

figure(1);
subplot(2,4,1);
h(1)=plot(TBLBEL(:,1),TBLBEL(:,2),'-ro','MarkerSize',4);hold on;
h(2)=plot(TBLBEL(:,1),TBLBEL(:,3),'-r^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('Individuals','interpreter','latex');
title('Belarus');
set(gca,'FontSize',25);

subplot(2,4,2);
h(1)=plot(TBLDEN(:,1),TBLDEN(:,2),'-go','MarkerSize',4);hold on;
h(2)=plot(TBLDEN(:,1),TBLDEN(:,3),'-g^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('Individuals','interpreter','latex');
title('Denmark');
set(gca,'FontSize',25);

subplot(2,4,3);
h(1)=plot(TBLISR(:,1),TBLISR(:,2),'-bo','MarkerSize',4);hold on;
h(2)=plot(TBLISR(:,1),TBLISR(:,3),'-b^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('Individuals','interpreter','latex');
title('Israel');
set(gca,'FontSize',25);

subplot(2,4,4);
h(1)=plot(TBLNED(:,1),TBLNED(:,2),'-mo','MarkerSize',4);hold on;
h(2)=plot(TBLNED(:,1),TBLNED(:,3),'-m^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('Individuals','interpreter','latex');
title('Netherlands');
set(gca,'FontSize',25);

subplot(2,4,5);
h(1)=semilogy(TBLBEL(:,1),TBLBEL(:,2),'-ro','MarkerSize',4);hold on;
h(2)=semilogy(TBLBEL(:,1),TBLBEL(:,3),'-r^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('$$\log$$ Individuals','interpreter','latex');
legend(h,'Vaccinated with at least one dose','Fully vaccinated','interpreter','latex','location','southoutside');
title('Belarus');
set(gca,'FontSize',25);

subplot(2,4,6);
h(1)=semilogy(TBLDEN(:,1),TBLDEN(:,2),'-go','MarkerSize',4);hold on;
h(2)=semilogy(TBLDEN(:,1),TBLDEN(:,3),'-g^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('$$\log$$ Individuals','interpreter','latex');
legend(h,'Vaccinated with at least one dose','Fully vaccinated','interpreter','latex','location','southoutside');
title('Denmark');
set(gca,'FontSize',25);

subplot(2,4,7);
h(1)=semilogy(TBLISR(:,1),TBLISR(:,2),'-bo','MarkerSize',4);hold on;
h(2)=semilogy(TBLISR(:,1),TBLISR(:,3),'-b^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('$$\log$$ Individuals','interpreter','latex');
legend(h,'Vaccinated with at least one dose','Fully vaccinated','interpreter','latex','location','southoutside');
title('Israel');
set(gca,'FontSize',25);

subplot(2,4,8);
h(1)=semilogy(TBLNED(:,1),TBLNED(:,2),'-mo','MarkerSize',4);hold on;
h(2)=semilogy(TBLNED(:,1),TBLNED(:,3),'-m^','MarkerSize',4);
xlabel('Time (days)','interpreter','latex');
ylabel('$$\log$$ Individuals','interpreter','latex');
legend(h,'Vaccinated with at least one dose','Fully vaccinated','interpreter','latex','location','southoutside');
title('Netherlands');
set(gca,'FontSize',25);

% estimate vaccination rates
%define the population size
%Nsize=[Belarusia,Denmark, Israel, Netherlands];
Nsize=[9449323,5792202,9053000,17134872];
% transmorm the vaccinated data to density
%1-V(t)/N=exp(-upsilon*t)
TBLBEL(:,4)=1-TBLBEL(:,2)/Nsize(1);
TBLDEN(:,4)=1-TBLDEN(:,2)/Nsize(2);
TBLISR(:,4)=1-TBLISR(:,2)/Nsize(3);
TBLNED(:,4)=1-TBLNED(:,2)/Nsize(4);
%Belarus
% recalculate plot
f = fit(TBLBEL(:,1),TBLBEL(:,4),'exp1','StartPoint',[1 0.0007178])
mBel=0.0007178;
%mBel=(TBLBEL(ind2,2)-TBLBEL(ind1,2))/(TBLBEL(ind2,1)-TBLBEL(ind1,1));
figure(1);subplot(2,4,1);
plot(TBLBEL(:,1),Nsize(1)*(1-exp(-mBel*TBLBEL(:,1))),'r','LineWidth',2);

%Denmark
f = fit(TBLDEN(:,1),TBLDEN(:,4),'exp1','StartPoint',[1 0.00392])
% recalculate plot
mDen=0.00392;
%mBel=(TBLBEL(ind2,2)-TBLBEL(ind1,2))/(TBLBEL(ind2,1)-TBLBEL(ind1,1));
figure(1);subplot(2,4,2);
plot(TBLDEN(:,1),Nsize(2)*(1-exp(-mDen*TBLDEN(:,1))),'g','LineWidth',4);

%Israel
f = fit(TBLISR(:,1),TBLISR(:,4),'exp1','StartPoint',[1 0.005293])
mIsr=0.005293;
figure(1);subplot(2,4,3);
plot(TBLISR(:,1),Nsize(3)*(1-exp(-mIsr*TBLISR(:,1))),'b','LineWidth',4);

%Netherlands
ind1=1;
f = fit(TBLNED(:,1),TBLNED(:,4),'exp1','StartPoint',[1 0.004648])
% recalculate plot
mNed=0.004648;
%mBel=(TBLBEL(ind2,2)-TBLBEL(ind1,2))/(TBLBEL(ind2,1)-TBLBEL(ind1,1));
figure(1);subplot(2,4,4);
plot(TBLNED(:,1),Nsize(4)*(1-exp(-mNed*TBLNED(:,1))),'m','LineWidth',4);


marr=[mBel,mDen,mIsr,mNed];
[m,ind]=max(marr);
countries={'Belarus','Denmark','Israel','The Netherlands','Supremum'};
disp(['Maximal vaccination rate was ',num2str(m,2),' individuals per day and was detected in ',countries{ind},'.']);

marr=[marr,1e-2];
%output vaccination curves corresponding to each rate
num=100;
t=linspace(0,200,num);
figure(2);
col='rgbmc';
N=1.7e7;
for counter=1:numel(marr)
    m=marr(counter);
    plot(t,100*(1-exp(-m*t)),'color',col(counter),'LineWidth',4);hold on;
end    
xlabel('Time (days)','interpreter','latex');
ylabel('Vaccinated individuals','interpreter','latex');
legend('Belarus','Denmark','Israel','Netherlands','Supremum');
set(gca,'FontSize',25);
