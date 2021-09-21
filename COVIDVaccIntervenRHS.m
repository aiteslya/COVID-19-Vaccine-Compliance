function dydt=COVIDVaccIntervenRHS(t,y,pars)
% this function defines the rhs of the ordinary differential equations
% system that models the dynamics of SARS-CoV-2 in a population where a
% vaccination rollout takes place and physical distancing measures are
% being imposed if the prevalence exceeds a threshold. As a result of the
% lockdown the contact rates of all individuals change

%parse parameters
%pars=[C,c,epsilon,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega,pc]
% C: pre-pandemic contact rates
% c: contact rate of non-compliant non-vaccinated individuals in the
% lockdown
% epsilon: probability of transmission upon a contact
%pc: critical level of prevalence above which a lockdown starts
C=pars(1);%strict lockdown contact rates
c=pars(2);
epsilon=pars(3);
r1=pars(4); %lockdown values
r2=pars(5); %lockdown values
delta=pars(6);
mu0=pars(7);
mu1=pars(8);
upsilon=pars(9);
alpha=pars(10);
gamma=pars(11);
k1=pars(12);
k2=pars(13);
omega=pars(14);
pc=pars(15);

%non-compliant
S=y(1);
E=y(2);
I=y(3);
R=y(4);
%compliant
Sc=y(5);
Ec=y(6);
Ic=y(7);
Rc=y(8);
%vaccinated
V=y(9);
Sv=y(10);
Ev=y(11);
Iv=y(12);
Rv=y(13);
Tv=y(14);

Expvacc=y(15);
Infvacc=y(16);
Recov=y(17);

%total population
N=sum(y(1:4));
Nc=sum(y(5:8));
NV=sum(y(9:13))+sum(y(15:17));

%calculate prevalence of infectious cases
prev=(I+Ic+Iv+Infvacc)/1.7e7;

c=(1-heaviside(prev-pc))*c+heaviside(prev-pc)*C;
r1=(1-heaviside(prev-pc))*r1+heaviside(prev-pc);
r2=(1-heaviside(prev-pc))*r2+heaviside(prev-pc);

%calculate beta
beta=c*epsilon;
%define infection rates
lambdaInf=beta*(I+Ic*r1+(Iv+Infvacc)*r2)/(N+r1*Nc+r2*NV);
lambdaCInf=beta*(I*r1+Ic*r1^2+(Iv+Infvacc)*r1*r2)/(N+r1*Nc+r2*NV);
lambdaVInf=beta*(I*r2+Ic*r1*r2+(Iv+Infvacc)*(r2^2))/(N+r1*Nc+r2*NV);
%define compliance rise rate
lambdaC=delta*alpha*(E+Ec+Ev+Expvacc);
mu=mu0+mu1*Tv;

%non-compliant
dydt(1,1)=-S*lambdaInf-S*lambdaC+mu*Sc-upsilon*S;
dydt(2,1)=S*lambdaInf-alpha*E-E*lambdaC+mu*Ec-upsilon*E;
dydt(3,1)=alpha*E-gamma*I-I*lambdaC+mu*Ic-k1*upsilon*I;
dydt(4,1)=gamma*I-R*lambdaC+mu*Rc-k2*upsilon*R;
%compliant
dydt(5,1)=-Sc*lambdaCInf+S*lambdaC-mu*Sc-upsilon*Sc;
dydt(6,1)=Sc*lambdaCInf-alpha*Ec+E*lambdaC-mu*Ec-upsilon*Ec;
dydt(7,1)=alpha*Ec-gamma*Ic+I*lambdaC-mu*Ic-k1*upsilon*Ic;
dydt(8,1)=gamma*Ic+R*lambdaC-mu*Rc-k2*upsilon*Rc;
%vaccinated
dydt(9,1)=upsilon*omega*(S+Sc);
dydt(10,1)=upsilon*(1-omega)*(S+Sc)-Sv*lambdaVInf;
dydt(11,1)=-alpha*Ev+upsilon*(E+Ec);
dydt(12,1)=alpha*Ev-gamma*Iv+k1*upsilon*(I+Ic);
dydt(13,1)=gamma*Iv+k2*upsilon*(R+Rc);
dydt(14,1)=upsilon*(S+E+k1*I+k2*R+Sc+Ec+k1*Ic+k2*Rc);

%infected while vaccinated
dydt(15,1)=Sv*lambdaVInf-alpha*Expvacc;
dydt(16,1)=alpha*Expvacc-gamma*Infvacc;
dydt(17,1)=gamma*Infvacc;
end
    