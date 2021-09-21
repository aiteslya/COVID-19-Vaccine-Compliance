function dydt=COVIDVaccineRHS2(t,y,pars)
% right hand-side of the system describing dynamics of COVID-19 in the
% population with compliance and vaccination process. Compliant individuals
% practice social distancing. Compliance raises as a function of incidence
% of infectious people and fades as a function of cumulative vaccinated

dydt=zeros(17,1);

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

%set parameters
%pars=[beta,r1,r2,delta,mu0,mu1,upsilon0,alpha,gamma,k1,k2,omega]
beta=pars(1);
r1=pars(2);
r2=pars(3);
delta=pars(4);
mu0=pars(5);
mu1=pars(6);
upsilon0=pars(7);
alpha=pars(8);%
gamma=pars(9);
k1=pars(10);
k2=pars(11);
omega=pars(12);

%define infection rates
lambdaInf=beta*(I+Ic*r1+(Iv+Infvacc)*r2)/(N+r1*Nc+r2*NV);
lambdaCInf=beta*(I*r1+Ic*r1^2+(Iv+Infvacc)*r1*r2)/(N+r1*Nc+r2*NV);
lambdaVInf=beta*(I*r2+Ic*r1*r2+(Iv+Infvacc)*(r2^2))/(N+r1*Nc+r2*NV);
%define compliance rise rate
lambdaC=delta*alpha*(E+Ec+Ev+Expvacc);
mu=mu0+mu1*Tv;
upsilon=upsilon0;
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