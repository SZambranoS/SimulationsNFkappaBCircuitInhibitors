
clear

hold off;

close all;

%VectorParamModel0; 
VectorParamModelNS_derivingfroma;

save vectorparamusedthreefeedbacksandA20NS.mat vectorparam; 

v0=[0 Ntot 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

v=v0;


[tODE_0,DataODE_0] = ode23(@FlowDetermNFkappaB_targetGene_threefeedbacks_A20_NS,[0 36000],v0);

[nframes,m]=size(DataODE_0)

v0def=DataODE_0(nframes,:); 

S=1;

vectorparam(41)=S; 

save vectorparamusedthreefeedbacksandA20NS.mat vectorparam; 
 

[tODE_R,DataODE_R] = ode23(@FlowDetermNFkappaB_targetGene_threefeedbacks_A20_NS,[0 36000],v0def);


figure(3)
plot(tODE_R/3600,DataODE_R(:,2)/vectorparam(40),'r','linewidth',2);

%VectorParamModel0; 

VectorParamModelNS_derivingfroma;


vectorparam=GetPAramsBfromR_NS(vectorparam); 

save vectorparamusedthreefeedbacksandA20NS.mat vectorparam; 

v0=[0 vectorparam(40) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

[tODE_0,DataODE_0] = ode23(@FlowDetermNFkappaB_targetGene_threefeedbacks_A20_NS,[0 36000],v0);

[nframes,m]=size(DataODE_0)

v0def=DataODE_0(nframes,:); 

S=1;

vectorparam(41)=S; 




save vectorparamusedthreefeedbacksandA20NS.mat vectorparam; 
 
[tODE_B,DataODE_B] = ode23(@FlowDetermNFkappaB_targetGene_threefeedbacks_A20_NS,[0 36000],v0def);


figure(3)
hold on;
plot(tODE_B/3600,DataODE_B(:,2)/vectorparam(40),'b','linewidth',2);
xlabel('t (h)')
ylabel('Nuc:Tot NF-\kappaB')
set(gca,'fontsize',20);
axis([0 5 0 1])


figure(4)
plot(tODE_R/3600,DataODE_R(:,6),'r','linewidth',2);
hold on
plot(tODE_B/3600,DataODE_B(:,6),'b','linewidth',2);
xlabel('t (h)')
ylabel('I\kappaB\alpha_{RNA}')
set(gca,'fontsize',20);


figure(5)
plot(tODE_R/3600,DataODE_R(:,10),'r','linewidth',2);
hold on
plot(tODE_B/3600,DataODE_B(:,10),'b','linewidth',2);
xlabel('t (h)')
ylabel('I\kappaB\beta_{RNA}')
set(gca,'fontsize',20);


figure(6)
plot(tODE_R/3600,DataODE_R(:,14),'r','linewidth',2);
hold on
plot(tODE_B/3600,DataODE_B(:,14),'b','linewidth',2);
xlabel('t (h)')
ylabel('I\kappaB\epsilon_{RNA}')

set(gca,'fontsize',20);


