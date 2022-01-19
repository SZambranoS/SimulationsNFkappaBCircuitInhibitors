

kaa=1.900000e-06; 
kda=8.400000e-04; 
dIa=6.700000e-05; 
dCa =1.340000e-05; 
dIaS= 5.000000e-09; 
dCaS =10*2.500000e-08; 
kRa =2.000000e-01;
dRa =7.500000e-04;
kIa =(2.500000e-01);
%k1Na =6.900000e-08; 
%k0Ia =4e4*1.400000e-08; %CORRECTED FOR NO STRIPPING NS

k1Na =6.900000e-08; 
k0Ia =3e4*1.400000e-08; %CORRECTED FOR NO STRIPPING NS

k0a =0.000000e+00;
k1a =0.000000e+00; %+13


paramsIa=[kaa, kda, dIa, dCa, dIaS, dCaS, kRa, dRa, kIa, k1Na, k0Ia, k0a, k1a];

paramsIb=paramsIa;
paramsIb(7)=0.018*paramsIa(7); 

paramsIe=paramsIa; 
paramsIe(7)=paramsIa(7)*0.116;


Ntot=3*10^4;

S=0;%+2
KK=(10^5)*2.1e-4; 
dK=2.1e-4; 
A0=3e4; 
n=1; 

paramsS=[S,KK,dK,A0,n];

dA1=paramsIa(3);
kA1=paramsIa(9); 
k1A1=paramsIa(10);
k0A1=paramsIa(11);
kRA1=paramsIa(7);
dRA1=paramsIa(8); 


paramsA=[dA1,kA1,k1A1,k0A1,kRA1,dRA1]; 

vectorparam0=[paramsIa,paramsIb,paramsIe, Ntot,paramsS,paramsA];

kNt=1*k1Na;

kIt=1*k0Ia;

kRt=kRa;

dRt=1*dRa; 

vectorgene=[kNt,kIt,kRt,dRt]; 

vectorparam=[vectorparam0,vectorgene]; 
