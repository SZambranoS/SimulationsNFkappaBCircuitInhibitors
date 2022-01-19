function y= FlowDetermNFkappaB_targetGene_threefeedbacks_A20_NS(t,x)



%TO FACILITATE THE FEEDBACK TUNING WE TAKE OUT THE STRIPPING MARCH 2021

load  vectorparamusedthreefeedbacksandA20NS.mat


%Parameters for the Ikba feedback loop
kaa=vectorparam(1);
kda=vectorparam(2);
dIa=vectorparam(3);
dCa= vectorparam(4);
dIaS= vectorparam(5);
dCSa= vectorparam(6);
kRa=vectorparam(7);
dRa= vectorparam(8);
kIa= vectorparam(9);
k1Na= vectorparam(10);
k0Ia= vectorparam(11);
k0a= vectorparam(12);
k1a= vectorparam(13);


%Parameters for the Ikbb feedback loop

kab=vectorparam(14);
kdb=vectorparam(15);
dIb=vectorparam(16);
dCb= vectorparam(17);
dIbS= vectorparam(18);
dCSb= vectorparam(19);
kRb=vectorparam(20);
dRb= vectorparam(21);
kIb= vectorparam(22);
k1Nb= vectorparam(23);
k0Ib= vectorparam(24);
k0b= vectorparam(25);
k1b= vectorparam(26);


%Parameters for the Ikbe feedback loop

kae=vectorparam(27);
kde=vectorparam(28);
dIe=vectorparam(29);
dCe= vectorparam(30);
dIeS= vectorparam(31);
dCSe= vectorparam(32);
kRe=vectorparam(33);
dRe= vectorparam(34);
kIe= vectorparam(35);
k1Ne= vectorparam(36);
k0Ie= vectorparam(37);
k0e= vectorparam(38);
k1e= vectorparam(39);


Ntot=vectorparam(40);
S=vectorparam(41);

KK=vectorparam(42);
dK=vectorparam(43);
A0=vectorparam(44);
n=vectorparam(45);
dA1=vectorparam(46); 
kA1=vectorparam(47); 
k1A1=vectorparam(48);
k0A1=vectorparam(49);
kRA1=vectorparam(50);
dRA1=vectorparam(51);


kNt=vectorparam(52);
kIat=vectorparam(53);
kRt=vectorparam(54);
dRt=vectorparam(55);

K=x(1); %Kinase

N=x(2); %Free NF-kB

Ia=x(3); %Ikba

Ca=x(4); %NFkB Ikba complex

Ga=x(5); %Gene Ikba

Ra=x(6); %RNA IkBa

Ib= x(7); %Ikbb 

Cb=x(8); %NFkB Ikbb complex

Gb=x(9); %Gene Ikbb

Rb=x(10); %RNA IkBb

Ie= x(11); %Ikbe 

Ce=x(12); %NFkB Ikbe complex

Ge=x(13); %Gene Ikbe

Re=x(14); %RNA IkBe

A=x(15); %A20

GA=x(16); %Gene A20

RA=x(17); %RNA A20

Gt=x(18); %Target gene

Rt=x(19); %Another target gene



fK=(KK*S/(1+((A/A0)^n)))-dK*K; 

fN=-kaa*N*Ia+kda*Ca+dCa*Ca+dCSa*K*Ca -kab*N*Ib+kdb*Cb+dCb*Cb+dCSb*K*Cb -kae*N*Ie+kde*Ce+dCe*Ce+dCSe*K*Ce;

fIa=-kaa*N*Ia+kda*Ca-dIa*Ia-dIaS*K*Ia+kIa*Ra;

fCa=kaa*N*Ia - kda*(Ca) - dCa*Ca - dCSa*K*Ca;

fGa=k1Na*N*(2-Ga)-k0Ia*Ga; %MODEL WITH NO STRIPPING

fRa=kRa*Ga-dRa*Ra; 

fIb=-kab*N*Ib+kdb*Cb-dIb*Ib-dIbS*K*Ib+kIb*Rb;

fCb=kab*N*Ib - kdb*(Cb) - dCb*Cb - dCSb*K*Cb;

fGb=k1Nb*N*(2-Gb)-k0Ib*Gb; %MODEL WITH NO STRIPPING 

fRb=kRb*Gb-dRb*Rb; 

fIe=-kae*N*Ie+kde*Ce-dIe*Ie-dIeS*K*Ie+kIe*Re;

fCe=kae*N*Ie - kde*(Ce) - dCe*Ce - dCSe*K*Ce;

fGe=k1Ne*N*(2-Ge)-k0Ie*Ge; %MODEL WITH NO STRIPPING

fRe=kRe*Ge-dRe*Re; 

fA=-dA1*A+kA1*RA;

fGA=k1A1*N*(2-GA)-k0A1*GA; 

fRA=kRA1*GA-dRA1*RA; 

fGt=kNt*N*(2-Gt)-kIat*Gt; %MODEL WITH NO STRIPPING

fRt=kRt*Gt-dRt*Rt; 


y=[fK,fN,fIa,fCa,fGa,fRa,fIb,fCb,fGb,fRb,fIe,fCe,fGe,fRe,fA,fGA,fRA,fGt,fRt]';
