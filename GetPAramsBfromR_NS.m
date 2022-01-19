function vectorparam = GetPAramsFfromK_NS(vectorparam0)

vectorparam=vectorparam0; 

%Now we apply the transformations we have seen


k1Na= vectorparam(10);
k0Ia= vectorparam(11);

NtotK=vectorparam(40); 

NtotF=1.4*NtotK;

factorF=(k1Na*NtotF)/(k1Na*NtotF+k0Ia);
factorK=(k1Na*NtotK)/(k1Na*NtotK+k0Ia);




vectorparam(7)=vectorparam(7)*0.809*(factorK/factorF); %Nfkbia smaller

vectorparam(20)=vectorparam(20)*(0.0248/0.018)*(factorK/factorF); %Nfkbib bigger

vectorparam(33)=vectorparam(33)*(0.089/0.116)*(factorK/factorF);%Nfkbie smaller

vectorparam(40)=NtotF; %Ntot bigger

vectorparam(42)=1.035*vectorparam(42); %KK bigger
