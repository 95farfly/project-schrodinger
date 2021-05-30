clc; 
close all;
clear ;
%=============================marlon keppitipola============================

%===========================Develop potential function ==================
tollarence = 1e-5;%since the end of the energy state only approache zero and is not exactly zero
V_upperboun = 10; V_lowerbound = 0 ;     
%===============================================defining potential well using arbitary well boundaries 
LengthofWell = 8e-10; Wellboundry = LengthofWell/2 ; N = 1000;
x = linspace(-LengthofWell,LengthofWell,N);  
dx = x(2)-x(1);%x and lenght of well in nanometers
V = zeros(1,length(x));
for i = 1:length(x)
    if -Wellboundry<=x(i) && x(i)<=Wellboundry
        V(i) = V_lowerbound;
    else
        V(i) = V_upperboun;
    end
end
figure(1);plot(x,V);
% =============================setting up Energy values for shooting method ==============================
EnergyState1 = 0;EnergyState2 = 10;divisions = 1/100000;  
%selecting minute energy values to test (step value)
E = EnergyState1:divisions:EnergyState2;
%===================== shooting method ==================================
for e = 1:length(E)
    PSI = wav_func1(E(e),V,dx,N);
    if(abs(PSI(end))<tollarence)%this is the shooting method to select the right energy value
        Ee(e) = E(e);
        figure(2)
        plot(x,PSI)
        grid on
        hold on
    end 
end
Eefirstvalue = unique(Ee)
for i = 2:length(Eefirstvalue)
    if Eefirstvalue(i)-Eefirstvalue(i-1)>0.5
        NewEvalue(i) = Eefirstvalue(i);
    end
end
unique(NewEvalue)
