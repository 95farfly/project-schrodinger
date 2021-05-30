function PSI = wav_func1(E,V,dx,N)
    H_bar = 1.055e-34;MassofElectron = 9.11e-31;Electronvolt = 1.602e-19;
    PSI = zeros(1,N); PSI(2) = 1;
    %================fdm shooting equation in a loop==========
    for n = 2:N-1   %finding every psi value for every x value
       kc = (2*MassofElectron*Electronvolt/H_bar^2)*(E-V(n))*dx^2;
       PSI(n+1) = (2 - kc)*PSI(n) - PSI(n-1); %this is the fdm shooting eqn
    end
    Area = trapz(PSI.*PSI); PSI = PSI /sqrt(Area); %calculating probability
end
