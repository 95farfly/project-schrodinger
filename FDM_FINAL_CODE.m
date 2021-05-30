clear all
close all
%====================================marlon===========================
%hello...simply uncomment each of the generated wells in the bottom or
%uncomment the entire electron atom section to generate the wells :)
%please note that this algorithm can only run one well at a time
% ------known constants------------------------------------------
Hbar = 1.055e-34;mass = 9.11e-31; electronVolt = 1.602e-19;
% ------varying constants-----------------------------
Lw = 8e-10;z_lower = -Lw;N = 1000;
z_upper = Lw;  z = linspace(z_lower,z_upper,N); d = z(2) - z(1); 
% ------potential energy matrix----------------------------
Wellheigth = 10; Vo = Wellheigth*electronVolt;                          %highest pottential of well
wellBoundary1 = -Lw/2; wellboundary2 = Lw/2; 
%=========================One-Electron Atom=============================
% Wellheigth = 0; epsi=8.854e-12; % Permittivity of free space and well height
% z_lower = 1e-10; % Minimum electron radius
% z = linspace(z_lower,z_upper,N);d = z(2) - z(1); 
% V = (-(electronVolt^2)/(4*pi*epsi))*(1./z); % Potential function for one electron
%=============================================================================================================
% -------generating different wells----------------------------------
%V((z<=wellBoundary1)) = Vo/(wellboundary2-Lw)^2*((z(z<=wellBoundary1)) + Lw).^2; V(z>= wellBoundary1 & z<=wellboundary2) = 0;V((z>=wellboundary2)) = Vo/(Lw-wellboundary2)^2*((z(z>=wellboundary2)) - Lw).^2; %Quadratic Well
%V = ones(1,N)*Vo; V(z>= wellBoundary1 & z<=wellboundary2) = 0; % Square Well
%V = ones(1,N)*Vo; V(z>= wellBoundary1 & z<=wellboundary2) = Vo/(wellboundary2^2)*(z((z>= wellBoundary1 & z<=wellboundary2)).^2); % Truncated Parabolic Well
%V((z<=wellBoundary1)) = Vo/(wellboundary2-Lw)^2*((z(z<=wellBoundary1)) - Lw).^2; V(z>= wellBoundary1 & z<=wellboundary2) = 0;V((z>=wellboundary2)) = Vo/(Lw-wellboundary2)^2*((z(z>=wellboundary2)) - Lw).^2; % unique Well
%--------------------------------------------------------------------------
Vn = eye(N,N);Vp = V'.*Vn;                            % ------kinetic energy matrix--------------------
Kinetic = eye(N,N);                                   %kinetic matrix is iniated using identity function
Kinetic = Kinetic*(-2);                                      % main diaganol matrix
for t = 1:N-1                                             % flanking the diagonals using a loop
   Kinetic(t,t+1) = 1;
   Kinetic(t+1,t) = 1;
end
                                                                  %  ------Hamiltonian----------------------------------------
Ham = (-(Hbar^2)/(2*mass*d^2))*Kinetic+Vp ;
                                                            % -----eigen values/states-------------------------
[phi,EnergyValue] = eig(Ham);                                            % using matlab function eig to obtain wavefunctions and the eigen energies
Prob = phi.^2;                                                % calculating the proability density 
E = diag(EnergyValue);                                                                         %designing column vector
Ec = E./(electronVolt);Ve = (1/electronVolt).*V;disp(E(1));disp(E(2));

% -------Plotting the Functions-------------------------------------------
%---------Plotting the allowable energy states-----------------------------
indices = find(Ec<Wellheigth); Ea = Ec(indices) 
figure
plot(z',Ve,'m','Linewidth',2);
hold on
for r = 1:length(Ea)
    yline(Ea(r),'k--','Linewidth',1.2);  
end
grid on
legend()

for int = 1:size(indices)

figure

yyaxis right
plot(z',Ve,'Linewidth',2);
title('wavefunction')
xlabel('well boundary') 

ylabel('potential function')
hold on
yyaxis left 
plot(z',[phi(:,int)],'g','Linewidth',1);

grid on
 
figure

 
yyaxis right 
plot(z',Ve,'Linewidth',2) ;
title('probability density')
xlabel('well boundary') 
ylabel('potential function')

hold on
yyaxis left 
plot(z',[Prob(:,int)],'g','Linewidth',1);

grid on

end
