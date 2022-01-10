%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Project           : Optimisation of a BLDC motor stator
% 
%  Program name      : stator.m
% 
%  Author            : Colin Laganier
% 
%  Date created      : 12/12/2020
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
close all 

% Variables
P = 4; %Num Rotor Poles
Z = 5; %Num Stator Teeth
Na = 60; %Num coil windings
Wm = 0.010; %Height of stator
r = 0.036; %Radius of stator

% Constants
I = 15; %Armature Current
A = 2; %Num of parallel paths - 2 because wave winding
Br = 1.2; %Residual Magnetism of NdFeB N35
Dm = 0.0025; %Magnet thickness
z = 0.002; %Rotor-stator gap
Pi = pi; %Value of Pi
D = 500*r; %Diameter of stator
rw = 0.00051; %Thickness (radius) of copper coil (18AWG) - m
Dw = 8960; %Density of copper windings - kg/m3
Ds = 7650; %Density of laminated steel
rb = 0.0125; %Stator bore hole radius
alpha = 0.5; %
beta = 0.8;
Wt = 0.003; %Stator tooth end thickness

% Calculate dimensions of magnets
if P < 4
    Lm = (2*Pi*(r))/4;
else
    Lm = ((2*Pi*(r))/P)*0.8;
end
% Permanent magnet flux density
Bm = (Br/Pi)*(atan((Lm*Wm)/(2*z*sqrt(4*z.^2+Lm.^2+Wm.^2)))-atan((Lm*Wm)/(2*(Dm+z)*sqrt(4*(Dm+z).^2+Lm.^2+Wm.^2))));
% Flux per pole
Fpp = (2*Bm*D*Wm)/Z; 
% Magnetic Torque
T = (P*Z*Na*Fpp*I)/(2*Pi*A);
% Stator mass
Ms = Ds*Wm*(((Pi*r.^2)/2)-Pi*rb.^2+Z*((((alpha*2*Pi*r)/Z)*(0.25*r-Wt))+((Wt*beta*2*Pi*r)/2)));
% Coil mass
Mc = Z*Na*Pi*Dw*(2*(0.001+((alpha*2*Pi*r)/Z))+2*(Wm+0.001))*rw.^2;

% Objective function
Fobj = T/(Ms+Mc);

disp(['Torque: ' num2str(T) 'Nm'])
disp(['Total Mass: ' num2str(Mc+Ms) 'kg'])
disp(['Torque to Mass ratio: ' num2str(Fobj) 'Nm/kg'])
