%% Stator Pattern Search

clc
close all
clear

x0 = [8,6,10,0.01,0.025];     % initial values
lb=[2,2,2,0.01,0.015];       % lower bounds
ub=[18,24,60,0.06,0.05];          % upper bounds
A =[]; b=[]; Aeq=[]; beq=[];

nlcon = @(x)nonlcon(x,false);
nlconInit = @(x)nonlcon(x,true);
Tobj = @(x)objective1(x);
Mobj = @(x)objective2(x,true);
 
% options = optimoptions('fmincon','Algorithm', 'interior-point');
% options = optimoptions('fmincon','Algorithm', 'active-set');
options = optimoptions('fmincon','Algorithm', 'sqp');

xSet=zeros(2,5);
xSet(1,:) = fmincon(Tobj,x0,A,b,Aeq,beq,lb,ub,nlconInit,options);
xSet(2,:) = fmincon(Mobj,x0,A,b,Aeq,beq,lb,ub,nlconInit,options);

opts_ps = optimoptions('paretosearch','PlotFcn','psplotparetof','ParetoSetSize',160,'InitialPoints',xSet);
[x_psx0,fval_psx0] = paretosearch(@objval,5,A,b,Aeq,beq,lb,ub,nlcon,opts_ps);

%% Pareto Front Processing 
%normalizing the values and identifying the minimum ratio and outputting
%the optimal parameters

fnorm = mapstd(fval_psx0');  

multifun = fnorm(1,:)/fnorm(2,:); 

minim = find(multifun==min(multifun)); 

X = x_psx0(minim,:);
F = fval_psx0(minim,:);

disp(table(X(1),X(2),X(3),X(4),X(5),'VariableNames',{'P', 'Z', 'Na', 'Wm', 'r'}))
disp(['Final Objective Torque: ' num2str(-F(1)) 'Nm'])
disp(['Final Objective Mass: ' num2str(F(2)) 'kg'])
disp(['Final Objective Ratio: ' num2str(-F(1)/F(2)) 'Nm/kg'])

%% Objective function for Pattern Search

function Fun = objval(x)
    Fun(1) = objective1(x);
    Fun(2) = objective2(x,false);
end

%% Objective Function for Fmincon

function T=objective1(x)
    D = 500*x(5); %Diameter of stator
    I = 15; %Armature Current
    A = 2; %Num of parallel paths - 2 because wave winding
    Br = 1.2; %Residual Magnetism of NdFeB N35
    Dm = 0.0025; %Magnet thickness
    z = 0.002; %Rotor-stator gap
    Pi = pi; %Value of Pi
   
    if x(1) < 4
        Lm = (2*Pi*(x(5)))/4;
    else
        Lm = ((2*Pi*(x(5)))/x(1))*0.8;
    end
    % Permanent magnet flux density
    Bm = (Br/Pi)*(atan((Lm*x(4))/(2*z*sqrt(4*z.^2+Lm.^2+x(4).^2)))-atan((Lm*x(4))/(2*(Dm+z)*sqrt(4*(Dm+z).^2+Lm.^2+x(4).^2))));
    % Flux per pole
    Fpp = (2*Bm*D*x(4))/x(2);
    % Magnetic Torque
    T = -(x(1)*x(3)*Fpp*I)/(2*Pi*A);
end

function M=objective2(x,init)
    Pi = pi; %Value of Pi
    rw = 0.00051; %Thickness (radius) of copper coil (18AWG) - m
    Dw = 8960; %Density of copper windings - kg/m3
    Ds = 7650; %Density of laminated steel
    rb = 0.0125; %Stator bore hole radius
    alpha = 0.4; %
    beta = 0.8;
    Wt = 0.003; %Stator tooth end thickness

    % Stator mass
    Ms = Ds*x(4)*(((Pi*x(5).^2)/2)-Pi*rb.^2+x(2)*((((alpha*2*Pi*x(5))/x(2))*(0.25*x(5)-Wt))+((Wt*beta*2*Pi*x(5))/2)));
    % Coil mass
    Mc = x(2)*x(3)*Pi*Dw*(2*(0.001+((alpha*2*Pi*x(5))/x(2)))+2*(x(4)+0.001))*rw.^2;

    if init
        M = -(Ms+Mc);
    else 
        M = Ms+Mc;
    end
end

%% Non Linear Constraints

function [c, ceq] = nonlcon(x,initial) 

    ceq = [];
    
    rx = 1000*(0.25*x(5)-0.003);
    Q = fix(x(3)/floor(rx));
    R = rem(x(3),floor(rx));
    if R > 0
        Q = Q+1;
    end
    g1 = Q-(((0.4*pi*x(5))/x(2))*1000);
    g2 = x(1)-(x(2)-1);
    if initial == true
        c = [g1,g2];
    else 
        T = objective1(x);
        g3 = 2+T;
        M = objective2(x,true);
        g4 = -M-1.5;
        c = [g1,g2,g3,g4];
    end

end