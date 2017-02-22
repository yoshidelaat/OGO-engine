%% OGO 7
clear all;
close all;
clc;

%%some units
mm = 10^-3;   bar = 10^5;

%% Engine parameters set in script EngineDimensions.m
EngineDimensions

%% starting conditions
p0 = 1.01*bar;                                                             % Atmospheric pressure [Pa]
T0 = 298.15;                                                               % Temperature at start [K]
pdifference = 0.0*bar;                                                     % Estimated pressure difference between the piston and the atmosphere
T_wall =T0;                                                                % Cylinder wall temperature [K];


%% Load Fuel
Fuel

%% Number of extra cycles
k=1;        

%% Set theta
dtheta = 0.001 * pi;
theta = 0 : dtheta : 4*pi + 4*pi*k;                                               %theta=0 on intake open, theta=4pi after one cycle

%% Kinematics
s    = (l+S/2) - (S/2*cos(theta)+sqrt(l^2-(S/2)^2*(sin(theta)).^2));       % distance piston-cylinderhead
V    = Vc + s * pi * (B/2)^2;                                              % total volume in cylinder                         
A    = 2*pi*B/2*s + 2*(pi*(B/2)^2);                                        % effective inner area cylinder           
RPM  = 3000;                                                               % rotational speed [rounds per minute]

%% Wiebe Function 
a = 5;                                                                      % constant (Handbook)
n = 3;                                                                      % constant (Handbook)
theta_s = 6.152286;                                                         % start combustion (=ignition angle)
theta_d = 15/720*4*pi;                                                      % duration combustion (estimation)
x_b = 1 - exp(-a * ((theta - theta_s)/theta_d) .^ n);                       % mass-fraction burned fuel



%% Timing                                                                   %the timing measured in degrees put into radians                  
IntakeStart = 0;                                                    %4*k*pi was used to run multiple cycles, where k is the number of extra cycles
IntakeValveOpen = 15/720*4*pi;
IntakeValveClose = 225/720*4*pi; 
IgnitionStart = 345/720*4*pi;
CombustionStop = IgnitionStart;
ExhaustValveOpen = 480/720*4*pi;
ExhaustValveClose = 705/720*4*pi;
ExhaustStop = 4*pi;

%% Starting conditions
p(1)= p0;
T(1)= T0;
pm(1)=p0;                                                                   %motored pressure (for Woschni model)
Tm(1)=T0;                                                                   %motored temperature (for Woschni model)
RefState= T(1)/p(1)*V(1);                                                   %reference state (for Woschni model)

m(1)=p(1) * V(1) / (Runiv/Mmix_before) * T(1);                                         
dQloss(1)=0;


%% Loop veriabele
theta_loop = 0;
for i= 2 : length(theta)   % i goes from 1 to the length of theta


%% Intake with intake valve closed (closed system, expansion)
if theta_loop <= IntakeValveOpen     
    dV = V(i) - V(i-1);    
    C1=2.28;                                                                
    C2=0;                                                                   %Woschni constants
    dQloss(i)= 0;%Woschni(p(i-1),pm(i-1),T(i-1),V(i-1),A(i-1),C2,RPM,T_wall,RefState); %Woschni model for heat loss to wall (for comments, see 'Woschni.m')
       
    dT=(-dQloss(i-1)-p(i-1)*dV)/cv_1/m(i-1);                                %1st law of closed system
    T(i)=T(i-1)+dT;
    m(i) = m(i-1);                                                          %mass conservation
    p(i)=m(i)*Runiv/Mmix_before*T(i)/V(i);                                  %ideal gas law
    
    %motored cycle values:
    pm(i)=p(i);                                                             
    Tm(i)=T(i);
end

%% Intake with intake valve open   (open system, no heat loss to the wall assumed) 
if  theta_loop <= IntakeValveClose && theta_loop > IntakeValveOpen
    T(i) = T(i-1);                                                          %constant temperature
    p(i) = p(i-1);                                                          %constant pressure
    m(i) = p(i) * V(i) / ((Runiv/Mmix_before) * T(i));                      %ideal gas law
    m_fuel = (p0 * V(length(find(theta < IntakeValveClose))) / ((Runiv/Mmix_before) * T0))  * (1/(AF_stoic + 1));
                                                                            %total mass of the fuel in the mixture    
    %motored cycle values:
    pm(i)=p(i);                                                             
    Tm(i)=T(i);             
end

    
%% Expansion and compression with closed valves (closed system)
if theta_loop > IntakeValveClose && theta_loop <= IgnitionStart
    dV = V(i) - V(i-1);
    
    C1=2.28;                                                                
    C2=0;                                                                   %Woschni constants
    dQloss(i)= 0; %Woschni(p(i-1),pm(i-1),T(i-1),V(i-1),A(i-1),C1,C2,RPM,T_wall,RefState); %Woschni model for heat loss to wall (for comments, see 'Woschni.m')
    
    dT=(-dQloss(i)-p(i-1)*dV)/cv_1/m(i-1);
    T(i)=T(i-1)+dT;
    m(i) = m(i-1);    
    p(i)=m(i)*Runiv/Mmix_before*T(i)/V(i);
    m_fuel      = m(i) * (1/(AF_stoic+1));                                       % Mass fuel
   
    %motored cycle values:
    pm(i)=p(i);                                                             
    Tm(i)=T(i);
end
      

    
%% Ignition and combustion (closed system)
if theta_loop > IgnitionStart && theta_loop <= CombustionStop
 
    for j = 1:NSp
        cv(j) = CvNasa(T(length(find(theta < IgnitionStart))),Sp(iSp(j)));   %loading cv and cp from nasa                                      % Determine cv for T0 for each element with Nasa
        cp(j) = CpNasa(T(length(find(theta < IgnitionStart))),Sp(iSp(j)));   % Determine cv for T0 for each element with Nasa
    end
    cv_2    = cv*Yi_before';                                                 % cv mix
    cp_2    = cp*Yi_before';                                                 % cp mix
    gamma_2 = cp_2/cv_2;                                                     % calculating gamma 
    
    dV = V(i) - V(i-1);                               

    C1=2.28;                                                                
    C2=3.28e-3;                                                              %Woschni constants
    dQloss(i)= Woschni(p(i-1),pm(i-1),T(i-1),V(i-1),A(i-1),C1,C2,RPM,T_wall,RefState); %Woschni model for heat loss to wall (for comments, see 'Woschni.m')
    
    dm(i) = m_fuel * (x_b(i) - x_b(i-1));                                   % burned mass fraction with Wiebefunction x_b
    dQcomp(i) = Q_LHV_gasoline * dm(i);                                     % released head from burned mass
  
    dT=(-dQloss(i)+dQcomp(i)-p(i-1)*dV)/cv_2/m(i-1);                       
    T(i)=T(i-1)+dT;                                                        
    m(i) = m(i-1);                                                          
    p(i)=m(i)*Runiv/Mmix_before*T(i)/V(i);
    
    %motored cycle values (dQcomp=0)   
    dTm=(-dQloss(i)-pm(i-1)*dV)/cv_1/m(i-1);
    Tm(i)=Tm(i-1)+dTm;
    pm(i)= m(i)*Runiv/Mmix_before*Tm(i)/V(i);
end


%% expansion with valves closed (closed system)
if theta_loop > CombustionStop && theta_loop <= ExhaustValveOpen
    
    dV = V(i) - V(i-1);
    
    C1=2.28;                                                                
    C2=3.28e-3;                                                              %Woschni constants   
    dQloss(i)= 0;%Woschni(p(i-1),pm(i-1),T(i-1),V(i-1),A(i-1),C1,C2,RPM,T_wall,RefState); %Woschni model for heat loss to wall (for comments, see 'Woschni.m')
       
    dT=(-dQloss(i)-p(i-1)*dV)/cv_1/m(i-1);                                      
    T(i)=T(i-1)+dT;
    m(i) = m(i-1);    
    p(i)=m(i)*Runiv/Mmix_before*T(i)/V(i); 
    
    %motored cycle values:
    pm(i)=p(i);                                                             
    Tm(i)=T(i);
end

%% expansion and compression with exhaust valve open (open system, no heat loss to the wall assumed)
if theta_loop > ExhaustValveOpen && theta_loop <= ExhaustValveClose
       T(i) = T(i-1);
       p(i) = p(i-1);
       m(i) = p(i) * V(i) / ((Runiv/Mmix_before) * T(i));
    
       %motored cycle values:
       pm(i)=p(i);                                                             
       Tm(i)=T(i);
end

%% compression with exhaust valve closed (closed system)
if theta_loop > ExhaustValveClose %&& i <= length(find(theta < ExhaustStop))
    dV = V(i) - V(i-1);                                         
    
    C1=2.28;                                                                
    C2=0;                                                                   %Woschni constants
    dQloss(i)= 0;%Woschni(p(i-1),pm(i-1),T(i-1),V(i-1),A(i-1),C1,C2,RPM,T_wall,RefState); %Woschni model for heat loss to wall (for comments, see 'Woschni.m')
       
    dT=(-dQloss(i-1)-p(i-1)*dV)/cv_1/m(i-1);                                %1st law of closed system
    T(i)=T(i-1)+dT;
    m(i) = m(i-1);                                                          %mass conservation
    p(i)=m(i)*Runiv/Mmix_before*T(i)/V(i);                                  %ideal gas law
    
    %motored cycle values:
    pm(i)=p(i);                                                             
    Tm(i)=T(i);
end
if theta_loop < 4*pi
    theta_loop = theta_loop + dtheta;
else
    theta_loop = 0;
end
end   


%% Plot pV-diagram and TV-diagram
figure(1)
plot(V,p)
xlabel('Volume [m^3]')
ylabel('Druk [Pascal]')
title('pV-diagram')
grid on

figure(2)
plot(V,T)
xlabel('Volume [m^3]')
ylabel('Temperature [K]')
title('TV-diagram')
grid on
