clear all, close all
% Author: Yipeng Li (y9li@ucsd.edu)
global kv1 kv2 k1 k2 k3 k4 Vcmax a b  Ep Ecw T Ti Thold  Qmax t_slope slope_VCV FRC setPvent exPvent exHoldv1 exHoldv2 exHoldp1 exHoldp2 t_instop 
global VD Do13 Dc13 Do Dc Solubility_O2 Solubility_CO2 Vpc Th L KT  KR  L2 h r2 AcceRate Ps Ts Tbody fo fc exHold1 exHold2 KT_co2 KR_co2

% Ventilation parameters£¬pressure (cmH2O)
kv1 = 0.53;             kv2  = 4.41;             k1 = 0.31;               k2 = 0.2;  
k3 = 0.2;               Vcmax = 0.1;             Ep = 10;                 Ecw = 10; 
T = 4;                  Ti = 0.9;                Thold = 0.1;             Qmax = 0.86;     
t_slope = 0.44;         t_instop = 0.87;         exPvent = 0;             slope_VCV = Qmax/t_slope;  
stop_VCV = 180;         stop_PCV = 2*stop_VCV;   h1 = stop_VCV/T;         h2 = stop_PCV/T;
exHoldv1 = (h1-1)*T;    exHoldv2 = h1*T;         exHoldp1 = (h2-1)*T;     exHoldp2 = h2*T;           
k4 = 47.5;                a = 1.05;                 b = 6;                  FRC = 3;      setPvent = 23.683; 

% Gas exchange parameters,pressure (mmHg) 
fo = 0.21;              fc = 0.0004;               Do = 3.5e-4;                    Dc = 7.08e-3;     
Do13 = 1.56e-5;         Dc13 = 3.16e-4;            Solubility_O2 = 1.4e-6;         Solubility_CO2 = 3.3e-5;
Vpc = 0.07;             Th = 2e-3;                 L = 171.2e6;                    KT = 10e3;
KT_co2 = 25.1;          KR_co2 = 29.5;
KR = 3.6e6;             L2 = 164e3;                h = 10^-7.4;                    r2 = 0.12;
Ps = 760;               Ts = 273;                  Tbody = 310;                    VD = 0.05;   
AcceRate = 10^1.9;      T_heart = 60/75;           n = 216;

Pcw0_VCV = 3.7;     Pel0_VCV = 3.7;     Pc0_VCV = 3;
[t_VCV,P_VCV] = ode15s(@(t,y)odeVCV_hold(t,y), [0,stop_VCV], [Pcw0_VCV Pel0_VCV Pc0_VCV]); % VCV model ODE
Pcw_VCV = P_VCV(:,1);   Pel_VCV = P_VCV(:,2);    Pc_VCV = P_VCV(:,3);

Vp_VCV = FRC + Pel_VCV./Ep;
fid_Vc = fopen('Vc_VCV.txt','a');
for i = 1:length(t_VCV)
    Vc_VCV = Vcmax./(1 + exp(-a.*(Pc_VCV(i) - b)));
    if (Vc_VCV  < 0.999*Vcmax) && (Vc_VCV > 0.001*Vcmax)
        Vc_VCV  = Vcmax./(1 + exp(-a.*(Pc_VCV(i)  - b)));
    elseif (Vc_VCV  <= 0.001*Vcmax)
        Vc_VCV  = 0.001*Vcmax;
    else
        Vc_VCV  = 0.999*Vcmax;
    end
    fprintf(fid_Vc,'%d\r\n',Vc_VCV);
end
fclose(fid_Vc);
Vc_VCV = load('Vc_VCV.txt');
delete('Vc_VCV.txt');

Vcw_VCV = Vp_VCV + Vc_VCV;
Rc_VCV = k3.*(Vcmax./Vc_VCV).^2;
Rs_VCV = k4./Vp_VCV;

fid_Q = fopen('Q_VCV.txt','a');
for i = 1:length(t_VCV)
if mod(t_VCV(i),T) <= Ti && (t_VCV(i) <= exHoldv1 || t_VCV(i) > exHoldv2)
    if mod(t_VCV(i),T) <= t_slope
        Qv = slope_VCV.*mod(t_VCV(i),T);
    elseif  mod(t_VCV(i),T) <= t_instop
        Qv = Qmax;
    else
        Qv = (Qmax./(t_instop-Ti)).*mod(t_VCV(i),T) + ( - (Qmax./(t_instop-Ti)).*Ti);
    end
    Q_VCV= Qv;
elseif mod(t_VCV(i),T) <= (Ti + Thold) && (t_VCV(i) <= exHoldv1 || t_VCV(i) > exHoldv2)
    Q_VCV = 0;
elseif t_VCV(i) >= exHoldv1 && t_VCV(i) <= exHoldv2
    Q_VCV = 0;
else
    Q_VCV = (( k1 + kv1 + Rc_VCV(i)) - sqrt(( k1 + kv1 + Rc_VCV(i)).^2 - 4.*( k2 + kv2).*(exPvent - Pc_VCV(i) - Pcw_VCV(i))))./(2.*( k2 + kv2));
end
fprintf(fid_Q,'%d\r\n',Q_VCV);
end
fclose(fid_Q);
Q_VCV = load('Q_VCV.txt');
delete('Q_VCV.txt');

Ru_VCV = k1 + k2.*abs(Q_VCV);
Rv_VCV = kv1 + kv2.*abs(Q_VCV);
Ppa_VCV = Pcw_VCV + Pel_VCV;  
Paw_VCV = (Ru_VCV + Rc_VCV).*Q_VCV + Pc_VCV + Pcw_VCV;
Pvent_VCV = (Rv_VCV + Ru_VCV + Rc_VCV).*Q_VCV + Pc_VCV + Pcw_VCV;
Qp_VCV = (Pc_VCV - Pel_VCV)./Rs_VCV;
Qc_VCV = Q_VCV - Qp_VCV;
Eaw_VCV = Vcmax./(a.*Vc_VCV.*(Vcmax - Vc_VCV));

% VCV gas exchange
t_vent = t_VCV;       Vp_vent = Vp_VCV;     Vc_vent = Vc_VCV;     Q_vent = Q_VCV;
Qp_vent = Qp_VCV;     Qc_vent = Qc_VCV;     Paw_vent = Paw_VCV;   Ppa_vent = Ppa_VCV;
exHold1 = exHoldv1;   exHold2 = exHoldv2;          

t_int = [];             P_int = []; 
for i = 1:n;
    t_begin = (i-1)*T_heart;         t_stop = i*T_heart;
    if i==1
        Pdead_O2_IC = 101;           Pcollap_O2_IC = 100;           PA_O2_IC = 97.5;
        Pdead_CO2_IC = 40;           Pcollap_CO2_IC = 38;           PA_CO2_IC = 37.5;
        Pblood_O2_IC = 40;           Pblood_CO2_IC = 46;            z0 = 46.*Solubility_CO2.*r2./(L2.*h); 
    end
    
    [t_Gas,P_Gas] = ode15s(@(t,y)odeVCV_exchange(t,y,t_vent,Vp_vent,Vc_vent,Q_vent,Qp_vent,Qc_vent,Paw_vent,Ppa_vent),(t_begin:0.01:t_stop),[Pblood_O2_IC Pblood_CO2_IC z0 Pdead_O2_IC Pcollap_O2_IC PA_O2_IC Pdead_CO2_IC Pcollap_CO2_IC PA_CO2_IC]);
    Pdead_O2_IC = P_Gas(length(t_Gas),4);           Pcollap_O2_IC = P_Gas(length(t_Gas),5);           PA_O2_IC = P_Gas(length(t_Gas),6);
    Pdead_CO2_IC = P_Gas(length(t_Gas),7);          Pcollap_CO2_IC = P_Gas(length(t_Gas),8);          PA_CO2_IC = P_Gas(length(t_Gas),9);   
    t_int = cat(1,t_int,t_Gas);
    P_int = cat(1,P_int,P_Gas);
end
t_VCV_Gas = t_int;
P_VCV_Gas = P_int;

% VCV tspan
scopev = find((t_VCV >= (stop_VCV-3*T)) & (t_VCV <= (stop_VCV-2*T)));
Vt_VCV = max(Vcw_VCV(scopev)) - min(Vcw_VCV(scopev));
scopev_control_2T = find((t_VCV >= (stop_VCV - 3*T)) & (t_VCV <= stop_VCV));
scopev_gas_2T = find((t_VCV_Gas >= (stop_VCV - 3*T)) & (t_VCV_Gas <= stop_VCV));
PA_O2_VCV = P_VCV_Gas(:,6);
PA_CO2_VCV = P_VCV_Gas(:,9);

% PCV Model
Pcw0_PCV = min(Pcw_VCV(scopev));     Pel0_PCV = min(Pel_VCV(scopev));     Pc0_PCV = min(Pc_VCV(scopev));
[t_PCV,P_PCV] = ode15s(@(t,y)odePCV_hold(t,y), [stop_VCV, stop_PCV], [Pcw0_PCV Pel0_PCV Pc0_PCV]);
Pcw_PCV = P_PCV(:,1);   Pel_PCV = P_PCV(:,2);    Pc_PCV = P_PCV(:,3);

Vp_PCV = FRC + Pel_PCV./Ep;
fid_Vc = fopen('Vc_PCV.txt','a');
for i = 1:length(t_PCV)
    Vc_PCV = Vcmax./(1 + exp(-a.*(Pc_PCV(i)  - b)));
    if (Vc_PCV  < 0.999*Vcmax) && (Vc_PCV > 0.001*Vcmax)
        Vc_PCV  = Vcmax./(1 + exp(-a.*(Pc_PCV(i) - b)));
    elseif (Vc_PCV  <= 0.001*Vcmax)
        Vc_PCV  = 0.001*Vcmax;
    else
        Vc_PCV  = 0.999*Vcmax;
    end
    fprintf(fid_Vc,'%d\r\n',Vc_PCV);
end
fclose(fid_Vc);
Vc_PCV = load('Vc_PCV.txt');
delete('Vc_PCV.txt');

Vcw_PCV = Vp_PCV + Vc_PCV;
Rc_PCV = k3.*(Vcmax./Vc_PCV).^2;
Rs_PCV = k4./Vp_PCV;

fid_Q = fopen('Q_PCV.txt','a');
for i = 1:length(t_PCV)
if mod(t_PCV(i),T) <= Ti && (t_PCV(i) <= exHoldp1 || t_PCV(i) > exHoldp2)
    inPvent = setPvent; 
    Q_PCV= ( - ( k1 + kv1 + Rc_PCV(i)) + sqrt(( k1 + kv1 + Rc_PCV(i)).^2 + 4.*( k2 + kv2).*(inPvent - Pc_PCV(i) - Pcw_PCV(i))))./(2.*( k2 + kv2));
elseif mod(t_PCV(i),T) <= (Ti + Thold) && (t_PCV(i) <= exHoldp1 || t_PCV(i) > exHoldp2)
    Q_PCV = 0;
elseif t_PCV(i) >= exHoldp1 && t_PCV(i) <= exHoldp2
    Q_PCV = 0;
else
    Q_PCV = (( k1 + kv1 + Rc_PCV(i)) - sqrt(( k1 + kv1 + Rc_PCV(i)).^2 - 4.*( k2 + kv2).*(exPvent - Pc_PCV(i) - Pcw_PCV(i))))./(2.*( k2 + kv2));
end
fprintf(fid_Q,'%d\r\n',Q_PCV);
end
fclose(fid_Q);
Q_PCV = load('Q_PCV.txt');
delete('Q_PCV.txt');

Ru_PCV = k1 + k2.*abs(Q_PCV);
Rv_PCV = kv1 + kv2.*abs(Q_PCV);
Ppa_PCV = Pcw_PCV + Pel_PCV;  
Paw_PCV = (Ru_PCV + Rc_PCV).*Q_PCV + Pc_PCV + Pcw_PCV;
Pvent_PCV = (Rv_PCV + Ru_PCV + Rc_PCV).*Q_PCV + Pc_PCV + Pcw_PCV;
Qp_PCV = (Pc_PCV - Pel_PCV)./Rs_PCV;
Qc_PCV = Q_PCV - Qp_PCV;
Eaw_PCV = Vcmax./(a.*Vc_PCV.*(Vcmax - Vc_PCV));

% PCV Gas Exchange
t_vent = t_PCV;       Vp_vent = Vp_PCV;     Vc_vent = Vc_PCV;     Q_vent = Q_PCV;
Qp_vent = Qp_PCV;     Qc_vent = Qc_PCV;     Paw_vent = Paw_PCV;   Ppa_vent = Ppa_PCV;
exHold1 = exHoldp1;   exHold2 = exHoldp2;          

t_int = [];             P_int = []; 
for i = 1:n;
    t_begin = (i-1)*T_heart + stop_VCV;         t_stop = i*T_heart + stop_VCV;
    if i==1
        Pdead_O2_IC = 101;           Pcollap_O2_IC = 100;           PA_O2_IC = 97.5;
        Pdead_CO2_IC = 40;           Pcollap_CO2_IC = 38;           PA_CO2_IC = 37.5;
        Pblood_O2_IC = 40;           Pblood_CO2_IC = 46;            z0 = 46.*Solubility_CO2.*r2./(L2.*h); 
    end
    
    [t_Gas,P_Gas] = ode15s(@(t,y)odePCV_exchange(t,y,t_vent,Vp_vent,Vc_vent,Q_vent,Qp_vent,Qc_vent,Paw_vent,Ppa_vent),(t_begin:0.01:t_stop),[Pblood_O2_IC Pblood_CO2_IC z0 Pdead_O2_IC Pcollap_O2_IC PA_O2_IC Pdead_CO2_IC Pcollap_CO2_IC PA_CO2_IC]);
    Pdead_O2_IC = P_Gas(length(t_Gas),4);           Pcollap_O2_IC = P_Gas(length(t_Gas),5);           PA_O2_IC = P_Gas(length(t_Gas),6);
    Pdead_CO2_IC = P_Gas(length(t_Gas),7);          Pcollap_CO2_IC = P_Gas(length(t_Gas),8);          PA_CO2_IC = P_Gas(length(t_Gas),9);   
    t_int = cat(1,t_int,t_Gas);
    P_int = cat(1,P_int,P_Gas);
end
t_PCV_Gas = t_int;
P_PCV_Gas = P_int;

% PCV tspan
scopep = find((t_PCV >= (stop_PCV - 3*T)) & (t_PCV <= (stop_PCV - 2*T)));
Vt_PCV = max(Vcw_PCV(scopep)) - min(Vcw_PCV(scopep));
scopep_control_2T = find((t_PCV >= (stop_PCV - 3*T)) & (t_PCV <= stop_PCV));
scopep_gas_2T = find((t_PCV_Gas >= (stop_PCV - 3*T)) & (t_PCV_Gas <= stop_PCV));
PA_O2_PCV = P_PCV_Gas(:,6);
PA_CO2_PCV = P_PCV_Gas(:,9);


figure;
% Subplot 1: Po
subplot(2,2,1);
plot(t_VCV_Gas(scopev_gas_2T), P_VCV_Gas(scopev_gas_2T, 1), 'r--'); % VCV with dashed orange line
hold on;
plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV, P_PCV_Gas(scopep_gas_2T, 1), 'b'); % PCV with solid blue line
xlabel('Time (s)');
ylabel('P_o (mmHg)');
grid on;
xlim([168 172]);
xticks(168:172);
xticklabels({'0', '1', '2', '3', '4'});
text(0.1, 0.9, 'a', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 2: Pc
subplot(2,2,2);
plot(t_VCV_Gas(scopev_gas_2T), P_VCV_Gas(scopev_gas_2T, 2), 'r--'); % VCV with dashed orange line
hold on;
plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV, P_PCV_Gas(scopep_gas_2T, 2), 'b'); % PCV with solid blue line
xlabel('Time (s)');
ylabel('P_c (mmHg)');
legend('VCV', 'PCV');
grid on;
xlim([168 172]);
ylim([34 50]);
xticks(168:172);
xticklabels({'0', '1', '2', '3', '4'});
text(0.1, 0.9, 'b', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 7: Pao
subplot(2,2,3);
plot(t_VCV_Gas(scopev_gas_2T), P_VCV_Gas(scopev_gas_2T, 6), 'r--'); % VCV with dashed orange line
hold on;
plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV, P_PCV_Gas(scopep_gas_2T, 6), 'b'); % PCV with solid blue line
xlabel('Time (s)');
ylabel('P_a_o (mmHg)');
grid on;
xlim([168 172]);
ylim([126 134]);
xticks(168:172);
xticklabels({'0', '1', '2', '3', '4'});
text(0.1, 0.9, 'c', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 8: Pac
subplot(2,2,4);
plot(t_VCV_Gas(scopev_gas_2T), P_VCV_Gas(scopev_gas_2T, 9), 'r--'); % VCV with dashed orange line
hold on;
plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV, P_PCV_Gas(scopep_gas_2T, 9), 'b'); % PCV with solid blue line
xlabel('Time (s)');
ylabel('P_a_c (mmHg)');
grid on;
xlim([168 172]);
ylim([34 40]);
xticks(168:172);
xticklabels({'0', '1', '2', '3', '4'});
text(0.1, 0.9, 'd', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold');

%Extract PA O2 (column 6) and PA CO2 (column 9)
PCV_PA_O2_NoCO2Binding = P_PCV_Gas(:, 6);
PCV_PA_CO2_NoCO2Binding = P_PCV_Gas(:, 9);
VCV_PA_O2_NoCO2Binding = P_VCV_Gas(:, 6);
VCV_PA_CO2_NoCO2Binding = P_VCV_Gas(:, 9);

% PCV_PA_O2_withCO2Binding = P_PCV_Gas(:, 6);
% PCV_PA_CO2_withCO2Binding = P_PCV_Gas(:, 9);
% VCV_PA_O2_withCO2Binding = P_VCV_Gas(:, 6);
% VCV_PA_CO2_withCO2Binding = P_VCV_Gas(:, 9);

% Save data
%save('PCV_PA_withCO2Binding.mat', 'PCV_PA_O2_withCO2Binding', 'PCV_PA_CO2_withCO2Binding');
%save('VCV_PA_withCO2Binding.mat', 'VCV_PA_O2_withCO2Binding', 'VCV_PA_CO2_withCO2Binding');

% All 13 figures
% figure,subplot(3,1,1),plot(t_VCV_Gas,P_VCV_Gas(:,1)),hold on, plot(t_VCV_Gas,P_VCV_Gas(:,2),'r'),legend('Pblood O2','Pblood CO2')
%        subplot(3,1,2),plot(t_VCV_Gas,P_VCV_Gas(:,4)),hold on, plot(t_VCV_Gas,P_VCV_Gas(:,5),'r'),plot(t_VCV_Gas,P_VCV_Gas(:,6),'g'),legend('Pdead O2','Pcollap O2','PA O2','Paw')
%                       plot(t_VCV,(Paw_VCV./1.36 + Ps).*fo,'m')
%        subplot(3,1,3),plot(t_VCV_Gas,P_VCV_Gas(:,7)),hold on, plot(t_VCV_Gas,P_VCV_Gas(:,8),'r'),plot(t_VCV_Gas,P_VCV_Gas(:,9),'g'),legend('Pdead CO2','Pcollap CO2','PA CO2')
% figure,subplot(2,1,1),plot(t_VCV_Gas,P_VCV_Gas(:,6)),hold on, plot(t_VCV_Gas,P_VCV_Gas(:,1),'g'),legend('PA O2','Pblood O2'),
%        subplot(2,1,2),plot(t_VCV_Gas,P_VCV_Gas(:,9)),hold on, plot(t_VCV_Gas,P_VCV_Gas(:,2),'g'),legend('PA CO2','Pblood CO2'),
% figure, subplot(3,1,1),plot(t_PCV_Gas,P_PCV_Gas(:,1)),hold on, plot(t_PCV_Gas,P_PCV_Gas(:,2),'r'),legend('Pblood O2','Pblood CO2')
%         subplot(3,1,2),plot(t_PCV_Gas,P_PCV_Gas(:,4)),hold on, plot(t_PCV_Gas,P_PCV_Gas(:,5),'r'),plot(t_PCV_Gas,P_PCV_Gas(:,6),'g'),legend('Pdead O2','Pcollap O2','PA O2')
%                        plot(t_PCV,(Paw_PCV./1.36 + Ps).*fo,'m')
%         subplot(3,1,3),plot(t_PCV_Gas,P_PCV_Gas(:,7)),hold on, plot(t_PCV_Gas,P_PCV_Gas(:,8),'r'),plot(t_PCV_Gas,P_PCV_Gas(:,9),'g'),legend('Pdead CO2','Pcollap CO2','PA CO2')
% figure,subplot(2,1,1),plot(t_PCV_Gas,P_PCV_Gas(:,6)),hold on, plot(t_PCV_Gas,P_PCV_Gas(:,1),'g'),legend('PA O2','Pblood O2'),
%        subplot(2,1,2),plot(t_PCV_Gas,P_PCV_Gas(:,9)),hold on, plot(t_PCV_Gas,P_PCV_Gas(:,2),'g'),legend('PA CO2','Pblood CO2'),
% 
% figure,subplot(4,1,1),plot(t_VCV,Paw_VCV,'k'),hold on, plot(t_VCV,Ppa_VCV,'g'), plot(t_VCV,Pcw_VCV,'k-.'),
%                       plot(t_PCV,Paw_PCV,'k'),hold on, plot(t_PCV,Ppa_PCV,'g'), plot(t_PCV,Pcw_PCV,'k-.'),
%                       title('Pressure'),legend('Paw','Ppa','Pcw'),grid on;
%        subplot(4,1,2),plot(t_VCV,Vcw_VCV,'k'),hold on, plot(t_VCV,Vp_VCV,'g'), plot(t_VCV,Vc_VCV,'b'),
%                       plot(t_PCV,Vcw_PCV,'k'),hold on, plot(t_PCV,Vp_PCV,'g'), plot(t_PCV,Vc_PCV,'b'), 
%                       title('volume'),legend('Vcw','Vp','Vc'),grid on;      
%        subplot(4,1,3),plot(t_VCV,Q_VCV,'k'),hold on, plot(t_VCV,Qp_VCV,'g'), plot(t_VCV,Qc_VCV,'b'),
%                       plot(t_PCV,Q_PCV,'k'),hold on, plot(t_PCV,Qp_PCV,'g'), plot(t_PCV,Qc_PCV,'b'),
%                       title('flow'),legend('Q','Qp','Qc'),grid on;
%        subplot(4,1,4),plot(t_VCV,Q_VCV,'k'),hold on, plot(t_PCV,Q_PCV,'k'), title('flow'),legend('Q'),grid on;
% figure,subplot(2,1,1),plot(t_VCV,Pcw_VCV+Pc_VCV,'k'), hold on, plot(t_VCV,Ppa_VCV),
%                       title('Pressure VCV'),legend('Pcw+Pc','Ppa'),grid on;
%        subplot(2,1,2),plot(t_PCV,Pcw_PCV+Pc_PCV,'k'), hold on, plot(t_PCV,Ppa_PCV), 
%                       title('Pressure PCV'),legend('Pcw+Pc','Ppa'),grid on;
% figure,subplot(2,2,1),plot(t_VCV,Rs_VCV), hold on, plot(t_PCV,Rs_PCV), 
%                       title('Rs'),legend('Rs'),grid on;
%        subplot(2,2,2),plot(t_VCV,Pel_VCV), hold on, plot(t_PCV,Pel_PCV),  
%                       title('Pel'),legend('Pel'),grid on;
%        subplot(2,2,3),plot(t_VCV,Pc_VCV), hold on,  plot(t_PCV,Pc_PCV), 
%                       title('Pc'),legend('Pc'),grid on;
%        subplot(2,2,4),plot(t_VCV,Rc_VCV), hold on, plot(t_PCV,Rc_PCV),
%                       title('Rc'),legend('Rc'),grid on; 
% figure,plot(t_PCV, Pvent_PCV)
% 
% 
% figure,subplot(4,2,1), plot(t_VCV(scopev_control_2T),Pvent_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Pvent_PCV(scopep_control_2T),'g'), title('Pvent'), grid on
%        subplot(4,2,2), plot(t_VCV(scopev_control_2T),Paw_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Paw_PCV(scopep_control_2T),'g'),   title('Paw'),   grid on
%        subplot(4,2,3), plot(t_VCV(scopev_control_2T),Pc_VCV(scopev_control_2T) + Pcw_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Pc_PCV(scopep_control_2T) + Pcw_PCV(scopep_control_2T),'g'),   title('Pc + Pcw'),   grid on                       
%        subplot(4,2,4), plot(t_VCV(scopev_control_2T),Ppa_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Ppa_PCV(scopep_control_2T),'g'),   title('Ppa'),   grid on
%        subplot(4,2,5), plot(t_VCV(scopev_control_2T),Pcw_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Pcw_PCV(scopep_control_2T),'g'),   title('Pcw'),   grid on
%        subplot(4,2,6), plot(t_VCV(scopev_control_2T),Pc_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Pc_PCV(scopep_control_2T),'g'),    title('Pc'),   grid on
%        subplot(4,2,7), plot(t_VCV(scopev_control_2T),Pel_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Pel_PCV(scopep_control_2T),'g'),   title('Pel'),   grid on
%        subplot(4,2,8), plot(Vc_VCV(scopev_control_2T),Eaw_VCV(scopev_control_2T)), hold on, 
%                        plot(Vc_PCV(scopep_control_2T),Eaw_PCV(scopep_control_2T),'g'),     title('Vc - Eaw'),   grid on
% figure,subplot(2,4,1), plot(t_VCV(scopev_control_2T),Rv_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Rv_PCV(scopep_control_2T),'g'),   title('Rv'),   grid on
%        subplot(2,4,2), plot(t_VCV(scopev_control_2T),Ru_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Ru_PCV(scopep_control_2T),'g'),   title('Ru'),   grid on
%        subplot(2,4,3), plot(t_VCV(scopev_control_2T),Rc_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Rc_PCV(scopep_control_2T),'g'),   title('Rc'),   grid on                       
%        subplot(2,4,4), plot(t_VCV(scopev_control_2T),Rs_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Rs_PCV(scopep_control_2T),'g'),   title('Rs'),   grid on                       
%        subplot(2,4,5), plot(t_VCV(scopev_control_2T),Rv_VCV(scopev_control_2T).*Q_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Rv_PCV(scopep_control_2T).*Q_PCV(scopep_control_2T),'g'),   title('Rv*Q'),   grid on
%        subplot(2,4,6), plot(t_VCV(scopev_control_2T),Ru_VCV(scopev_control_2T).*Q_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Ru_PCV(scopep_control_2T).*Q_PCV(scopep_control_2T),'g'),   title('Ru*Q'),   grid on                       
%        subplot(2,4,7), plot(t_VCV(scopev_control_2T),Rc_VCV(scopev_control_2T).*Q_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Rc_PCV(scopep_control_2T).*Q_PCV(scopep_control_2T),'g'),   title('Rc*Q'),   grid on                 
%        subplot(2,4,8), plot(t_VCV(scopev_control_2T),Rs_VCV(scopev_control_2T).*Qp_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Rs_PCV(scopep_control_2T).*Qp_PCV(scopep_control_2T),'g'),   title('Rs*Qp'),   grid on                       
% figure,subplot(2,3,1), plot(t_VCV(scopev_control_2T),Q_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Q_PCV(scopep_control_2T),'g'),   title('Q'),   grid on     
%        subplot(2,3,2), plot(t_VCV(scopev_control_2T),Qp_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Qp_PCV(scopep_control_2T),'g'),   title('Qp'),   grid on                  
%        subplot(2,3,3), plot(t_VCV(scopev_control_2T),Qc_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Qc_PCV(scopep_control_2T),'g'),   title('Qc'),   grid on                  
%        subplot(2,3,4), plot(t_VCV(scopev_control_2T),Vcw_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Vcw_PCV(scopep_control_2T),'g'),   title('Vcw'),   grid on      
%        subplot(2,3,5), plot(t_VCV(scopev_control_2T),Vp_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Vp_PCV(scopep_control_2T),'g'),   title('Vp'),   grid on      
%        subplot(2,3,6), plot(t_VCV(scopev_control_2T),Vc_VCV(scopev_control_2T)), hold on, 
%                        plot(t_PCV(scopep_control_2T)-stop_VCV,Vc_PCV(scopep_control_2T),'g'),   title('Vc'),   grid on      
% figure,subplot(4,2,1), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),1)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),1),'g-.'),title('Pblood O2'), grid on                    
%        subplot(4,2,2), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),2)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),2),'g-.'),title('Pblood CO2'), grid on   
%        subplot(4,2,3), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),4)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),4),'g-.'),title('Pdead O2'), grid on
%        subplot(4,2,4), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),7)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),7),'g-.'),title('Pdead CO2'), grid on 
%        subplot(4,2,5), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),5)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),5),'g-.'),title('Pcollapse O2'), grid on
%        subplot(4,2,6), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),8)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),8),'g-.'),title('Pcollapse CO2'), grid on
%        subplot(4,2,7), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),6)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),6),'g-.'),title('PA O2'), grid on
%        subplot(4,2,8), plot(t_VCV_Gas(scopev_gas_2T),P_VCV_Gas((scopev_gas_2T),9)), hold on, 
%                        plot(t_PCV_Gas(scopep_gas_2T)-stop_VCV,P_PCV_Gas((scopep_gas_2T),9),'g-.'),title('PA CO2'), grid on  


