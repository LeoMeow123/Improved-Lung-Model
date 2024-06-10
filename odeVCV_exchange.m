function dy = odeVCV_exchange(t,y,t_vent,Vp_vent,Vc_vent,Q_vent,Qp_vent,Qc_vent,Paw_vent,Ppa_vent)
global VD Do13 Dc13 Do Dc Solubility_O2 Solubility_CO2 Vpc Th L KT KR L2 h r2 AcceRate Ps Ts Tbody fo fc exHold1 exHold2  T Ti Thold

dy = zeros(9,1);

Vp = interp1(t_vent,Vp_vent,t);
Vc = interp1(t_vent,Vc_vent,t);
Q = interp1(t_vent,Q_vent,t);
Qp = interp1(t_vent,Qp_vent,t);
Qc = interp1(t_vent,Qc_vent,t);
Paw = interp1(t_vent,Paw_vent,t);
Ppa = interp1(t_vent,Ppa_vent,t);

Paw = Paw./1.36 + Ps;
Ppa = Ppa./1.36 + Ps;

Paw_O2 = Paw.*fo;
Paw_CO2 = Paw.*fc;

dVc_dt = Qc;
dVpa_dt = Qp - (Ps.*Tbody./(Ts.*Ppa))*(Do.*(y(6)-y(1)) + Dc.*(y(9)-y(2)));

Af = KT.*Solubility_O2.*y(1);
Bf = KR.*Solubility_O2.*y(1);
df_dPo = (KR*Solubility_O2*(Bf + 1)^3 + 3*KR*Solubility_O2*Bf*(Bf + 1)^2 + KT*L*Solubility_O2*(Af + 1)^3 + 3*KT*L*Solubility_O2*Af*(Af + 1)^2)/((Bf + 1)^4 + L*(Af + 1)^4)...
          - ((Bf*(Bf + 1)^3 + L*Af*(Af + 1)^3)*(4*KR*Solubility_O2*(Bf + 1)^3 + 4*KT*L*Solubility_O2*(Af + 1)^3))/((Bf + 1)^4 + L*(Af + 1)^4)^2;

dy(1) = Do13.*(1 + 4.*Th.*df_dPo./Solubility_O2)^(-1).*(y(6) - y(1))./(Solubility_O2.*Vpc);                          % O2 partial pressure in blood
dy(2) = Dc13.*(y(9) - y(2))./(Solubility_CO2.*Vpc) + AcceRate.*L2.*h.*y(3)./Solubility_CO2 - AcceRate.*r2.*y(2);     % CO2 partial pressure in blood
dy(3) = AcceRate.*r2.*Solubility_CO2.*y(2) - AcceRate.*L2.*h.*y(3);        % z

if mod(t,T) <= (Ti + Thold) && (t <= exHold1 || t > exHold2) 
    dy(4) = Q.*(Paw_O2 - y(4))./VD;                                                    % O2 partial pressure in dead space
    dy(5) = (Q.*y(4) - Qp.*y(5) - y(5).*dVc_dt)./Vc;                                   % O2 partial pressure in collapsible airways
    dy(6) = (Qp.*y(5) - y(6).*dVpa_dt - Ps.*Tbody.*Do.*(y(6) - y(1))./Ts)./Vp;         % O2 partial pressure in alveoli
    
    dy(7) = Q.*(Paw_CO2 - y(7))./VD;                                                   % CO2 partial pressure in dead space
    dy(8) = (Q.*y(7) - Qp.*y(8) - y(8).*dVc_dt)./Vc;                                   % CO2 partial pressure in collapsible airways
    dy(9) = (Qp.*y(8) - y(9).*dVpa_dt - Ps.*Tbody.*Dc.*(y(9) - y(2))./Ts)./Vp;         % CO2 partial pressure in alveoli
else
    dy(4) = Q.*(y(4) - y(5))./VD;                                                      %dPD_O2,单位是mmHg
    dy(5) = (Q.*y(5) - Qp.*y(6) - y(5).*dVc_dt)./Vc;                                   %dPc_O2,单位是mmHg
    dy(6) = (Qp.*y(6) - y(6).*dVpa_dt - Ps.*Tbody.*Do.*(y(6) - y(1))./Ts)./Vp;         %dPA_O2,单位是mmHg
    
    dy(7) = Q.*(y(7) - y(8))./VD;                                                      %dPD_CO2,单位是mmHg
    dy(8) = (Q.*y(8) - Qp.*y(9) - y(8).*dVc_dt)./Vc;                                   %dPc_CO2,单位是mmHg
    dy(9) = (Qp.*y(9) - y(9).*dVpa_dt - Ps.*Tbody.*Dc.*(y(9) - y(2))./Ts)./Vp;         %dPA_CO2,单位是mmHg
end

% fid_odet = fopen('odet.txt','a');   fprintf(fid_odet,'%d\r\n',t);      fclose(fid_odet);
% fid_Vp = fopen('Vp.txt','a');   fprintf(fid_Vp,'%d\r\n',Vp);      fclose(fid_Vp);
% fid_Vc = fopen('Vc.txt','a');   fprintf(fid_Vc,'%d\r\n',Vc);      fclose(fid_Vc);
% fid_Q = fopen('Q.txt','a');   fprintf(fid_Q,'%d\r\n',Q);      fclose(fid_Q);
% fid_Qp = fopen('Qp.txt','a');   fprintf(fid_Qp,'%d\r\n',Qp);      fclose(fid_Qp);
% fid_Qc = fopen('Qc.txt','a');   fprintf(fid_Qc,'%d\r\n',Qc);      fclose(fid_Qc);
% fid_Paw = fopen('Paw.txt','a');   fprintf(fid_Paw,'%d\r\n',Paw);      fclose(fid_Paw);
% fid_Ppa = fopen('Ppa.txt','a');   fprintf(fid_Ppa,'%d\r\n',Ppa);      fclose(fid_Ppa);
% fid_PawO2 = fopen('PawO2.txt','a');   fprintf(fid_PawO2,'%d\r\n',Paw_O2);      fclose(fid_PawO2);
% fid_PawCO2 = fopen('PawCO2.txt','a');   fprintf(fid_PawCO2,'%d\r\n',Paw_CO2);      fclose(fid_PawCO2);
















% elseif mod(t,T) <= (Ti+Thold) && (t <= exHoldv1 || t > exHoldv2)
%     dy(1) = 0;
%     dy(2) =   Ep.*(y(3) - y(2))./Rs;         % Pel
%     dy(3) = -Eaw.*(y(3) - y(2))./Rs;         % Pc
%     Q = 0;
%     Qp = (y(3) - y(2))./Rs;
%     Ppa = y(1) + y(2);
%     y(4) = Q - Qp;
%     y(5) = Qp - Ps.*Tbody.*(Do.*(y(9) - y(4)) + Dc.*(y(12) - y(5)))./(Ppa.*Ts);
%     Ru = k1 + k2.*abs(Q);
%     Paw = (Ru + Rc).*Q + y(1) + y(3);
%     Paw_O2 = Paw.*0.21 + 149;
%     Paw_CO2 = Paw.*0.0004 + 0.3;
%     dy(6) = Do.*(1 + 4.*Th.*df_dPo./Solubility_O2)^(-1).*(y(11) - y(6))./(Solubility_O2.*Vpc);                         % O2 partial pressure in blood
%     dy(7) = Dc.*(y(14) - y(7))./(Solubility_CO2.*Vpc) + AcceRate.*L2.*h.*y(8)./Solubility_CO2 - AcceRate.*r2.*y(7);    % CO2 partial pressure in blood
%     dy(8) = AcceRate.*r2.*Solubility_CO2.*y(7) - AcceRate.*L2.*h.*y(8);                % z
%     dy(9) = Q.*(Paw_O2 - y(9))./VD;                                                    % O2 partial pressure in dead space
%     dy(10) = (Q.*y(9) - Qp.*y(10) - y(10).*y(4))./Vc;                                  % O2 partial pressure in collapsible airways
%     dy(11) = (Qp.*y(10) - y(11).*y(5) - Ps.*Tbody.*Do.*(y(11)- y(6))./Ts)./Vp;         % O2 partial pressure in alveoli
%     dy(12) = Q.*(Paw_CO2 - y(12))./VD;                                                 % CO2 partial pressure in dead space
%     dy(13) = (Q.*y(12) - Qp.*y(13) - y(13).*y(4))./Vc;                                 % CO2 partial pressure in collapsible airways
%     dy(14) = (Qp.*y(13) - y(14).*y(5) - Ps.*Tbody.*Dc.*(y(14)- y(7))./Ts)./Vp;         % CO2 partial pressure in alveoli
% elseif (t >= exHoldv1 && t <= exHoldv2)
%     dy(1) = 0;
%     dy(2) =   Ep.*(y(3) - y(2))./Rs;         % Pel
%     dy(3) = -Eaw.*(y(3) - y(2))./Rs;         % Pc
%     Q = 0;
%     Qp = (y(3) - y(2))./Rs;
%     Ppa = y(1) + y(2);
%     y(4) = -Q - Qp;
%     y(5) = -Qp - Ps.*Tbody.*(Do.*(y(9) - y(4)) + Dc.*(y(12) - y(5)))./(Ppa.*Ts);
%     dy(6) = Do.*(1 + 4.*Th.*df_dPo./Solubility_O2)^(-1).*(y(11) - y(6))./(Solubility_O2.*Vpc);                         % O2 partial pressure in blood
%     dy(7) = Dc.*(y(14) - y(7))./(Solubility_CO2.*Vpc) + AcceRate.*L2.*h.*y(8)./Solubility_CO2 - AcceRate.*r2.*y(7);    % CO2 partial pressure in blood
%     dy(8) = AcceRate.*r2.*Solubility_CO2.*y(7) - AcceRate.*L2.*h.*y(8);                % z
%     dy(9) = Q.*(y(9) - y(10))./VD;                                                     %dPD_O2,单位是mmHg
%     dy(10) = (Q.*y(10) - Qp.*y(11) - y(10).*y(4))./Vc;                                 %dPc_O2,单位是mmHg
%     dy(11) = (Qp.*y(11) - y(11).*y(5) - Ps.*Tbody.*Do.*(y(11)-y(6))./Ts)./Vp;          %dPA_O2,单位是mmHg
%     dy(12) = Q.*(y(12) - y(13))./VD;                                                   %dPD_CO2,单位是mmHg
%     dy(13) = (Q.*y(13) - Qp.*y(14) - y(13).*y(4))./Vc;                                 %dPc_CO2,单位是mmHg
%     dy(14) = (Qp.*y(14) - y(14).*y(5) - Ps.*Tbody.*Dc.*(y(14)-y(7))./Ts)./Vp;          %dPA_CO2,单位是mmHg

% 
% fid_Pao = fopen('Pao_VCV.txt','a');
% fprintf(fid_Pao,'%d\r\n',y(9));
% fclose(fid_Pao);
% fid_Po = fopen('Po_VCV.txt','a');
% fprintf(fid_Po,'%d\r\n',y(4));
% fclose(fid_Po);
% fid_df_dPo = fopen('df_dPo_VCV.txt','a');
% fprintf(fid_df_dPo,'%d\r\n',df_dPo);
% fclose(fid_df_dPo);
% fid_f = fopen('f_VCV.txt','a');
% fprintf(fid_f,'%d\r\n',f);
% fclose(fid_f);

