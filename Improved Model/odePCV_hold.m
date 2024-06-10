function dy = odePCV_hold(t,y)
global  kv1 kv2 k1 k2 k3 k4 Vcmax a b  Ep Ecw T Ti Thold FRC setPvent exPvent exHoldp1 exHoldp2 

dy = zeros(3,1);
Vp = FRC + y(2)./Ep;
Vc = Vcmax./(1 + exp(-a.*(y(3) - b)));
if (Vc < 0.999*Vcmax) && (Vc> 0.001*Vcmax)
    Vc = Vcmax./(1 + exp(-a.*(y(3) - b)));
elseif (Vc <= 0.001*Vcmax)
    Vc = 0.001*Vcmax;
else
    Vc = 0.999*Vcmax;
end
Rc = k3.*(Vcmax./Vc).^2;
Rs = k4./Vp;
Eaw = Vcmax./(a.*Vc.*(Vcmax - Vc));

if mod(t,T) <= Ti  && (t <= exHoldp1 || t > exHoldp2)
    inPvent = setPvent; 
    dy(1) = Ecw.*(( - ( k1 + kv1 + Rc) + sqrt(( k1 + kv1 + Rc).^2 + 4.*( k2 + kv2).*(inPvent - y(1) - y(3))))./(2.*( k2 + kv2)));
    dy(2) = Ep.*(y(3) - y(2))./Rs;    % Pel
    dy(3) = Eaw.*(((- ( k1 + kv1 + Rc) + sqrt(( k1 + kv1 + Rc).^2 + 4.*( k2 + kv2).*(inPvent - y(1) - y(3))))./(2.*( k2 + kv2))) - (y(3) - y(2))./Rs);   % Pc
elseif mod(t,T) <= (Ti+Thold) && (t <= exHoldp1 || t > exHoldp2)
    dy(1) = 0;
    dy(2) =   Ep.*(y(3) - y(2))./Rs;         % Pel
    dy(3) = -Eaw.*(y(3) - y(2))./Rs;         % Pc
elseif (t >= exHoldp1 && t <= exHoldp2)
    dy(1) = 0;
    dy(2) =   Ep.*(y(3) - y(2))./Rs;         % Pel
    dy(3) = -Eaw.*(y(3) - y(2))./Rs;         % Pc
else
    dy(1) = Ecw.*(( ( k1 + kv1 + Rc) - sqrt(( k1 + kv1 + Rc).^2 - 4.*( k2 + kv2).*(exPvent - y(1) - y(3))))./(2.*( k2 + kv2)));
    dy(2) = Ep.*(y(3) - y(2))./Rs;            % Pel
    dy(3) = Eaw.*(((( k1 + kv1 + Rc) - sqrt(( k1 + kv1 + Rc).^2 - 4.*( k2 + kv2).*(exPvent - y(1) - y(3))))./(2.*( k2 + kv2))) - (y(3) - y(2))./Rs);   % Pc
end