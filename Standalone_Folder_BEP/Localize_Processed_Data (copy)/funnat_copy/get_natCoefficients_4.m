% This function fits the gamma coefficients given fitted aberration bead
% data where this first 5 aberrations are fitted (n + |m| <= 4).
% TODO: make universal script for arbittrary number of aberrations.

function [RAstig3,RComa3,RCurv6] = get_natCoefficients_4(xn,yn,zers,params)

Ncfg = params.Ncfg;
tol = 1e-6;

xn2 = xn.^2;
yn2 = yn.^2;

% B = ones(Ncfg,1);
Z5 = zers(:,1);
Z6 = zers(:,2);
Z7 = zers(:,3);
Z8 = zers(:,4);
Z22 = zers(:,12);

MAstig3 = zeros(2*Ncfg,14);
MComa3 = zeros(2*Ncfg,10);
MCurv6 = zeros(Ncfg,4);
for nn = 1:Ncfg
    mm = 2*nn;
    
    x = [xn(nn) xn2(nn) xn3(nn) xn4(nn)];
    y = [yn(nn) yn2(nn) yn3(nn) yn4(nn)];
    
    Z56(mm-1:mm,1) = [Z5(nn); Z6(nn)];
    Z78(mm-1:mm,1) = [Z7(nn); Z8(nn)];
    
    %LS Matrix for Astig3 (third order astigmatism)
    M1 = [2*(x(3)*y(1)+x(1)*y(3)) x(2)*y(1)+y(3) x(3)+x(1)*y(2) x(2)+y(2) 0 ...
        x(1)*y(2)-x(3) 2*x(1)*y(2) y(1) -x(1) 2*x(1)*y(1) ...
        x(1) y(1) 1 0];
    M2 = [y(4)-x(4) -x(3)-x(1)*y(2) x(2)*y(1)+y(3) 0 x(2)+y(2) ...
        -2*x(2)*y(1) y(3)-x(2)*y(1) x(1) y(1) y(2)-x(2) ...
        y(1) -x(1) 0 1];
    MAstig3(mm-1:mm,:) = [M1;M2];
    
    %LS Matrix for Coma3 (field cubed coma: W331)
    M1 = [x(3)+x(1)*y(2) x(2) x(1)*y(1) x(1) x(2)+y(2) 0 y(1) x(1) 1 0];
    M2 = [y(3)+x(2)*y(1) x(1)*y(1) y(2) y(1) 0 x(2)+y(2) x(1) -y(1) 0 1];
    MComa3(mm-1:mm,:) = [M1;M2];
    
    % (Z60) LS Matrix for 5th order spherical
    MCurv6(nn,:) = [x(2)+y(2) x(1) y(1) 1];
  
end

RAstig3 = pinv(MAstig3,tol)*Z56;
fAstig3 = MAstig3*RAstig3;

RComa3 = pinv(MComa3,tol)*Z78;
fComa3 = MComa3*RComa3;

RCurv6 = pinv(MCurv6,tol)*Z22;
fCurv6 = MCurv6*RCurv6;

% R^2 values
[r25,rmse5] = rsquare(Z5,fAstig3(1:2:end));
[r26,rmse6] = rsquare(Z6,fAstig3(2:2:end));
[r27,rmse7] = rsquare(Z7,fComa3(1:2:end));
[r28,rmse8] = rsquare(Z8,fComa3(2:2:end));
[r222,rmse22] = rsquare(Z22,fCurv6);

% display R^2
disp(' ')
disp(['Astig3:  R2/RMSE = ' num2str(r25,2) '/' num2str(rmse5,2) ', ' num2str(r26,2) '/' num2str(rmse6,2)])
disp(['Coma3:   R2/RMSE = ' num2str(r27,2) '/' num2str(rmse7,2) ', ' num2str(r28,2) '/' num2str(rmse8,2)])
disp(['Curv6:   R2/RMSE = ' num2str(r222,2) '/' num2str(rmse22,2)])
