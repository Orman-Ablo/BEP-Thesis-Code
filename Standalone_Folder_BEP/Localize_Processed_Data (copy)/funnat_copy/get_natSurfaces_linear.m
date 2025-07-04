function Zsurface = get_natSurfaces(Xim,Yim,RAstig3,RComa3,RTrefoil,RCurv5,RAstig5,RComa5,RCurv6)

% Z(5,6) first order astigmatism
Zsurface(:,:,1) = RAstig3(4)*(Xim.^2+Yim.^2) +...
    RAstig3(8)*Yim - RAstig3(9).*Xim + RAstig3(10)*(2*Xim.*Yim) +...
    RAstig3(11)*Xim+RAstig3(12)*Yim+RAstig3(13);
Zsurface(:,:,2) = RAstig3(5)*(Xim.^2+Yim.^2)+...
    RAstig3(8).*Xim+...
    RAstig3(9)*Yim+RAstig3(10)*(Yim.^2-Xim.^2)+RAstig3(11)*Yim-...
    RAstig3(12)*Xim+RAstig3(14);

% Z(7,8) first order coma
Zsurface(:,:,3) = RComa3(4)*Xim+RComa3(7)*Yim+RComa3(8)*Xim+RComa3(9);
Zsurface(:,:,4) = RComa3(4)*Yim+RComa3(7)*Xim-RComa3(8)*Yim+RComa3(10);

% Z(9,10) trefoil
Zsurface(:,:,5) = RTrefoil(1)*(3*Yim.^2.*Xim-Xim.^3)+RTrefoil(2)*(Yim.^2-Xim.^2)+...
    RTrefoil(3)*(2*Xim.*Yim)+RTrefoil(4)*Xim+RTrefoil(5)*Yim+RTrefoil(6);
Zsurface(:,:,6) = RTrefoil(1)*(Yim.^3-3*Xim.^2.*Yim)-RTrefoil(2)*(2*Xim.*Yim)+...
    RTrefoil(3)*(Yim.^2-Xim.^2)+RTrefoil(4)*Yim-RTrefoil(5)*Xim+RTrefoil(7);

% Z(11) first order spherical
Zsurface(:,:,7) = RCurv5(1)*(Xim.^2+Yim.^2)+RCurv5(2)*Xim+RCurv5(3)*Yim+RCurv5(4);

% Z(12,13)
Zsurface(:,:,8) = RAstig5(1)*(2*Xim.*Yim)+RAstig5(2)*Yim+RAstig5(3)*Xim+RAstig5(4);
Zsurface(:,:,9) = RAstig5(1)*(Yim.^2-Xim.^2)-RAstig5(2)*Xim+RAstig5(3)*Yim+RAstig5(5);

% Z(16,17)
Zsurface(:,:,10) = RComa5(1)*Xim+RComa5(2);
Zsurface(:,:,11) = RComa5(1)*Yim+RComa5(3);

% Z(6,0)
Zsurface(:,:,12) = RCurv6(1)*(Xim.^2+Yim.^2)+RCurv6(2)*Xim+RCurv6(3)*Yim+RCurv6(4);
