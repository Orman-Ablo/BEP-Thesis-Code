
input1 ="C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\CRLB_Loc_Combined_initial.mat";
input2 = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\CRLB_Loc_Combined.mat";
input2 = "C:\Users\Naam\Documents\BEP\BEP Code Roman\Combining Results\CRLB_Loc_Combined_Second.mat";

initial = load(input1);

locations1 = initial.locations;
CRLB1 = initial.CRLB;

second = load(input2);

CRLB2 = second.CRLB;
locations2 =second.locations;

Npatch = 53;

Npatch_z = 53;

minxyz = min(locations1);
maxxyz = max(locations1);


Xpatch = minxyz(1)+(0:Npatch)*(maxxyz(1)-minxyz(1))/Npatch;
Ypatch = minxyz(2)+(0:Npatch)*(maxxyz(2)-minxyz(2))/Npatch;
Zpatch = minxyz(3)+(0:Npatch_z)*(maxxyz(3)-minxyz(3))/Npatch_z;
patches = [Xpatch; Ypatch; Zpatch];

patch_mid = (patches(:,1:end-1)+ patches(:,2:end))/2;

Zpatch = linspace(-503.5, 503.5,Npatch);

zz =(Zpatch(:,1:end-1)+ Zpatch(:,2:end))/2;





numdims= 3;
dimstore = zeros(numdims,Npatch-1,2);

% for ii = 1:numdims
%     for jj = 1:Npatch
%         filter1 = (locations1(:,ii)>patches(ii,jj))&(locations1(:,ii)<patches(ii,jj+1));
%         filter2 = (locations2(:,ii)>patches(ii,jj))&(locations2(:,ii)<patches(ii,jj+1));
%         CRLB_val1 = mean(CRLB1(filter1,ii));
%         CRLB_val2 = mean(CRLB2(filter2,ii));
% 
%         dimstore(ii,jj,1) = CRLB_val1;
%         dimstore(ii,jj,2) = CRLB_val2;
% 
% 
%     end
% end

Npatch = 52;

for ii = 1:numdims
    for jj = 1:Npatch
        filter1 = (locations1(:,3)>Zpatch(jj))&(locations1(:,3)<Zpatch(jj+1));
        filter2 = (locations2(:,3)>Zpatch(jj))&(locations2(:,3)<Zpatch(jj+1));
        CRLB_val1 = mean(CRLB1(filter1,ii));
        CRLB_val2 = mean(CRLB2(filter2,ii));

        dimstore(ii,jj,1) = CRLB_val1;
        dimstore(ii,jj,2) = CRLB_val2;


    end
end


str = ['X' 'Y' 'Z'];
figure;


hold on;



plot(zz, dimstore(1,:,1),'LineWidth',1.5,'DisplayName',[str(1),' intitial'],'Color','#ff0000')
plot(zz, dimstore(1,:,2),'LineWidth',1.5,'DisplayName',[str(1),' second'],'Color','#ff9e00')
plot(zz, dimstore(2,:,1),'LineWidth',1.5,'DisplayName',[str(2),' intitial'],'LineStyle','--','Color','#0013ff')
plot(zz, dimstore(2,:,2),'LineWidth',1.5,'DisplayName',[str(2),' second'],'LineStyle','--','Color','#18dafc')
plot(zz, dimstore(3,:,1),'LineWidth',1.5,'DisplayName',[str(3),' intitial'],'Color','#1bff00')
plot(zz, dimstore(3,:,2),'LineWidth',1.5,'DisplayName',[str(3),' second'],'Color','#048500')



legend
xlim([-503.5 503.5])  
xlabel('Z-Position (nm)')
ylabel('CRLB (nm)')
title('CRLB over z-positions')
fontsize(16,'points')
hold off;


% xx2 = (xx(1:end-1)+xx(2:end))/2;
% 
% count = zeros(numdims, intervals);
% 
% 
% for ii = 1:numdims
%     for jj = 1:intervals
%     positionmask =  (locations(:,ii)>xx(jj)) & (locations(:,ii)< xx(jj+1));
%     % disp(sum(positionmask))
%     dimstore(ii,jj) = mean(CRLB(:,ii));
%     end
%     figure;
%     plot(xx2, dimstore(ii,:))
%     xlabel('range (nm)')
%     ylabel('CRLB (nm)')
% 
% end
