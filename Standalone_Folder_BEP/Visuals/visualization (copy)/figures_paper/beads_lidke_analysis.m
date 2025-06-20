%path_output_data = '/home/idroste/Desktop/TUDelft/Data/Lidke/data_analysis/beads/With_Astigmatism/results_zstack_00*';
path_output_data = 'U:\NATaberrations\ExperimentalData\Lidke\data_analysis\DATA1\beads\Without_Astigmatism\results_zstack_00*';

all_output_files = dir(path_output_data);
nr_of_outfiles = size(all_output_files,1);
fprintf('Found %i files\n',nr_of_outfiles);

theta_combined = [];
roixy_combined = [];
mu_combined = [];
allspots_combined = [];

%for file_i=1:nr_of_outfiles
for file_i=6:9

    filename = all_output_files(file_i).name
    foldername = all_output_files(file_i).folder;
    path_output_data_full = fullfile(foldername, filename);    
    
    dataout = load(path_output_data_full,'theta','roixy','mu','allspots','outliers','params');
    
    theta_temp = dataout.theta.local;
    roixy_temp = dataout.roixy;
    mu_temp = dataout.mu;
    allspots_temp = dataout.allspots;
    outliers = dataout.outliers;
    %outliers = [];
    no_outliers = setdiff(1:size(theta_temp,2),outliers);
    params = dataout.params;

    theta_combined = cat(2,theta_combined,theta_temp(:,no_outliers));
    roixy_combined = cat(2,roixy_combined,roixy_temp(:,no_outliers));
    mu_combined = cat(4,mu_combined,mu_temp(:,:,:,no_outliers));
    allspots_combined = cat(4,allspots_combined,allspots_temp(:,:,:,no_outliers));

end


errM = get_fiterror(mu_combined,allspots_combined,params);

idx_err = find(errM(1,:)>50 | errM(2,:)>50 | errM(3,:)>5e4); % Lidke 2D

small_err = setdiff(1:size(roixy_combined,2),idx_err);

theta_combined = theta_combined(:,small_err);
roixy_combined = roixy_combined(:,small_err);
mu_combined = mu_combined(:,:,:,small_err);
allspots_combined = allspots_combined(:,:,:,small_err);

%Convert to mlambda
zernike_coeffs = theta_combined(6:end,:)*1e3/params.lambda;


%% Lidke beads: file 1-5 vs file 6-9

pixelsize = params.pixelsize;
xsize = pixelsize*params.imgSizeX;
ysize = pixelsize*params.imgSizeY;

[Xq,Yq] = meshgrid(1:xsize*1e-3);

[~,FOV_coordinates] = get_fov_coordinates(roixy_combined,theta_combined(1,:),theta_combined(2,:),params);
%X_1 = FOV_coordinates(1,:);
%Y_1 = FOV_coordinates(2,:);
%X_2 = FOV_coordinates(1,:);
%Y_2 = FOV_coordinates(2,:);
%X_3 = FOV_coordinates(1,:);
%Y_3 = FOV_coordinates(2,:);
%X_4 = FOV_coordinates(1,:);
%Y_4 = FOV_coordinates(2,:);
%X_5 = FOV_coordinates(1,:);
%Y_5 = FOV_coordinates(2,:);
%X_6 = FOV_coordinates(1,:);
%Y_6 = FOV_coordinates(2,:);
%X_7 = FOV_coordinates(1,:);
%Y_7 = FOV_coordinates(2,:);
%X_8 = FOV_coordinates(1,:);
%Y_8 = FOV_coordinates(2,:);
X_9 = FOV_coordinates(1,:);
Y_9 = FOV_coordinates(2,:);

% sample tilt
%z_1 = theta_combined(3,:);
%z_2 = theta_combined(3,:);
%z_3 = theta_combined(3,:);
%z_4 = theta_combined(3,:);
%z_5 = theta_combined(3,:);
%z_6 = theta_combined(3,:);
%z_7 = theta_combined(3,:);
%z_8 = theta_combined(3,:);
z_9 = theta_combined(3,:);

%% 3D plot per file
% Plot sample tilt
Zq_1 = griddata(X_1,Y_1,z_1,Xq,Yq,"cubic");
% Zq_2 = griddata(X_2,Y_2,z_2,Xq,Yq,"cubic");
% Zq_3 = griddata(X_3,Y_3,z_3,Xq,Yq,"cubic");
% Zq_4 = griddata(X_4,Y_4,z_4,Xq,Yq,"cubic");
% Zq_5 = griddata(X_5,Y_5,z_5,Xq,Yq,"cubic");
% Zq_6 = griddata(X_6,Y_6,z_6,Xq,Yq,"cubic");
% Zq_7 = griddata(X_7,Y_7,z_7,Xq,Yq,"cubic");
% Zq_8 = griddata(X_8,Y_8,z_8,Xq,Yq,"cubic");
% Zq_9 = griddata(X_9,Y_9,z_9,Xq,Yq,"cubic");

figure
surf(Xq,Yq,Zq_1,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',"#0072BD");
hold on
% surf(Xq,Yq,Zq_2,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',"#D95319");
% hold on
% surf(Xq,Yq,Zq_3,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',"#EDB120");
% hold on
% surf(Xq,Yq,Zq_4,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',"#4DBEEE");
% hold on
% surf(Xq,Yq,Zq_5,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',"#77AC30");
% hold on
% surf(Xq,Yq,Zq_6,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',	"#7E2F8E");
% hold on
% surf(Xq,Yq,Zq_7,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',	"#A2142F");
% hold on
% surf(Xq,Yq,Zq_8,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',	"#000000");
% hold on
% surf(Xq,Yq,Zq_9,'FaceAlpha',0.4,'EdgeColor','none','FaceColor',	"#FF00FF");
% hold on

plot3(X_1,Y_1,z_1,"o",'LineWidth',2,'Color',"#0072BD");
% hold on
% plot3(X_2,Y_2,z_2,"o",'LineWidth',2,'Color',"#D95319");
% hold on
% plot3(X_3,Y_3,z_3,"o",'LineWidth',2,'Color',"#EDB120");
% hold on
% plot3(X_4,Y_4,z_4,"o",'LineWidth',2,'Color',	"#4DBEEE");
% hold on
% plot3(X_5,Y_5,z_5,"o",'LineWidth',2,'Color',"#77AC30");
% hold on
% plot3(X_6,Y_6,z_6,"o",'LineWidth',2,'Color',"#7E2F8E");
% hold on
% plot3(X_7,Y_7,z_7,"o",'LineWidth',2,'Color',	"#A2142F");
% hold on
% plot3(X_8,Y_8,z_8,"o",'LineWidth',2,'Color',	"#000000");
% hold on
% plot3(X_9,Y_9,z_9,"o",'LineWidth',2,'Color',	"#FF00FF");
% hold on
axis square
zlim([-200 400])

xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('z (nm)')
title('Sample tilt')
%legend('file 1','file 2','file 3','file 4','file 5')
legend('file 6','file 7','file 8','file 9')
ax = gca;
ax.FontSize = 32;

%% 2D plot per file (YZ)

p1 = polyfit(Y_1,z_1,1);
p2 = polyfit(Y_2,z_2,1);
p3 = polyfit(Y_3,z_3,1);
p4 = polyfit(Y_4,z_4,1);
p5 = polyfit(Y_5,z_5,1);
p6 = polyfit(Y_6,z_6,1);
p7 = polyfit(Y_7,z_7,1);
p8 = polyfit(Y_8,z_8,1);
p9 = polyfit(Y_9,z_9,1);

l1 = p1(1)*Y_1 + p1(2);
l2 = p2(1)*Y_2 + p2(2);
l3 = p3(1)*Y_3 + p3(2);
l4 = p4(1)*Y_4 + p4(2);
l5 = p5(1)*Y_5 + p5(2);
l6 = p6(1)*Y_6 + p6(2);
l7 = p7(1)*Y_7 + p7(2);
l8 = p8(1)*Y_8 + p8(2);
l9 = p9(1)*Y_9 + p9(2);

figure
a1 = plot(Y_1,z_1,'o','MarkerSize',12,'LineWidth',4,'Color',"#0072BD");
hold on
plot(Y_1,l1,'--','Linewidth',4,'Color',"#0072BD")
hold on
a2 = plot(Y_2,z_2,'o','MarkerSize',12,'LineWidth',4,'Color',"#D95319");
hold on
plot(Y_2,l2,'--','Linewidth',4,'Color',"#D95319")
hold on
a3 = plot(Y_3,z_3,'o','MarkerSize',12,'LineWidth',4,'Color',"#EDB120");
hold on
plot(Y_3,l3,'--','Linewidth',4,'Color',"#EDB120")
hold on
a4 = plot(Y_4,z_4,'o','MarkerSize',12,'LineWidth',4,'Color',"#4DBEEE");
hold on
plot(Y_4,l4,'--','Linewidth',4,'Color',"#4DBEEE")
hold on
a5 = plot(Y_5,z_5,'o','MarkerSize',12,'LineWidth',4,'Color',"#77AC30");
hold on
plot(Y_5,l5,'--','Linewidth',4,'Color',"#77AC30")

axis square
xlabel('y(\mum)')
ylabel('z(nm)')
ylim([-175,475])
ax = gca;
ax.FontSize = 34;
legend([a1,a2,a3,a4,a5],{'loc. 1','loc. 2','loc. 3','loc. 4','loc. 5'})

figure
a6 = plot(Y_6,z_6,'o','MarkerSize',12,'LineWidth',4,'Color',"#7E2F8E");
hold on
plot(Y_6,l6,'--','Linewidth',4,'Color',"#7E2F8E")
hold on
a7 = plot(Y_7,z_7,'o','MarkerSize',12,'LineWidth',4,'Color',"#A2142F");
hold on
plot(Y_7,l7,'--','Linewidth',4,'Color',"#A2142F")
hold on
a8 = plot(Y_8,z_8,'o','MarkerSize',12,'LineWidth',4,'Color',"#000000");
hold on
plot(Y_8,l8,'--','Linewidth',4,'Color',"#000000")
hold on
a9 = plot(Y_9,z_9,'o','MarkerSize',12,'LineWidth',4,'Color',"#ff6929");
hold on
plot(Y_9,l9,'--','Linewidth',4,'Color',"#ff6929")

axis square
xlabel('y(\mum)')
ylabel('z(nm)')
ylim([-175,475])
ax = gca;
ax.FontSize = 34;
legend([a6,a7,a8,a9],{'loc. 6','loc. 7','loc. 8','loc. 9'})


%% Lidke beads: file 1-5 vs file 6-9
%Separate into two parts

pixelsize = params.pixelsize;
xsize = pixelsize*params.imgSizeX;
ysize = pixelsize*params.imgSizeY;

[Xq,Yq] = meshgrid(1:xsize*1e-3);

%zernike_coeffs_15 = theta_combined(6:end,:)*1e3/params.lambda;%% Plot aberration surfaces
zernike_coeffs_69 = theta_combined(6:end,:)*1e3/params.lambda;%% Plot aberration surfaces

[~,FOV_coordinates] = get_fov_coordinates(roixy_combined,theta_combined(1,:),theta_combined(2,:),params);
%X_15 = FOV_coordinates(1,:);
%Y_15 = FOV_coordinates(2,:);

X_69 = FOV_coordinates(1,:);
Y_69 = FOV_coordinates(2,:);

% sample tilt
%z_15 = theta_combined(3,:);
z_69 = theta_combined(3,:);

%% coma surface plot
labels = ["A(2,-2)" "A(2,2)" "A(3,-1)" "A(3,1)" "A(4,0)"];
figure
for izer = 3%1:size(params.aberrations,1)
    %hsub = subplot(2,3,izer);
    V_15 = zernike_coeffs_15(izer,:);
    Vq_15 = griddata(X_15,Y_15,V_15,Xq,Yq,"cubic");

    V_69 = zernike_coeffs_69(izer,:);
    Vq_69 = griddata(X_69,Y_69,V_69,Xq,Yq,"cubic");

    surf(Xq,Yq,Vq_15,'FaceAlpha',0.3,'EdgeColor','none','FaceColor','#00F')
    hold on
    surf(Xq,Yq,Vq_69,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','#F80')

    hold on
    plot3(X_15,Y_15,V_15,"o",'LineWidth',3,'MarkerSize',12)
    hold on
    plot3(X_69,Y_69,V_69,"o",'LineWidth',3,'MarkerSize',12)
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    zlabel('m\lambda')
    %zlim([-100 100])
    title(labels(izer))
    ax = gca;
    ax.FontSize = 16;

end
legend('Beads file 1-5','Beads file 6-9','Fontsize',32)

%% coma line plot
labels = ["A(2,-2)" "A(2,2)" "A(3,-1)" "A(3,1)" "A(4,0)"];


figure
for izer = 3%1:size(params.aberrations,1)
    %hsub = subplot(2,3,izer);
    V_15 = zernike_coeffs_15(izer,:);

    V_69 = zernike_coeffs_69(izer,:);

    p15 = polyfit(Y_15,V_15,1);
    p69 = polyfit(Y_69,V_69,1);

    l15 = p15(1)*Y_15 + p15(2);
    l69 = p69(1)*Y_69 + p69(2);

    plot(Y_15,V_15,"o",'LineWidth',4,'MarkerSize',12,'Color',"#0072BD")
    hold on
    plot(Y_69,V_69,"o",'LineWidth',4,'MarkerSize',12,'Color',"#D95319")
    hold on

    plot(Y_15(:,[1,end]),l15(:,[1,end]),'Linewidth',4,'LineStyle','--','Color',"#0072BD")
    hold on
    plot(Y_69(:,[1,end]),l69(:,[1,end]),'Linewidth',4,'LineStyle','--','Color',"#D95319")

    % surf(Xq,Yq,Vq_15,'FaceAlpha',0.3,'EdgeColor','none','FaceColor','#00F')
    % hold on
    % surf(Xq,Yq,Vq_69,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','#F80')
    % hold on

    xlabel('y (\mum)')
    ylabel('A(3,-1) (m\lambda)')
    %zlim([-100 100])
    %title(labels(izer))
    ax = gca;
    ax.FontSize = 34;

end
legend('Sample 1 (loc. 1-5)','Sample 2 (loc. 6-9)')