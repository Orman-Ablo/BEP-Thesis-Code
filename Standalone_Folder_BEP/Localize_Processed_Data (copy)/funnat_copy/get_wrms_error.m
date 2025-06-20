function [wrms_error_avg,wrms_error] = get_wrms_error(theta,groundtruth,roixy,params,make_plot)
%get_wrms_error Computes the average Wrms error over the spots. It also
%makes a plot of the distribution of the errors.

pixelsize = params.pixelsize;
lambda = params.lambda;

groundtruth_local = groundtruth.local;
groundtruth_global = groundtruth.global;
theta_local = theta.local;
theta_global = theta.global;

% Calculate the groundtruth zernike coefficients.
%gt_field_coordinates_x = (roixy(1,:)*pixelsize + groundtruth_local(1,:))*1e-03;
%gt_field_coordinates_y = (roixy(2,:)*pixelsize + groundtruth_local(2,:))*1e-03;

fov_coordinates = get_fov_coordinates(roixy,groundtruth_local(1,:),groundtruth_local(1,:),params);
zernikeCoefficients = get_zernike_coefficients(fov_coordinates(1,:), fov_coordinates(2,:), groundtruth_global, params);
gt_zernikeCoefficients = (1e3/lambda)*zernikeCoefficients;

% Calculate the fitted zernike coefficients.
%fitted_field_coordinates_x = (roixy(1,:)*pixelsize + theta_local(1,:))*1e-03;
%fitted_field_coordinates_y = (roixy(2,:)*pixelsize + theta_local(2,:))*1e-03;

fov_coordinates = get_fov_coordinates(roixy,theta_local(1,:),theta_local(1,:),params);
zernikeCoefficients = get_zernike_coefficients(fov_coordinates(1,:), fov_coordinates(2,:), theta_global, params);

fitted_zernikeCoefficients = (1e3/lambda)*zernikeCoefficients;

% Calculate the average Wrms value.
wrms_error = sqrt(sum((gt_zernikeCoefficients-fitted_zernikeCoefficients).^2,2));
mean(wrms_error,1);
wrms_error_avg = mean(wrms_error,1);
wrms_error_std = std(wrms_error(:));

% Plot Wrms error distribution
if make_plot
    y_axis = wrms_error;
    x_axis = 0.4*std(y_axis)*randn(numel(y_axis),1);
    
    axes(figure, 'NextPlot', 'add', 'XColor', 'none');
    points = scatter(x_axis,y_axis,75,'Marker','.','MarkerEdgeColor','black');
    hold on
    stdev = errorbar(0,wrms_error_avg,wrms_error_std,'LineWidth',3,'MarkerSize',10,'Marker','x','MarkerEdgeColor','blue','Color','blue','CapSize',50);
    
    axis equal
    legend([points stdev],{'error', '\sigma(error)'});
    ylabel('Wrms (m\lambda)','FontSize',20);
    title('Wrms error');
    ax=gca;
    ax.FontSize = 20;
end

end