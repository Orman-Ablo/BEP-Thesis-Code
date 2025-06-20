function outliers = get_outliers(theta,meritstore,numiters,mu,allspots,params)
% outlier removal based on statistics of maximum log-likelihood, on spot
% size, and on found position

Nitermax = params.max_local_iterations;

if contains(params.fitmodel,'xyz')
    box_size = 5;
else
    box_size = 3;
end

outliers_type1 = (abs(theta(1,:))>box_size*params.pixelsize)|(abs(theta(2,:))>box_size*params.pixelsize);
outliers_type3 = (numiters==Nitermax+1);
%errM = get_fiterror(mu,allspots,params);
%outliers_type4 = []; %errM(1,:)>1.5;

meritfinal = squeeze(meritstore(:,end))';
no_outliers_type3 = setdiff(1:numel(meritfinal),outliers_type3);
meanmeritfinal = mean(meritfinal(no_outliers_type3));
stdmeritfinal = std(meritfinal(no_outliers_type3));
outliers_type2 = (meritfinal<meanmeritfinal-3*stdmeritfinal); % or 2x std is better?
outliers = outliers_type1 | outliers_type2 | outliers_type3;

fprintf(['\nNumber of spots outside ROI center: ' num2str(sum(outliers_type1)) '\n']);
fprintf(['Number of spots with merit too small: ' num2str(sum(outliers_type2)) '\n']);
fprintf(['Number of unconverged spots: ' num2str(sum(outliers_type3)) '\n']);
%fprintf(['Number spots with large fiterror: ' num2str(sum(outliers_type4)) '\n']);
fprintf(['\nNumber of outliers: ' num2str(sum(outliers)) '\n']);
