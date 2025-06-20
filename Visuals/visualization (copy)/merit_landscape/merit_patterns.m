
meritstore_last = meritstore_local(:,end);
merit_outliers = meritstore_last(outliers);
no_outliers = setdiff(1:params.Ncfg,outliers);
merit_no_outliers = meritstore_last(no_outliers);

mean(merit_outliers)
mean(merit_no_outliers)

figure
bar(merit_outliers)

figure
bar(merit_no_outliers)

%%
figure
scatter(roixy(1,1:50),roixy(2,1:50))

figure
scatter(roixy(1,:),roixy(2,:))