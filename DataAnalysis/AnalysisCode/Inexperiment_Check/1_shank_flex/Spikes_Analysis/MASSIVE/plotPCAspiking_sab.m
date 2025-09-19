function plotPCAspiking_sab(CHNinterest)
% used to plot PCA of first 2 and 3 principal components
% input is the channel of interest
name = pwd;
name = strsplit(name,'\');
name = name{end};
name = name(1:end-14);
load([name '.sp.mat'])

L=length(sp{CHNinterest}(:,1));
if L<1000
 figure
    for count=1:L
        hold on
        plot(sp{CHNinterest}(count,2:end))
        xlabel('Samples')
        ylabel('Voltage (uV)')
        title('Ovelayed spikes')
    end
else
    figure
    for count=1:1000
        hold on
        plot(sp{CHNinterest}(count,2:end))
        xlabel('Samples')
        ylabel('Voltage (uV)')
        title('Ovelayed spikes')
    end
end

[coeff,score,latent] = pca(sp{CHNinterest}(:,2:end));

figure
scatter(score(:,1),score(:,2),5,'filled')
xlim([-200 200])
ylim([-200 200])
xlabel('PC1')
ylabel('PC2')

figure
scatter3(score(:,1),score(:,2),score(:,3),5,'filled')
axis equal
xlabel('PC1')
ylabel('PC2')
zlabel('PC3')
end

