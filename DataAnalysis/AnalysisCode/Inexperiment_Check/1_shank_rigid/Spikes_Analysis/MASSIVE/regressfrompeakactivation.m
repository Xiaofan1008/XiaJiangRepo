function B=equationRegression

load('overallsigstim.mat','numsigchn','stimelectrodes')
avgnum=numsigchn./stimelectrodes;

elect=find(~isnan(avgnum));
maxelect=find(avgnum==nanmax(avgnum));
endfrompeak=avgnum(maxelect:end);
beginfrompeak=[NaN; flipud(avgnum(1:maxelect-1))];
if length(endfrompeak)>length(beginfrompeak)
    diff=length(endfrompeak)-length(beginfrompeak);
    beginfrompeak=[beginfrompeak; NaN(diff,1)];
elseif length(endfrompeak)<length(beginfrompeak)
    diff=length(beginfrompeak)-length(endfrompeak);
    endfrompeak=[endfrompeak; NaN(diff,1)];
end
distfrompeak=nansum([beginfrompeak endfrompeak],2);
equal=find((~isnan(beginfrompeak)==~isnan(endfrompeak)));
equal=intersect(find(~isnan(endfrompeak)==1),equal);
distfrompeak(equal,:)=distfrompeak(equal,:)./2;
avgnum(isnan(avgnum))=[];
x=[0 50.*(1:(length(distfrompeak)-1))];
figure
plot(elect,avgnum)
xlabel('Stimulating Electrode')
ylabel('Number of electrodes with significant spiking')

% figure
% plot(x, distfrompeak)
% title(['Distance from electrode ' num2str(maxelect) ' with peak activation'])


f = @(b,x) b(1).*exp(b(2).*x);                                     % Objective Function
B = fminsearch(@(b) (sum((distfrompeak' - f(b,x)).^2)), [-10; 0]);                % Estimate Parameters
figure
scatter(x, distfrompeak,'MarkerEdgeColor',[0 .5 .5],'MarkerFaceColor',[0 .7 .7],'LineWidth',1.5)
hold on
plot(x, f(B,x'), '-b')
hold off
xlabel('Distance (um)')
ylabel('Number of electrodes with significant spiking')
text(400, 20, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))
title(['Stimulating electrode distance from electrode ' num2str(maxelect) ' with peak activation'])

end

