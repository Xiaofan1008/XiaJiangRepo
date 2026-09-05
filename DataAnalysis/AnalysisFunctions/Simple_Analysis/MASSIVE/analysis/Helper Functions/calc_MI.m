function MI = calc_MI(P,N,T)
%% Takes a phase distribution and calculates the Tort Modulation Index
MI = zeros(size(P(1),1),1);
if (T == 2)
    for i = 1:size(P,1)
        %P(i,:) = P(i,:) ./ mean(P(i,:));
        MI(i) = (log(N) - (-nansum(P(i,:).*log(P(i,:))))) ./ log(N);
    end
end
%% Alternate method for calculation
if (T == 1)
    for i = 1:size(P,1)
        MI(i) = diff([min(P(i,:)),max(P(i,:))]);
    end
end
end