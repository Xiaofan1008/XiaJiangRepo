function delLast
% Deletes the last line added to a plot
chH = get(gca,'Children');
delete(chH(1));
end