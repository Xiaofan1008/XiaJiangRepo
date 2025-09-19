function [rate, nTr] = psth(spikecell, bn, smoothing, maxrate, marks, flag, multicolour)
%
%  psth(spike, bn, smoothing, maxrate, marks, flag)
%
%	Spike is a cell array of spike times in ms.  This is assumed to start from 0ms.
%	bn is [start,stop] to give the start and stop times
%		of the spike times in Spike.  This realigns the spike times
%	smoothing is the degree of smoothing in the psth in ms.
%	maxrate is used for display - it specifies the maximum firing rate
%		on the y-axis
%

if nargin < 3 || isempty(smoothing); smoothing = 50; end
if nargin < 4; maxrate = 50; end
if nargin < 5; marks = []; end
if nargin < 6; flag = 0; end
if nargin < 7; multicolour = 0; end

nTr = length(spikecell);
dT = maxrate./nTr;

Start = bn(1); Stop = bn(2);

X = [];
Y = [];
for iTr = 1:nTr
    x = spikecell{iTr};
    x = (x+Start)';
    y = ones(1,length(x))*(iTr*dT);
    
    if isempty(X)
        X = x;
    else
        X = [X x];
    end
    if isempty(Y)
        Y = y;
    else
        Y = [Y y];
    end
end

if (multicolour)
    figure; hold on;
    for n = 1:length(Y)
        chk = randi(4,1);
        if nargout == 0 || flag == 1
            if X(n) > -1 && X(n) < 10
                if Y(n) < 185.7046
                    if nTr > 5
                        plot(X(n),Y(n),'b.','MarkerSize',5);
                    else
                        plot(X(n),Y(n),'b.','MarkerSize',5);
                    end
                end
                if Y(n) > 185.7046 && Y(n) < 371.4092
                    if nTr > 5
                        plot(X(n),Y(n),'r.','MarkerSize',5);
                    else
                        plot(X(n),Y(n),'r.','MarkerSize',5);
                    end
                end
                if Y(n) > 371.4092 && Y(n) < 557.1138
                    if nTr > 5
                        plot(X(n),Y(n),'g.','MarkerSize',5);
                    else
                        plot(X(n),Y(n),'g.','MarkerSize',5);
                    end
                end
                if Y(n) > 557.1138
                    if nTr > 5
                        plot(X(n),Y(n),'y.','MarkerSize',5);
                    else
                        plot(X(n),Y(n),'y.','MarkerSize',5);
                    end
                end
                continue;
            end
            if (chk == 1)
                if nTr > 5
                    plot(X(n),Y(n),'b.','MarkerSize',5);
                else
                    plot(X(n),Y(n),'b.','MarkerSize',5);
                end
            end
            if (chk == 2)
                if nTr > 5
                    plot(X(n),Y(n),'r.','MarkerSize',5);
                else
                    plot(X(n),Y(n),'r.','MarkerSize',5);
                end
            end
            if (chk == 3)
                if nTr > 5
                    plot(X(n),Y(n),'g.','MarkerSize',5);
                else
                    plot(X(n),Y(n),'g.','MarkerSize',5);
                end
            end
            if (chk == 4)
                if nTr > 5
                    plot(X(n),Y(n),'y.','MarkerSize',5);
                else
                    plot(X(n),Y(n),'y.','MarkerSize',5);
                end
            end
        end
    end
else
    if nargout == 0 || flag == 1
        if nTr > 5
            plot(X,Y,'r.','MarkerSize',5);
        else
            plot(X,Y,'r.','MarkerSize',5);
        end
    end
end

if ~isempty(marks)
    hold on;
    plot(marks,(1:nTr).*dT,'k.','Markersize',5)
end


Z = hist(X,Start:Stop);
window = normpdf((-3*smoothing:3*smoothing),0,smoothing);

rate = (1000/nTr)*conv(Z,window);
rate = rate(3*smoothing:end-3*smoothing-1);

if nargout==0 || flag == 1
hold on;
plot(Start:Stop,rate,'k','LineWidth',2);
plot([0 0],[0 maxrate],'b')
axis([Start Stop 0 maxrate])
end
