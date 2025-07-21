function SpikePlot(FS,nChn,mode)
%% Spike Plot
% FS - Sampling Rate (30kHz normally)
% nChn - number of channels (assume 64 for now, update more later)
% mode - 0:Horizontal version (rows:shanks columns:electrodes) 1: Vertical
% version (rows:electrodes, columns: shanks)
% Spike plot across all channels (Horizonal version or vertical version)


% find max and min values
files = dir('*.sp.mat'); 
if ~isempty(files)
    load(files(1).name);  % Load the first matching file
else
    error('No .sp.mat files found in the current folder.');
end
d = Depth;
nSamples = size(sp{find(~cellfun(@isempty,sp),1)},2)-1;
time_axis = (0:nSamples-1)*(1000/FS); % time axis in ms

global_max = -inf; global_min = inf;
% find average spikes
for i = 1:nChn
    if ~isempty(sp{i})
       ave_spike(i,:) = mean(sp{i}(:,2:end),1); % average spikes 
       ymin = min(ave_spike(i,:));
       ymax = max(ave_spike(i,:));
       global_max = max(global_max,ymax);
       global_min = min(global_min,ymin);        
    end
end

if mode == 0 % Horizontal mode 
    nRows = 4; nCols = 16;
    left_margin = 0.02; right_margin = 0.01; top_margin = 0.04; bottom_margin = 0.04;
    h_spacing = 0.01; v_spacing = 0.01;
    plot_width = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
    plot_height = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
    figure('Units','normalized','Position',[0 0 1 1],'Color','w'); % fullscreen
    for i = 1:nChn
        if ~isempty(sp{i})
           ave_spike(i,:) = mean(sp{i}(:,2:end),1); % average spikes 
           ymin = min(ave_spike(i,:));
           ymax = max(ave_spike(i,:));
           global_max = max(global_max,ymax);
           global_min = min(global_min,ymin);        
        end
    end
    column = 16:-1:1;
    row = [4,1,3,2]; 
    for ch = 1:nChn
        row_idx = ceil(ch/16);
        col_idx = mod(ch-1,16)+1;
    
        shank = column(col_idx); % colomn position 
        elc = row(row_idx);% row position 
    
        left = left_margin + (shank-1)*(plot_width + h_spacing);
        bottom = 1 - top_margin - elc*plot_height - (elc-1)*v_spacing;
    
        axes('Position', [left bottom plot_width plot_height]);
        
        thisChn = d(ch);
         if ~isempty(sp{thisChn})
           plot(time_axis,ave_spike(thisChn,:),'k','LineWidth',2);
           xlim([time_axis(1),time_axis(end)]);
           xlim([0 1.5]);
           ylim([global_min,global_max]);
         end
    
        if elc == nRows    
            xlabel('Time(ms)','FontSize',7);
        else 
            set(gca,'XTickLabel',[]);
        end
    
        if shank == 1
            ylabel('µV','FontSize',7);
        else 
            % set(gca,'YTickLabel',[]);
        end
        text(time_axis(2), ymin, num2str(ch), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
        set(gca,'FontSize',6,'Box','off');
    end
elseif mode == 1
        nRows = 16;
        nCols = 4;
        left_margin = 0.02;right_margin = 0.01;top_margin = 0.04;bottom_margin = 0.04;
        h_spacing = 0.01;v_spacing = 0.01;
        plot_width = (1 - left_margin - right_margin - (nCols-1)*h_spacing) / nCols;
        plot_height = (1 - top_margin - bottom_margin - (nRows-1)*v_spacing) / nRows;
        figure('Units','normalized','Position',[0 0 1 1],'Color','w'); % fullscreen
        global_max = -inf; global_min = inf;
        % find average spikes
        for i = 1:nChn
            if ~isempty(sp{i})
               ave_spike(i,:) = mean(sp{i}(:,2:end),1); % average spikes 
               ymin = min(ave_spike(i,:));
               ymax = max(ave_spike(i,:));
               global_max = max(global_max,ymax);
               global_min = min(global_min,ymin);        
            end
        end
        row = 16:-1:1;
        column = [1,4,2,3]; 
        for ch = 1:nChn    
            temp = mod(ch,16);
            if temp == 0
                loc_col = fix(ch/16);
                loc_row = 16; 
            else 
                loc_col = fix(ch/16)+1;
                loc_row = temp;
            end
            
            shank = column(loc_col); % colomn position 
            elc = row(loc_row);% row position 
        
            left = left_margin + (shank-1)*(plot_width + h_spacing);
            bottom = 1 - top_margin - elc*plot_height - (elc-1)*v_spacing;
        
            axes('Position', [left bottom plot_width plot_height]);
            
            thisChn = d(ch);
             if ~isempty(sp{thisChn})
               plot(time_axis,ave_spike(thisChn,:),'k','LineWidth',2);
               xlim([time_axis(1),time_axis(end)]);
               xlim([0 1.5]);
               ylim([global_min,global_max]);
             end
        
            % if row == nRows
            if elc == nRows    
                xlabel('Time(ms)','FontSize',7);
            else 
                set(gca,'XTickLabel',[]);
            end
        
            if shank == 1
                ylabel('µV','FontSize',7);
            else 
                % set(gca,'YTickLabel',[]);
            end
            text(time_axis(2), ymin, num2str(ch), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 12);
            set(gca,'FontSize',6,'Box','off');
        end
end