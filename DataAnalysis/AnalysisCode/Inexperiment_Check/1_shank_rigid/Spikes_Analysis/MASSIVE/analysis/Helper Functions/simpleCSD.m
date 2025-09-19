%% Runs a spline CSD
% variables
gauss_sigma = 1e-4; % Gaussian Filter std. dev.
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
% electrical parameters:
cond = 0.3; % Siemens/metre
cond_top = 0.3; % Siemens/metre

% size, potential (m1 has to equal number of electrode contacts)
[m1,m2] = size(LFP);

% geometrical parameters:
diam = 5e-5; % diameter of electrodes in [m]

el_pos = (0.1:0.1:3.2)*1e-3; % Vector of positions of each electrode
el_pos(DEAD) = [];
if cond_top~=cond && (el_pos~=abs(el_pos) || length(el_pos)~=length(nonzeros(el_pos)))
    disp('Electrode contact positions must be positive when top cond. is different from ex. cond.')
    return;
end
if m1~=length(el_pos)
    disp(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
        num2str(length(el_pos)),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.']) 
    return
end

% compute spline iCSD:
Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);
[zs,CSD_cs] = make_cubic_splines(el_pos,LFP,Fcs);
%[pos1,my_CSD_spline]=new_CSD_range(zs,CSD_cs,0,2.4e-3);
if gauss_sigma~=0 %filter iCSD
  [zs,CSD_cs]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
%  [new_positions,gfiltered_spline_CSD]=gaussian_filtering(zs,CSD_cs,gauss_sigma,filter_range);
end
 
% %   plot_CSD_with_axes(new_positions,delta_t,gfiltered_spline_CSD,1)
%   [gpot_pos,gfiltered_spline_CSD_short]=new_CSD_range(new_positions,gfiltered_spline_CSD,zstart_plot,zstop_plot);

% plot CSD
unit_scale = 1e-3; % A/m^3 -> muA/mm^3
CSD_cs = CSD_cs*unit_scale;
%if nargin<3; max_plot=max(abs(CSD_matrix(:))); else; max_plot = max_plot*unit_scale; end;
%if nargin<2; scale_plot = 1; end;
max_plot=max(abs(CSD_cs(:)));

%figure();
clim=max_plot;

le = length(zs);
first_z = zs(1)-(zs(2)-zs(1))/2; %plot starts at z1-h/2;
last_z = zs(le)+(zs(le)-zs(le-1))/2; %ends at zN+h/2;
nzs = first_z:(last_z-first_z)/npoints:last_z;
zs(le+1) = zs(le)+(zs(le)-zs(le-1)); % need this in for loop
j=1; %counter
new_CSD_matrix = zeros(length(nzs),m2);
for i=1:length(nzs) % all new positions
    if nzs(i)>(zs(j)+(zs(j+1)-zs(j))/2) % > el_pos(j) + h/2
        j = min(j+1,le);
    end
    new_CSD_matrix(i,:)=CSD_cs(j,:);
end

% %increase z-resolution:
% inc_factor = 10; % 10 times the resolution
% for i = 1:length(CSD_matrix(:,1))*inc_factor
%     new_CSD_matrix(i,:) = CSD_matrix(round((i-1)/inc_factor+0.5),:);
% end;
imagesc(new_CSD_matrix,[-clim clim]);
colormap(jet);c = colorbar();
c.Label.String = 'm.\muA/mm^3';
set(gca,'FontSize',10)
[~,ntime] = size(CSD_cs);