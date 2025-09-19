function X_line=best_fit_3D_line(X)

%https://au.mathworks.com/matlabcentral/answers/406619-3d-coordinates-line-of-fit
% Generate sample dataset
% -------------------------------------------------------------------------
% % Standard deviation of noise
% s=1;
% % Random direction in space
% r=randn(3,1);
% r=r/(norm(r)+eps);
% % Random data points along r
% N=1E3; % number of point samples
% t=(10*s)*(2*rand(N,1)-1);
% Xo=bsxfun(@times,t,r');
% % Add (isotropic) Gaussian noise to Xo
% X=Xo+s*randn(N,3);
% % Offset X (relative to the origin) by a random amount
% X=bsxfun(@plus,X,5*s*randn(1,3));
% % Find line of best fit (in least-squares sense) through X
% -------------------------------------------------------------------------
N=length(X(:,1));
X_ave=mean(X,1);            % mean; line of best fit will pass through this point  
dX=bsxfun(@minus,X,X_ave);  % residuals
C=(dX'*dX)/(N-1);           % variance-covariance matrix of X
[R,D]=svd(C,0);             % singular value decomposition of C; C=R*D*R'
% NOTES:
% 1) Direction of best fit line corresponds to R(:,1)
% 2) R(:,1) is the direction of maximum variances of dX 
% 3) D(1,1) is the variance of dX after projection on R(:,1)
% 4) Parametric equation of best fit line: L(t)=X_ave+t*R(:,1)', where t is a real number
% 5) Total variance of X = trace(D)
% Coefficient of determineation; R^2 = (explained variance)/(total variance)
D=diag(D);
R2=D(1)/sum(D);
% Visualize X and line of best fit
% -------------------------------------------------------------------------
% End-points of a best-fit line (segment); used for visualization only 
x=dX*R(:,1);    % project residuals on R(:,1) 
x_min=min(x);
x_max=max(x);
dx=x_max-x_min;
Xa=(x_min-0.05*dx)*R(:,1)' + X_ave;
Xb=(x_max+0.05*dx)*R(:,1)' + X_ave;
X_line=[Xa;Xb];
% figure('color','w')
% axis equal 
% hold on
% plot3(X_line(:,1),X_line(:,2),X_line(:,3),'-r','LineWidth',3) % best fit line 
% plot3(X(:,1),X(:,2),X(:,3),'.k','MarkerSize',13)           % simulated noisy data
% set(get(gca,'Title'),'String',sprintf('R^2 = %.3f',R2),'FontSize',25,'FontWeight','normal')
% xlabel('X','FontSize',20,'Color','k')
% ylabel('Y','FontSize',20,'Color','k')
% zlabel('Z','FontSize',20,'Color','k')
% view([20 20])
% drawnow
end