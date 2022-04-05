function mkHDPtreeplot(PP)

nlevels = length(PP)-1;

% Plot-related params
ms = 8;
lw = 1;
c0 = [0 0 0]; c1 = [1 0 0];
% Make plots
for j = 1:nlevels
    pprev = repmat(PP{j}(:)',2,1);
    ppnext = PP{j+1};
    h = plot([pprev(:)';ppnext(:)'],[j-1;j],'-o','linewidth',lw,'markersize',ms);hold on;
    % If desired: set colors
    nlines = length(h);
    wts = linspace(0,1,nlines);
    for jj = 1:nlines
	set(h(jj),'color',c0*(1-wts(jj))+c1*(wts(jj)));
    end
end
axis ij;
ylabel('depth');
xlabel('p(1|context)');
set(gca,'xlim', [0, 1],'ytick',0:nlevels,'xtick',0:0.5:1,'tickdir','out');
box off; hold off;
