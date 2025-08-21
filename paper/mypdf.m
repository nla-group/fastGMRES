function mypdf(fname,r,s)

disp('mypdf disabled');
return

if nargin < 2
    r = .6; % height/width ratio
end
if nargin < 3
    s = 1; % scaling of font size
end

%return % save to PDF disabled

set(gcf,'PaperPositionMode','auto');
set(gcf,'PaperSize',s*[13,r*13]);
set(gcf,'PaperPosition',s*[0,0,13,r*13]);
print(gcf,'-dpdf', ['fig/' fname]);
%print(gcf,'-depsc2', ['fig/' fname]);
saveas(gcf,['fig/' fname '.fig'])

