function [Coeffsout, gof2] = fit_3rd_degree_poly_no_constant2(xdata,ydata,x0)
xdata = xdata/x0;

%converting to column vectors
if isrow(xdata)
    xdata = xdata';
end
if isrow(ydata)
    ydata = ydata';
end
a0 = abs((max(ydata)-min(ydata)) / (max(xdata)-min(xdata)).^3);
b0 = abs((max(ydata)-min(ydata)) / (max(xdata)-min(xdata)).^2);

a0 = max([2e-8,a0]);
b0 = max([2e-8,b0]);

%set weights proportional to values
weights = (ydata + (ydata == 0).*min(ydata(ydata ~= 0))).^2;  %assigning weight of 1 to points with value zero
s = fitoptions('Method','NonlinearLeastSquares',...
               'Display','off',...
               'Lower',[-10*a0,-10*b0],...
               'Upper',[10*a0,10*b0],...
               'Startpoint',[a0,b0],...
               'Weights',weights);
f = fittype('a*x.^3 + b*x.^2+ (1-a-b)*x','options',s);
[c2,gof2] = fit(xdata,ydata,f);
Coeffsout = coeffvalues(c2);

Coeffsout = [Coeffsout,(1-Coeffsout(1) - Coeffsout(2))/x0];
Coeffsout(1) = Coeffsout(1)/x0^3;
Coeffsout(2) = Coeffsout(2)/x0^2;
end