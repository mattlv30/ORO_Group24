% this function is aims show graphs that may be of interest to give a
% context to the developed part5 function, in particular to the
% first version solution "solution, monotonic decreasing err"

function [] = pre_part5(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)

%%
addpath("Functions\General_Functions")
addpath("Functions\Function_Part2")
addpath("Functions\Function_Part3")
addpath("Functions\Function_Part5")

%%

coeff_ad=1/2*rho*CD*A/m *r2;
tspan=[0 Delta_T]; % fixed ToF

range=norm(deltaV1,2)/200;
steps1=linspace(deltaV1(1)-abs(range),deltaV1(1)+abs(range),8);
steps2=linspace(deltaV1(2)-abs(range),deltaV1(2)+abs(range),50);

c=["k*","b*","c*","r*","ko","bo","co","ro"]; %4
 figure;
 hold on
minima2=zeros(length(steps1),1);
delta_eta_studied=zeros(length(steps1),1);
temporary_mem=zeros(length(steps2),1);
A=zeros(length(steps1),length(steps2)); %
for i=1:length(steps1)
    for j=1:length(steps2)
        deltaV1_opt=[steps1(i);steps2(j);0];
        [~,csi,eta,~,~,~,~] = solveHCW_drag(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV1_opt],tspan);
        plot(deltaV1_opt(2)*r2*n2,sqrt(csi(end)^2+eta(end)^2)*r2*1000,c(i))
        temporary_mem(j)=sqrt(csi(end)^2+eta(end)^2);
        A(i,j)=sqrt(csi(end)^2+eta(end)^2); %
    end
    [min2,idx]=min(temporary_mem);
    minima2(i)=min2;
    delta_eta_studied(i)=steps2(idx);
end
xlabel("\Delta V_{\eta} [m/s]")
ylabel("err [m]")
title("Final position offset for different first impulse")
hold off

% select minima from previous plot and graph them
figure;
hold on
plot(delta_eta_studied*r2*n2*1000,minima2*r2*1000,"-o")
xline(deltaV1(2)*r2*n2*1000)
xlabel("\Delta V_{\eta} [m/s]")
ylabel("err [m]")
title("Minimum offset for each associated impulse")
hold off

% minima plot associated to their impulse

[X,Y] = meshgrid(steps2,steps1); % create Matrixes for impulse components
figure;
extracted=sort(A(:))*r2*1000;
% the sorted values are sampled for readibility
levels=[extracted(1:2);extracted(10:10:50) ;extracted(55:20:95); extracted(100:50:end)];
contourf(X*r2*n2*1000,Y*r2*n2*1000,A*r2*1000,levels,"ShowText",true)
title('Position offset in the rendezvous');
xlabel('\DeltaV_{\eta} [m/s]');
ylabel('\DeltaV_{\xi} [m/s]');

figure;
surf(X, Y, A);  % Create 3D surface
colorbar;       % show colorbar
colormap jet;   % colormap style

title('Position offset in the rendezvous for different first impulse components');
xlabel('\DeltaV_{\eta} [-]');
ylabel('\DeltaV_{\xi} [-]');
zlabel('err [-]');
end