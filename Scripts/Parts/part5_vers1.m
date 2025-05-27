% this function is the first version of a possible solution

function []=part5_vers1(rho,CD,A,m,Delta_T,r2,n2,rH_err,rH_dot_err,deltaV1)
%% 

addpath("Functions\General_Functions")
addpath("Functions\Function_Part2")
addpath("Functions\Function_Part3")
addpath("Functions\Function_Part5")


%% solution, monotonic decreasing err

% small perturbation-->assuming there is a solution in a certain range of the ideal solution
% ToF is fixed as a design parameter already optimezed 
% with fixed ToF the solution of the rendezvous is described by 2 parameters csi and eta (for plane motion)
% for a fixed starting point the different trajectories are governed by 2 parameters
% the solution is supposed unique since there are 2 free parameters (chosen impulse) and 2 constraints (final position)
% to avoid loosing the solution from the range, it is implemented a control 
% in order to have a monotonically decreasing error each accepted iteration
% A better estimate on the solution range would be beneficial, in this case
% N is set equal to 2

coeff_ad=1/2*rho*CD*A/m *r2;
tspan=[0 Delta_T]; % fixed ToF
N=2;
exit=false;

Nsteps=11;
iter=1;
itermax=50;
err=zeros(itermax,1);

figure;
hold on

err_pre=inf; % just for first cycle
while exit==false & iter<itermax
    range=deltaV1/N;
    % different values of impulse in the estimated range
    steps1=linspace(deltaV1(1)-abs(range(1)),deltaV1(1)+abs(range(1)),Nsteps);
    steps2=linspace(deltaV1(2)-abs(range(2)),deltaV1(2)+abs(range(2)),Nsteps);

    minima2=zeros(length(steps1),1);
    delta_eta_studied=zeros(length(steps1),1);
    temporary_mem=zeros(length(steps2),1);

    for i=1:length(steps1)
        for j=1:length(steps2)
            deltaV1_opt=[steps1(i);steps2(j);0]; % trying different impulse every cycle
            [~,csi,eta,~,~,~,~] = solveHCW_drag(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV1_opt],tspan);
            min_pos=sqrt(csi(end)^2+eta(end)^2); % min pos offset for impulse ij
            if min_pos*r2*1000<10^(-6) % rendezvous at ToF threshold
                % fprintf('Final position offset: %.4f [m]\n', min_pos*r2*1000);
                % fprintf('with new first impulse 10^(-3): %.4f [m/s]\n', [steps1(i);steps2(j);0]*r2*n2*1000);
                exit = true;
            elseif sum((abs(csi)+abs(eta))*r2*1000<10^(-7))>0
                fprintf('At a certain time<ToF there is a rendezvous: %.4f [-]\n', true);
            end
            temporary_mem(j)=min_pos;
        end
        [min2,idx]=min(temporary_mem); % min pos offset for impulse_csi i
        minima2(i)=min2;
        delta_eta_studied(i)=steps2(idx); % associated impulse_eta for min pos offset
    end

    [abs_minima,idx_abs]=min(minima2); % min pos offset trying all analyzed impulses
    err(iter)=abs_minima*r2*1000;
    if iter>1 % this condition is used to avoid error in the first iteration
        err_pre=err(iter-1);
    end
    if err(iter)>err_pre % check if err is decreasing
        Nsteps=25; % use more values in the range to get monotonic decrease
        plot(iter,err(iter),"r*") % plotted in a different color since the found value is not used in the iteration
    else
        deltaV1(1)=steps1(idx_abs);
        deltaV1(2)=delta_eta_studied(idx_abs);
        N=N*2; % reducing the investigated range
        Nsteps=10; % default value
        plot(iter,err(iter),"b*") % point actually used in the iteration
        iter=iter+1;
    end
    
    %iter=iter+1;
end
hold off

if iter==itermax
    warning("Maximum chosen number of iteration reached")
end

% integrate with new impulse for a check and to find second impulse needed
[~,csi,eta,~,v_csi,v_eta,v_zeta] = solveHCW_drag(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV1],tspan);
deltaV2=-[v_csi(end);v_eta(end);v_zeta(end)]; % second impulse
deltaVtot=norm(deltaV1,2)+norm(deltaV2,2); % total impulse

disp(" ==================== Part 5 - Optimisation of Delta V and ToF ==================== ")
disp(" ")

fprintf('Final position offset with drag, relative CCS: %.4f m \n', norm([csi(end);eta(end)],2)*r2*1000);
fprintf('Time of Flight (ToF): %.4f s \n', Delta_T/n2);
fprintf('Total impulse: %.4f m/s \n', deltaVtot*r2*n2*1000);
disp(" ")
end