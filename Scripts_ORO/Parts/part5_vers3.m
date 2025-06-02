% this function is the second version of a possible solution
% note that in this case ToF is not fixed

function []=part5_vers3(rho,CD,A,m,r2,n2,rH_err,rH_dot_err,deltaV1,deltaVtot_vers2)


%%

coeff_ad=1/2*rho*CD*A/m *r2;
    N=1000;
    % attention to both the span and the number of points N
    ToF=linspace(0,2*pi,N);
    ToF=ToF(2:end); % exclude ToF=0
    bet=deltaV1; % first guess
    
    % hide fsolve output in command window
    options = optimoptions("fsolve", "Display", "off");
    warning("off", "all");
    
    T=zeros(1,length(ToF)); % initialization
    tot_impulse=zeros(1,length(ToF));
    for i=1:length(ToF)
        t=ToF(i); % considered ToF for the iteration
        tspan=[0 t];

        % find 1st impulse such that zero of the position at the end of the rendezvous
        [res,~,~,~]=fsolve(@(deltaV_opt) pos_err(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+deltaV_opt],tspan),bet,options);
        [~,csi,eta,~,v_csi,v_eta,~]=solveHCW_drag(coeff_ad,[rH_err/r2;rH_dot_err/r2/n2+res],tspan);
        deltaVtot=norm(res,2)+norm([-v_csi(end);-v_eta(end);0],2);

        % max position offset deeming the rendezvous solution acceptable
        if sqrt(csi(end)^2+eta(end)^2)*r2*1000<10^(-5)
            T(i)=t;
            tot_impulse(i)=deltaVtot;
        end
    end
    warning("on", "all");
    
    idx_0=tot_impulse==0;
    tot_impulse(idx_0)=[];
    T(idx_0)=[];
    [impulse_min, idx_min]=min(tot_impulse);
    T_min=T(idx_min);
    figure;
    hold on
    plot(T,tot_impulse,"b-o")
    plot(T_min,impulse_min,"r*")
    title("Values of total impulse for different ToF")
    %ylim([0 0.2])
    xlabel("ToF [-]")
    ylabel("\Delta V [-]")
    hold off

    percentage=impulse_min/deltaVtot_vers2*100;

fprintf('Rendezvous in ToF: %.4f s\n', T_min/n2); 
fprintf('using total impulse: %.4f m/s\n', impulse_min*r2*n2*1000);   
fprintf('This impulse is %4f%% of the one calculated with fixed ToF. \n', percentage);

end