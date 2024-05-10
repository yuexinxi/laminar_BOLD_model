%{
% Steady state example
K = 6;                       % Number of depths
P = LBR_parameters(K);       % Get parameter structure with default values
cbf      = kron([1.6,1.6,1.3,1.3,1.6,1.6],ones(P.T/P.dt,1)); % Define model input (Relative blood flow across depths)
[LBR,Y] = LBR_model(P,cbf);  % Generate the LBR
figure(1), 
plot(P.l,flipud(LBR(end,:)')); % plot the LBR profile as a function of normalized cortical depth
xlim([0 100]); ylim([0 6]); xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
%}


% Dynamic response example
K = 20;  
P.N = neuronal_NVC_parameters(K);  % consider default parameters
P.N.T = 60; % (in seconds)
P.H = LBR_parameters(K);
P.H.T  = P.N.T;
P.H.dt = P.N.dt;
P.H.alpha_v   = 0.35;
P.H.alpha_d   = 0.2;
P.H.tau_d_de  = 30;
U.u = zeros(P.N.T/P.N.dt,K);
dur = 30/P.N.dt; % 2 sec stimulus
onset     = 3/P.N.dt; % when K<0, P.N.dt = 0.01
offset    = onset + dur;
% input stimuli
for i = 1:10
    U.u(onset:offset,i) = 0.1*i;
    if i > 5
        U.u(onset:offset,i) = U.u(onset:offset,i-1)-0.05;
    end
    U.u(onset:offset,21-i) = U.u(onset:offset,i);
end
% get neuronal respone and cbf
P.s_d = 0.4;
P.alpha_v = 0.25;
p.alpha_d = 0.1;
[neuro, cbf]  = neuronal_NVC_model(P.N,U);
% get LBR
[LBR,Y]       = LBR_model(P.H,cbf);
time_axis = [0:P.H.dt:P.H.T-P.H.dt];
%figure(1), 
%subplot(141), plot(time_axis,U.u(:,4)); xlim([time_axis(1), time_axis(end)]); ylim([-1 1]);
%xlabel('Time (s)'); ylabel('Relative neuronal response (%)'); axis square;
%subplot(142), plot(time_axis,neuro); xlim([time_axis(1), time_axis(end)]); ylim([-1 1]);
%xlabel('Time (s)'); ylabel('Relative neuronal response (%)'); axis square;
%subplot(143), plot(time_axis,cbf); xlim([time_axis(1), time_axis(end)]); ylim([0.5 2]); 
%xlabel('Time (s)'); ylabel('Relative CBF (%)'); axis square;
%subplot(144), plot(time_axis,LBR); xlim([time_axis(1), time_axis(end)]); ylim([-1 4]);  
%xlabel('Time (s)'); ylabel('LBR (%)'); axis square;
%legend('1','2','3','4','5','6')