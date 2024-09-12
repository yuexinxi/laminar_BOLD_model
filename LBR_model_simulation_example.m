% Modified by Yuexin 2024-05
% To replicate Havlicek et Uludag, 2020. Figure 8
clear all;

% Steady state example
%{
K = 6;                       % Number of depths
P = LBR_parameters(K);       % Get parameter structure with default values
cbf      = kron([1.6,1.6,1.3,1.3,1.6,1.6],ones(P.T/P.dt,1)); % Define model input (Relative blood flow across depths)
[LBR,Y] = LBR_model(P,cbf);  % Generate the LBR
figure(1), 
plot(P.l,flipud(LBR(end,:)')); % plot the LBR profile as a function of normalized cortical depth
xlim([0 100]); ylim([0 6]); xlabel('1 - Cortical depth (%)'); ylabel('LBR (%)'); axis square;
%}


% Dynamic response example

K = 21;  
P.N = neuronal_NVC_parameters(K);  % consider default parameters
P.N.T = 60; % (in seconds)
P.H = LBR_parameters(K);
P.H.T  = P.N.T;
P.H.dt = P.N.dt;
P.N.mu = 1.5;
U.u = zeros(P.N.T/P.N.dt,K);
dur = 30/P.N.dt; % 30 sec stimulus
onset     = 3/P.N.dt; % when K<0, P.N.dt = 0.01
offset    = onset + dur;
% input stimuli
for i = 1:11
    U.u(onset:offset,i) = 0.24 + 0.06*i;
    if i > 6
        U.u(onset:offset,i) = U.u(onset:offset,i-1)-0.03;
    end
   U.u(onset:offset,22-i) = U.u(onset:offset,i);
end
%U.u(onset:offset,:) = 1;
% get neuronal respone and cbf
P.H.s_d = 0.4;
P.H.alpha_v = 0.25;
P.H.alpha_d = 0.1;
[neuro, cbf]  = neuronal_NVC_model(P.N,U);
% get LBR
P.H.tau_v_in = 10;
P.H.tau_v_de = 10;
P.H.tau_d_in = 40;
P.H.tau_d_de = 100;
[LBR,Y]       = LBR_model(P.H,cbf);
time_axis = [0:P.H.dt:P.H.T-P.H.dt];
cbf_response = (cbf-1).*100;
LBR_response = LBR;

time_axis = [0:P.H.dt:P.H.T-P.H.dt] - onset*P.N.dt; % time axis in seconds

[Peak_Amp_cbf,Peak_Pos_cbf] = max(cbf_response(onset:end,:));
[PSU_Amp_cbf,PSU_Pos_cbf]   = min(cbf_response(offset:end,:));
TTP_cbf = time_axis(onset+Peak_Pos_cbf)';
TTU_cbf = time_axis(offset+PSU_Pos_cbf)'; 

figure(1)
imagesc(time_axis,P.H.l,cbf_response');
hold on; xlim([-3 57]); ylim([0 100]);
xlabel('Time (s)'); ylabel('1 - Cortical depth (%)'); axis square; title('Laminar CBF response');
cb = colorbar;
cb.Title.String = "CBF (%)";
plot(TTP_cbf, P.H.l,'.');
hold on;
plot(TTU_cbf, P.H.l,'.');

[Peak_Amp_LBR, Peak_Pos_LBR] = max(LBR_response(onset:end,:));
[PSU_Amp_LBR,PSU_Pos_LBR]   = min(LBR_response(offset:end,:));
TTP_LBR = time_axis(onset+Peak_Pos_LBR)';
TTU_LBR = time_axis(offset+PSU_Pos_LBR)'; 

figure(2)
imagesc(time_axis,P.H.l,LBR_response');
hold on; xlim([-3 57]); ylim([0 100]);
xlabel('Time (s)'); ylabel('1 - Cortical depth (%)'); axis square; title('Laminar BOLD response');
cb = colorbar;
cb.Title.String = "BOLD (%)";
plot(TTP_LBR, P.H.l,'.');
hold on;
plot(TTU_LBR, P.H.l,'.');

pos_res_LBR = mean(LBR_response(onset+10/P.N.dt:onset+30/P.N.dt,:),1);
PSU_res_LBR = mean(LBR_response(onset+35/P.N.dt:onset+42/P.N.dt,:),1);
pos_res_cbf = mean(cbf_response(onset+10/P.N.dt:onset+30/P.N.dt,:),1);
PSU_res_cbf = mean(cbf_response(onset+35/P.N.dt:onset+42/P.N.dt,:),1);

%figure(1), 
%subplot(141), plot(time_axis,U.u(:,4)); xlim([time_axis(1), time_axis(end)]); ylim([-1 1]);
%xlabel('Time (s)'); ylabel('Exogenous input(%)'); axis square;
%subplot(142), plot(time_axis,neuro); xlim([time_axis(1), time_axis(end)]); ylim([-1 1]);
%xlabel('Time (s)'); ylabel('Relative neuronal response (%)'); axis square;
%subplot(143), plot(time_axis,cbf(:,5)); xlim([time_axis(1), time_axis(end)]); ylim([0.5 3]); 
%xlabel('Time (s)'); ylabel('Relative CBF (%)'); axis square;
%subplot(144), plot(time_axis,LBR); xlim([time_axis(1), time_axis(end)]); ylim([-2 8]);  
%xlabel('Time (s)'); ylabel('LBR (%)'); axis square;
%legend('1','2','3','4','5','6')
figure(3)
subplot(141),
plot(P.H.l,flipud(pos_res_cbf'),'-'); hold on; xlim([0 100]); ylim([20 100]); 
xlabel('1 - Cortical depth (%)'); ylabel('CBF (%)'); axis square; title('Positive Response')
subplot(142),
plot(P.H.l,flipud(pos_res_LBR'),'-'); hold on; xlim([0 100]); ylim([0 9]); 
xlabel('1 - Cortical depth (%)'); ylabel('BOLD signal (%)'); axis square; title('Positive Response')
subplot(143),
plot(P.H.l,flipud(PSU_res_cbf'),'-'); hold on; xlim([0 100]); ylim([-12 0]);
xlabel('1 - Cortical depth (%)'); ylabel('CBF (%)'); axis square; title('PSU')
subplot(144),
plot(P.H.l,flipud(PSU_res_LBR'),'-'); hold on; xlim([0 100]); ylim([-3 0]);
xlabel('1 - Cortical depth (%)'); ylabel('BOLD signal (%)'); axis square; title('PSU')

