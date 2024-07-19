clc,close all,clear all

load('res_w_P_L2.mat')
load('res_w_P_L10.mat')
load('res_w_P_L20.mat')
load('res_wo_P.mat')


fS = 20; lW = 3;
x2_lim = 3*ones(1,size(tVec_P,2));
dx_lim = [-3;3]*ones(1,size(tVec_P,2));
u_lim = [-0.5;0.5]*ones(1,size(tVec_P,2));
x1 = [DATA_w_Param.x(1,:);DATA2_w_Param.x(1,:);DATA_wo_Param.x(1,:)];
x2 = [DATA_w_Param.x(2,:);DATA2_w_Param.x(2,:);DATA_wo_Param.x(2,:)];
x3 = [DATA_w_Param.x(3,:);DATA2_w_Param.x(3,:);DATA_wo_Param.x(3,:)];
x4 = [DATA_w_Param.x(4,:);DATA2_w_Param.x(4,:);DATA_wo_Param.x(4,:)];
x5 = [DATA_w_Param.x(5,:);DATA2_w_Param.x(5,:);DATA_wo_Param.x(5,:)];
x6 = [DATA_w_Param.x(6,:);DATA2_w_Param.x(6,:);DATA_wo_Param.x(6,:)];
u1 = [DATA_w_Param.u(1,:);DATA2_w_Param.u(1,:);DATA_wo_Param.u(1,:)];
u2 = [DATA_w_Param.u(2,:);DATA2_w_Param.u(2,:);DATA_wo_Param.u(2,:)];
u3 = [DATA_w_Param.u(3,:);DATA2_w_Param.u(3,:);DATA_wo_Param.u(3,:)];

 

figure(1)
subplot(3,1,1);plot(tVec_P,x1,'linewidth',lW);grid on; set(gca,'fontweight','bold','fontsize', fS)
title('angles [rad]')
ylabel('\phi',"FontSize",fS);legend('off')
subplot(3,1,2);plot(tVec_P,x2,'linewidth',lW);grid on; grid on;legend('off')
set(gca,'fontweight','bold','fontsize', fS)
ylabel('\theta')
subplot(3,1,3);plot(tVec_P,x3,'linewidth',lW);grid on; xlabel('time [s]');ylabel('\psi')
set(gca,'fontweight','bold','fontsize', fS)
legend('off')
figure(2)
subplot(3,1,1);plot(tVec_P,x4,'linewidth',lW);grid on; hold on; legend('off')
set(gca,'fontweight','bold','fontsize', fS)
title('angular velocities [rad/s]')
ylabel('\omega_1')
subplot(3,1,2);plot(tVec_P,x5,'linewidth',lW);grid on; set(gca,'fontweight','bold','fontsize', fS); legend('off')
ylabel('\omega_2')
subplot(3,1,3);plot(tVec_P,x6,'linewidth',lW);grid on; set(gca,'fontweight','bold','fontsize', fS)
xlabel('time [s]');ylabel('\omega_3')
legend('off')
figure(3)
subplot(3,1,1);plot(tVec_P,u1,'linewidth',lW);grid on; hold on; plot(tVec_P,u_lim,'--k','LineWidth',lW);legend('off')
title('input [N/kg]')
set(gca,'fontweight','bold','fontsize', fS)
ylabel('u1')
subplot(3,1,2);plot(tVec_P,u2,'linewidth',lW);grid on; hold on; plot(tVec_P,u_lim,'--k','LineWidth',lW);legend('off')
set(gca,'fontweight','bold','fontsize', fS)
ylabel('u2')
subplot(3,1,3);plot(tVec_P,u3,'linewidth',lW);grid on; hold on; plot(tVec_P,u_lim,'--k','LineWidth',lW);
set(gca,'fontweight','bold','fontsize', fS)
xlabel('time [s]');ylabel('u3')
legend('off')
 tComp = [DATA_w_Param.tComp(1,:);DATA2_w_Param.tComp(1,:);DATA_wo_Param.tComp(1,:)];


 
 
for i = 1:size(tComp,1)
    avg_tComp(i) = mean(tComp(i,:));
    sum_tComp(i) = sum(tComp(i,:));
end


avg_tComp
sum_tComp

cost_fcn_res




