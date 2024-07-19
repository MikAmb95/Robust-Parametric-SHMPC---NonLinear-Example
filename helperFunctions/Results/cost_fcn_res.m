
% clearvars
% clc,close all,

Q = diag([1,1,1,0.01,0.01,0.01]);
R = diag([1,1,1])*0.01;
x_ref = [0 0 0 0 0 0]';


xL2 = DATA_w_Param.x; uL2 = DATA_w_Param.u;

xL10 = DATA2_w_Param.x; uL10 = DATA2_w_Param.u;

xL20 = DATA3_w_Param.x; uL20 = DATA3_w_Param.u;

xL0 = DATA_wo_Param.x; uL0 = DATA_wo_Param.u;



J = zeros(4,1);

for i = 1:size(xL2,2)
    J(1) = J(1) + (xL2(:,i)-x_ref)'*Q*(xL2(:,i)-x_ref) + uL2(:,i)'*R*uL2(:,1);
    J(2) = J(2) + (xL10(:,i)-x_ref)'*Q*(xL10(:,i)-x_ref) + uL10(:,i)'*R*uL10(:,1);
    J(3) = J(3) + (xL20(:,i)-x_ref)'*Q*(xL20(:,i)-x_ref) + uL20(:,i)'*R*uL20(:,1);
    J(4) = J(4) + (xL0(:,i)-x_ref)'*Q*(xL0(:,i)-x_ref) + uL0(:,i)'*R*uL0(:,1);
end   

J'*1e-4

sum_nVar_vec = [sum(nVar) sum(nVar2) sum(nVar3) sum(nVar_wo)]

% nVar_vec = [nVar nVar2 nVar3 nVar_wo];
% 
% 
% plot(tVec_P*2,nVar_vec,'linewidth',lW)
% set(gca,'fontweight','bold','fontsize', fS)
% xlabel('time-step k');ylabel('# Decision Variables'); grid on