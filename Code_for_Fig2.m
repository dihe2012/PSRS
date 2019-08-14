%Code_for_Fig2.m, (parallel stochastic resonance)
%计算SR前后的信噪比，分析信噪比增益

clear all;
close all;
clc;
%假设信号服从零均值正态分布,方差为ds^2,噪声服从零均值正态分布，方差为dn^2 
%阵元个数
antenna_number=[1,8,64,256,512,1024,4096,8192];
%样本长度
N=200;
%假定AWGN方差为1,改变信号方差进行仿真,dn为AWGN的标准差
dn=1;

montecarlo_time=1000;

for index3=1:8
    
    N_antenna=antenna_number(index3);

    index3/8,
    
for index2=1:31
    
temp2=0;

for montecarlo=1:montecarlo_time

SNR=-31+index2;    %SNR从-30dB开始增加，间隔1dB

input_SNR(index2)=SNR;


Rsn1=power(10,SNR/10);

ds=sqrt(dn*dn*Rsn1);
t=1:N;

sn=sqrt(2)*ds*cos(2.0*pi*0.02*t);

carrier=cos(2.0*pi*0.02*t);


for index1=1:N_antenna
    
    nn(index1,:)=dn*randn(1,N);
    xn(index1,:)=sn+nn(index1,:);
    
    %compute the SNR(dB) at each receiving antenna
    estimate_amplitude_of_signal(index1)=2*mean(xn(index1,:).*carrier);
    estimate_signal_power(index1)=0.5*power(estimate_amplitude_of_signal(index1),2);
    estimate_total_power(index1)=mean(xn(index1,:).*xn(index1,:));
    estimate_noise_power(index1)=estimate_total_power(index1)-estimate_signal_power(index1);
    estimate_receiving_SNR_dB(index1)=10*log10(estimate_signal_power(index1)/estimate_noise_power(index1));

end

%mean_estimate_amplitude_of_signal=mean(estimate_amplitude_of_signal);
mean(estimate_receiving_SNR_dB);

%introduce the SR system, the input signal is xn
for index1=1:N_antenna
    std_xn(index1)=std(xn(index1,:));
    std_nn(index1)=std(nn(index1,:));
end

%Normalization
for index1=1:N_antenna
    xn2(index1,:)=xn(index1,:)./std_xn(index1);
    nn2(index1,:)=nn(index1,:)./std_nn(index1);
end

delta_t=0.0195;

for index1=1:N_antenna
    x(index1,1)=2.0*rand-1.0;
    y(index1,1)=x(index1,1);
    for i=2:N
        x(index1,i)=x(index1,i-1)+delta_t*(x(index1,i-1)-power(x(index1,i-1),3)+70.0*xn2(index1,i-1));
        y(index1,i)=y(index1,i-1)+delta_t*(y(index1,i-1)-power(y(index1,i-1),3)+70.0*nn2(index1,i-1));
    end
    
    %compute the SNR(dB) of the SR output system
    estimate_amplitude_of_x(index1)=2*mean(x(index1,:).*carrier);
    estimate_signal_power_in_x(index1)=0.5*power(estimate_amplitude_of_x(index1),2);
    estimate_total_power_in_x(index1)=mean(x(index1,:).*x(index1,:));
    estimate_noise_power_in_x(index1)=estimate_total_power_in_x(index1)-estimate_signal_power_in_x(index1);
    estimate_output_SNR_dB(index1)=10*log10(estimate_signal_power_in_x(index1)/estimate_noise_power_in_x(index1));
    
end

x_mean=mean(x,1);
y_mean=mean(y,1);
xn_mean=mean(xn,1);

%compute the SNR(dB) of the SR average output system
    estimate_amplitude_of_x_average=2*mean(x_mean.*carrier);
    estimate_signal_power_in_x_average=0.5*power(estimate_amplitude_of_x_average,2);
    estimate_total_power_in_x_average=mean(x_mean.*x_mean);
    estimate_noise_power_in_x_average=estimate_total_power_in_x_average-estimate_signal_power_in_x_average;
    estimate_output_average_SNR_dB=10*log10(estimate_signal_power_in_x_average/estimate_noise_power_in_x_average);

    
    mean(estimate_output_SNR_dB);
    output_SNR(index2)=mean(estimate_output_average_SNR_dB);

    multi_output_SNR(index3,index2)=output_SNR(index2);

    temp2=temp2+multi_output_SNR(index3,index2);

end

%求平均
sum_output_SNR(index3,index2)=temp2/montecarlo_time;

end

end


% 对 sum_output_SNR 的结果进行差值处理
s=size(sum_output_SNR);
ind=find(~isnan(sum_output_SNR));
[i j]=ind2sub(s,ind);
v=sum_output_SNR(ind);
[ii jj]=ndgrid(1:s(1),1:s(2));
ib=griddata(i,j,v,ii,jj);
sum_output_SNR_interpolation=ib;

%补点
sum_output_SNR_interpolation(8,1)=sum_output_SNR_interpolation(8,3);
sum_output_SNR_interpolation(8,2)=sum_output_SNR_interpolation(8,3);

%画图
figure(1);
plot(input_SNR,sum_output_SNR_interpolation(1,:),'b*-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(input_SNR,sum_output_SNR_interpolation(2,:),'rs-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(input_SNR,sum_output_SNR_interpolation(3,:),'gd-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(input_SNR,sum_output_SNR_interpolation(4,:),'y^-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(input_SNR,sum_output_SNR_interpolation(5,:),'kv-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(input_SNR,sum_output_SNR_interpolation(6,:),'mo-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(input_SNR,sum_output_SNR_interpolation(7,:),'cx-','LineWidth',1.5,'MarkerSize',8);
hold on;
plot(input_SNR,sum_output_SNR_interpolation(8,:),'b+-','LineWidth',1.5,'MarkerSize',8);
hold on;
legend('PSRS Unit Number=1','PSRS Unit Number=8','PSRS Unit Number=64','PSRS Unit Number=256',...
       'PSRS Unit Number=512','PSRS Unit Number=1024','PSRS Unit Number=4096','PSRS Unit Number=8192','location','southeast');
xlabel('Input SNR (dB)');
ylabel('Output SNR (dB)');
grid on;
hold on;
