%  Title :��BP�㷨���״�Ŀ����г���
%  Date��2019-12-10
%  Author:  ������
clear all;clc;close all;
fc=77e9;
j=sqrt(-1);
B=1e9;   
c=3e8;
Fs=1e9;
T=500e-9;        % PRT 
PRF=1/T;
K=B/T;
R=[6,6;7,7;8,10;9,10;4,6;8,8];            % ����Ŀ����������ͺ�������
R_zong=[sqrt(sum(R(1,:).^2)),sqrt(sum(R(2,:).^2)),sqrt(sum(R(3,:).^2))];

% Info_Target=[
%      5 0; 5 -20;   
% %      110 0;  110 20; 110 -20;
% %      150 0;   150 20;  150 -20;                     
%     ];        %   position of target ;
% R=
% R_zong=[];
% for i=1:size(Info_Target,1)
%     R(i,:)=[Info_Target(i,1)*cosd(Info_Target(i,2)),Info_Target(i,1)*sind(Info_Target(i,2))];
%     R_temp=sqrt(sum(R(i,:).^2));
%     R_zong=[R_zong,R_temp];
% end
lambda=c/fc;                                 %  �źŲ���
d=lambda/2;
SNR=10;                           %       �����
L=128;                                    %  ���У�L��ʾ��Ԫ����
N=Fs*T;                                 %   ��ʾ��������
t_n=[1:N]/Fs;              %  ��ʱ�䷽��
length_Antenna=d*L;
X_Antenna=linspace(-length_Antenna/2,length_Antenna/2,L);
Q=size(R,1);
IF_data=zeros(L,N);      %   ��Ƶ�ز�����
Pc_data=zeros(size(IF_data));
DeltaR=0.1425;
Range=DeltaR*[0:N-1];
% index=1
Nsinc=11;
f=[0:N-1]/N*Fs;
% f=0:
for l=0:L-1
%     Tau_Range(index)=sqrt((R(1,1)-X_Antenna(l+1))^2+R(1,2)^2);
    x=zeros(Q,N);
    for q=1:Q
        Tau_Range(q)=sqrt((R(q,1)-X_Antenna(l+1))^2+R(q,2)^2);
        x(q,:)=exp(2*1j*pi*(((2*K*Tau_Range(q))/c)*t_n+2*fc*Tau_Range(q)/c));
%         pc_temp(q,:)=1e13*sinc(B*(f-2*K*Tau_Range(q)/c))*exp(1j*2*pi*2*fc*Tau_Range(q)/c);
    end
%     pc(l+1,:)=sum(pc_temp);
    IF_data(l+1,:)=sum(x);
end
dn=awgn(IF_data,SNR);
% pc=pc+awgn(pc,SNR);
Pc_data=fft(IF_data,N,2);           %  ����ѹ�����ź�
figure;
plot(f,abs(Pc_data(2,:)))
% figure
% plot(f,abs(pc(100,:)))
%%   ��Ŀ��������г���
V =100;   %   �������
H=100;    %   �������
V_grid=2+linspace(0,1,V)*10;
H_grid=2+linspace(0,1,H)*10;
Matrix_Imag=zeros(V,H);
for v=1:V
    for h=1:H
        Index=1;
        for l=1:L
            t_mn(Index)=2*sqrt((V_grid(v)-X_Antenna(l))^2+H_grid(h)^2)/c;   %  �����ʱ����ʱ����
%             if v==22 & h==24
%                 sel=1
%             end
            [out,x_out] = sinc_interp(Pc_data(l,:),f,K*t_mn(Index),Nsinc,1);    %  ������Ĳ�ֵ
            signal_echo1(Index,:)=out*exp(-1j*2*pi*fc*(t_mn(Index)));
            Index=Index+1;

        end
        Matrix_Imag(v,h)=abs(sum(signal_echo1));

    end
end

figure;imagesc(V_grid,H_grid,Matrix_Imag)    %  ����V ��ʾ���Ǻ����꣬��H��ʾ������
set(gca,'YDir','normal');colormap('gray');
title('BP Imaging Algorithm For mmWave Radar');
xlabel('������');ylabel('��λ��');
figure
contour(Matrix_Imag,15)