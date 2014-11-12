clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%发动机外特性曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=600:1:4000;
Tq=-19.313+295.27*(n/1000)-165.44*((n/1000).^2)+40.874*((n/1000).^3)-3.8445*((n/1000).^4);
plot(n,Tq);
hold on;
title('发动机外特性曲线');
ylabel('发动机转矩(N·m)');
xlabel('发动机转速(r/min)');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%五档变速箱%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%驱动力曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i0=5.83;
nT=0.85;
r=0.367;
ig4_1=5.56;
ig4_2=2.769;
ig4_3=1.644;
ig4_4=1.00;
ig4_5=0.793;
Ft4_1=Tq.*ig4_1*i0*nT/r;
Ft4_2=Tq.*ig4_2*i0*nT/r;
Ft4_3=Tq.*ig4_3*i0*nT/r;
Ft4_4=Tq.*ig4_4*i0*nT/r;
Ft4_5=Tq.*ig4_5*i0*nT/r;
ua4_1=0.377*r.*n.*(1./ig4_1).*(1/i0);
ua4_2=0.377*r.*n.*(1./ig4_2).*(1/i0);
ua4_3=0.377*r.*n.*(1./ig4_3).*(1/i0);
ua4_4=0.377*r.*n.*(1./ig4_4).*(1/i0);
ua4_5=0.377*r.*n.*(1./ig4_5).*(1/i0);
figure;
plot(ua4_1,Ft4_1,'g');
hold on;
plot(ua4_2,Ft4_2,'m');
hold on;
plot(ua4_3,Ft4_3,'b');
hold on;
plot(ua4_4,Ft4_4,'r');
hold on;
plot(ua4_5,Ft4_5,'r');
hold on;
title('五档变速箱驱动力-行驶阻力平衡图');
ylabel('Ft(N)');
xlabel('ua(km/h)');
legend('Ft1','Ft2','Ft3','Ft4','Ft5');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%空载时的最高车速%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CD=0.68;
%A=4.07;
CD_A=2.77;
f=0.013;
ua=0:0.1:120;
Fw=CD_A*ua.^2/21.15;
Wk=1800;
g=9.8;
Gk=Wk*g;
Ffk=Gk*f;
Ftk=Ffk+Fw;
figure;
plot(ua4_1,Ft4_1,'g');
hold on;
plot(ua4_2,Ft4_2,'m');
hold on;
plot(ua4_3,Ft4_3,'b');
hold on;
plot(ua4_4,Ft4_4,'r');
hold on;
plot(ua4_5,Ft4_5,'r');
hold on;
plot(ua,Ftk,'k');
hold on;
title('五档变速箱空载时的最高车速');
ylabel('Ft(N)');
xlabel('ua(km/h)');
legend('Ft1','Ft2','Ft3','Ft4','Ft5','Ff+Fw');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%空载时的最大爬坡度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fw4_1=CD_A*ua4_1.^2/21.15;
Fw4_2=CD_A*ua4_2.^2/21.15;
Fw4_3=CD_A*ua4_3.^2/21.15;
Fw4_4=CD_A*ua4_4.^2/21.15;
Fw4_5=CD_A*ua4_5.^2/21.15;
D4_1=(Ft4_1-Fw4_1)/Gk;
a4_1=asin((D4_1-f*sqrt((1-D4_1.^2+f*f)))/(1+f*f));
a4_2=asin((Ft4_2-(Ffk+Fw4_2))/Gk);
a4_3=asin((Ft4_3-(Ffk+Fw4_3))/Gk);
a4_4=asin((Ft4_4-(Ffk+Fw4_4))/Gk);
a4_5=asin((Ft4_5-(Ffk+Fw4_5))/Gk);
i4_1=tan(a4_1);
i4_2=tan(a4_2);
i4_3=tan(a4_3);
i4_4=tan(a4_4);
i4_5=tan(a4_5);
figure;
plot(ua4_1,i4_1,'g');
hold on;
plot(ua4_2,i4_2,'m');
hold on;
plot(ua4_3,i4_3,'b');
hold on;
plot(ua4_4,i4_4,'r');
hold on;
plot(ua4_5,i4_5,'r');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%空载时的最大附着率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[i4_1_max,ua4_1_max]=max(i4_1);%求最大值及其位置
plot(ua4_1(ua4_1_max),i4_1_max,'.');%画点
hold on;
title('五档变速箱空载时的最大爬坡度');
ylabel('i');
xlabel('ua(km/h)');
legend('i1','i2','i3','i4','i5');
grid on;
%由于缺少数据，未计算附着系数
%L=3.2;%轴距
%a=1.947;%质心至前轴距离（满载）
%hg=0.9;%质心高（满载）
%q=i4_1_max;%包含加速阻力在内的等效坡度
%CO2=q/(a/L+hg/L*q);
disp('空载时克服最大坡度的附着率');
disp('由于缺少数据，未计算附着系数');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%空载时加速度倒数曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iw1=1.798;%二前轮转动惯量
Iw2=3.598;%四后轮转动惯量
Iw=2*Iw1+4*Iw2;
If=0.218;%飞轮转动惯量
delta4_1=1+(1/Wk)*(Iw/r^2)+(1/Wk)*(If*ig4_1^2*i0^2*nT)/(r*r);
delta4_2=1+(1/Wk)*(Iw/r^2)+(1/Wk)*(If*ig4_2^2*i0^2*nT)/(r*r);
delta4_3=1+(1/Wk)*(Iw/r^2)+(1/Wk)*(If*ig4_3^2*i0^2*nT)/(r*r);
delta4_4=1+(1/Wk)*(Iw/r^2)+(1/Wk)*(If*ig4_4^2*i0^2*nT)/(r*r);
delta4_5=1+(1/Wk)*(Iw/r^2)+(1/Wk)*(If*ig4_5^2*i0^2*nT)/(r*r);
aj4_1=(1/(delta4_1*Wk))*(Ft4_1-(Ffk+Fw4_1));
aj4_2=(1/(delta4_2*Wk))*(Ft4_2-(Ffk+Fw4_2));
aj4_3=(1/(delta4_3*Wk))*(Ft4_3-(Ffk+Fw4_3));
aj4_4=(1/(delta4_4*Wk))*(Ft4_4-(Ffk+Fw4_4));
aj4_5=(1/(delta4_5*Wk))*(Ft4_5-(Ffk+Fw4_5));
figure;
plot(ua4_1,aj4_1,'g');
hold on;
plot(ua4_2,aj4_2,'m');
hold on;
plot(ua4_3,aj4_3,'b');
hold on;
plot(ua4_4,aj4_4,'r');
hold on;
plot(ua4_5,aj4_5,'r');
hold on;
title('五档变速箱空载时的加速度曲线');
ylabel('aj(m/g^2)');
xlabel('ua(km/h)');
legend('aj1','aj2','aj3','aj4','aj5');
grid on;
ajd4_1=1./aj4_1;
ajd4_2=1./aj4_2;
ajd4_3=1./aj4_3;
ajd4_4=1./aj4_4;
ajd4_5=1./aj4_5;
figure;
plot(ua4_1,ajd4_1,'g');
hold on;
plot(ua4_2,ajd4_2,'m');
hold on;
plot(ua4_3,ajd4_3,'b');
hold on;
plot(ua4_4,ajd4_4,'r');
hold on;
ylim([0,4]);%限制Y轴取值，防止图像变形
plot(ua4_5,ajd4_5,'r');
hold on;
title('五档变速箱空载时的加速度倒数曲线');
ylabel('1/aj');
xlabel('ua');
legend('1/aj1','1/aj2','1/aj3','1/aj4','1/aj5');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%空载时用2挡起步加速到70km/h的最短加速时间%%%%%%%%%%%%%
%t1=quad(@intg4_1,ua4_1(1,1),ua4_1(1,3401));
%t2=quad(@intg4_2,ua4_1(1,3401),ua4_2(1,3401));
%t3=quad(@intg4_3,ua4_2(1,3401),ua4_3(1,3401));
%t4=quad(@intg4_4,ua4_3(1,3401),70);
t1=trapz(ua4_1,ajd4_1);
t2=trapz(ua4_2,ajd4_2);
t3=trapz(ua4_3,ajd4_3);
for i=1:length(ua4_4)
    if (ua4_4(i)>=70)
        ua4_4_70=i;
        break;
    end
end
if (isempty(ua4_4_70))
    disp('无法达到70km/h');
end
t4=trapz(ua4_4(1:ua4_4_70),ajd4_4(1:ua4_4_70));
t=(t1+t2+t3+t4)/3.6;
disp('空载时用2挡起步加速到70km/h的最短加速时间为');
disp(t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%满载时的最高车速%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wm=3880;    %总质量
Gm=Wm*g;
Ffm=Gm*f;
Ftm=Ffm+Fw;
figure;
plot(ua4_1,Ft4_1,'g');
hold on;
plot(ua4_2,Ft4_2,'m');
hold on;
plot(ua4_3,Ft4_3,'b');
hold on;
plot(ua4_4,Ft4_4,'r');
hold on;
plot(ua4_5,Ft4_5,'r');
hold on;
plot(ua,Ftm,'k');
hold on;
title('五档变速箱满载时的最高车速');
ylabel('Ft(N)');
xlabel('ua(km/h)');
legend('Ft1','Ft2','Ft3','Ft4','Ft5','Ff+Fw');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%满载时的最大爬坡度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D4_1=(Ft4_1-Fw4_1)/Gm;
a4_1=asin((D4_1-f*sqrt((1-D4_1.^2+f*f)))/(1+f*f));
a4_2=asin((Ft4_2-(Ffk+Fw4_2))/Gm);
a4_3=asin((Ft4_3-(Ffk+Fw4_3))/Gm);
a4_4=asin((Ft4_4-(Ffk+Fw4_4))/Gm);
a4_5=asin((Ft4_5-(Ffk+Fw4_5))/Gm);
i4_1=tan(a4_1);
i4_2=tan(a4_2);
i4_3=tan(a4_3);
i4_4=tan(a4_4);
i4_5=tan(a4_5);
figure;
plot(ua4_1,i4_1,'g');
hold on;
plot(ua4_2,i4_2,'m');
hold on;
plot(ua4_3,i4_3,'b');
hold on;
plot(ua4_4,i4_4,'r');
hold on;
plot(ua4_5,i4_5,'r');
hold on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%满载时的最大附着率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[i4_1_max,ua4_1_max]=max(i4_1);%求最大值及其位置
plot(ua4_1(ua4_1_max),i4_1_max,'.');%画点
hold on;
title('五档变速箱满载时的最大爬坡度');
ylabel('i');
xlabel('ua(km/h)');
legend('i1','i2','i3','i4','i5');
grid on;
L=3.2;%轴距
a=1.947;%质心至前轴距离（满载）
hg=0.9;%质心高（满载）
q=i4_1_max;%包含加速阻力在内的等效坡度
CO2=q/(a/L+hg/L*q);
disp('满载时克服最大坡度的附着率');
disp(CO2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%满载时加速度倒数曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iw1=1.798;%二前轮转动惯量
Iw2=3.598;%四后轮转动惯量
Iw=2*Iw1+4*Iw2;
If=0.218;%飞轮转动惯量
delta5_1=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig4_1^2*i0^2*nT)/(r*r);
delta5_2=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig4_2^2*i0^2*nT)/(r*r);
delta5_3=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig4_3^2*i0^2*nT)/(r*r);
delta5_4=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig4_4^2*i0^2*nT)/(r*r);
delta5_5=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig4_5^2*i0^2*nT)/(r*r);
aj5_1=(1/(delta5_1*Wm))*(Ft4_1-(Ffm+Fw4_1));
aj5_2=(1/(delta5_2*Wm))*(Ft4_2-(Ffm+Fw4_2));
aj5_3=(1/(delta5_3*Wm))*(Ft4_3-(Ffm+Fw4_3));
aj5_4=(1/(delta5_4*Wm))*(Ft4_4-(Ffm+Fw4_4));
aj5_5=(1/(delta5_5*Wm))*(Ft4_5-(Ffm+Fw4_5));
figure;
plot(ua4_1,aj5_1,'g');
hold on;
plot(ua4_2,aj5_2,'m');
hold on;
plot(ua4_3,aj5_3,'b');
hold on;
plot(ua4_4,aj5_4,'r');
hold on;
plot(ua4_5,aj5_5,'k');
hold on;
title('五档变速箱满载时的加速度曲线');
ylabel('aj(m/g^2)');
xlabel('ua(km/h)');
legend('aj1','aj2','aj3','aj4','aj5');
grid on;
ajd5_1=1./aj5_1;
ajd5_2=1./aj5_2;
ajd5_3=1./aj5_3;
ajd5_4=1./aj5_4;
ajd5_5=1./aj5_5;
figure;
plot(ua4_1,ajd5_1,'g');
hold on;
plot(ua4_2,ajd5_2,'m');
hold on;
plot(ua4_3,ajd5_3,'b');
hold on;
plot(ua4_4,ajd5_4,'r');
hold on;
plot(ua4_5,ajd5_5,'k');
hold on;
ylim([0,4]);%限制Y轴取值，防止图像变形
title('五档变速箱满载时的加速度倒数曲线');
ylabel('1/aj');
xlabel('ua');
legend('1/aj1','1/aj2','1/aj3','1/aj4','1/aj5');
grid on;
%axis([0 100 0 20]);
%%%%%%%%%%%%%%%%%%%%%%%%满载时用2挡起步加速到70km/h的最短加速时间%%%%%%%%%%%%%
%t1=quad(@intg4_1,ua4_1(1,1),ua4_1(1,3401));
%t2=quad(@intg4_2,ua4_1(1,3401),ua4_2(1,3401));
%t3=quad(@intg4_3,ua4_2(1,3401),ua4_3(1,3401));
%t4=quad(@intg4_4,ua4_3(1,3401),70);
t1=trapz(ua4_1,ajd5_1);
t2=trapz(ua4_2,ajd5_2);
t3=trapz(ua4_3,ajd5_3);
for i=1:length(ua4_4)
    if (ua4_4(i)>=70)
        ua4_4_70=i;
        break;
    end
end
if (isempty(ua4_4_70))
    disp('无法达到70km/h');
end
t4=trapz(ua4_4(1:ua4_4_70),ajd4_4(1:ua4_4_70));
t=(t1+t2+t3+t4)/3.6;
disp('满载时用2挡起步加速到70km/h的最短加速时间为');
disp(t);