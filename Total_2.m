clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%发动机外特性曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_min=600;
n_max=4000;
n=n_min:1:n_max;
Tq=-19.313+295.27*(n/1000)-165.44*((n/1000).^2)+40.874*((n/1000).^3)-3.8445*((n/1000).^4);
plot(n,Tq);
hold on;
title('发动机外特性曲线');
ylabel('发动机转矩(N·m)');
xlabel('发动机转速(r/min)');
grid on;

%%%%求最大车速
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%驱动力曲线%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i0=5.83;%主减速器传动比
nT=0.85;%传动系机械效率
r=0.367;%车辆半径
ig5_1=5.56;
ig5_2=2.769;
ig5_3=1.644;
ig5_4=1.00;
ig5_5=0.793;
Ft5_1=Tq.*ig5_1*i0*nT/r;%驱动力
Ft5_2=Tq.*ig5_2*i0*nT/r;
Ft5_3=Tq.*ig5_3*i0*nT/r;
Ft5_4=Tq.*ig5_4*i0*nT/r;
Ft5_5=Tq.*ig5_5*i0*nT/r;
ua5_1=0.377*r.*n.*(1./ig5_1).*(1/i0);%车速
ua5_2=0.377*r.*n.*(1./ig5_2).*(1/i0);
ua5_3=0.377*r.*n.*(1./ig5_3).*(1/i0);
ua5_4=0.377*r.*n.*(1./ig5_4).*(1/i0);
ua5_5=0.377*r.*n.*(1./ig5_5).*(1/i0);

CD_A=2.77;%空气阻力系数*迎风面积
f=0.013;%滚动阻力系数
%Wk=1800;%整车装备质量
g=9.8;
ua=0:0.1:120;%车速

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%空载时的最大爬坡度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fw=CD_A*ua.^2/21.15;%空气阻力
Fw5_1=CD_A*ua5_1.^2/21.15;%空气阻力
Fw5_2=CD_A*ua5_2.^2/21.15;
Fw5_3=CD_A*ua5_3.^2/21.15;
Fw5_4=CD_A*ua5_4.^2/21.15;
Fw5_5=CD_A*ua5_5.^2/21.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%满载时的最高车速%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wm=3880;%总质量
Gm=Wm*g;
Ffm=Gm*f;
Ftm=Ffm+Fw;%行驶阻力
Ftm5_1=Ffm+Fw5_1;
Ftm5_2=Ffm+Fw5_2;
Ftm5_3=Ffm+Fw5_3;
Ftm5_4=Ffm+Fw5_4;
Ftm5_5=Ffm+Fw5_5;
for i=1:3401
    if(Ft5_5(i)<Ftm5_5(i))
        ua_F_max=ua5_5(i-1);%满载最高车速
        break;
    end
end
disp('5档变速箱满载最高车速');
disp(ua_F_max);
figure;
plot(ua5_1,Ft5_1,'g');
hold on;
plot(ua5_2,Ft5_2,'m');
hold on;
plot(ua5_3,Ft5_3,'b');
hold on;
plot(ua5_4,Ft5_4,'r');
hold on;
plot(ua5_5,Ft5_5,'k');
hold on;
plot(ua5_5,Ftm5_5,'y');
hold on;
plot(ua5_5(i-1),Ftm5_5(i-1),'.');%画交点
hold on;
title('5档变速箱驱动力-行驶阻力平衡图');
ylabel('Ft(N)');
xlabel('ua(km/h)');
legend('Ft1','Ft2','Ft3','Ft4','Ft5','Ff+Fw');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%满载时的最大爬坡度%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D5_1=(Ft5_1-Fw5_1)/Gm;
D5_2=(Ft5_2-Fw5_2)/Gm;
D5_3=(Ft5_3-Fw5_3)/Gm;
D5_4=(Ft5_4-Fw5_4)/Gm;
D5_5=(Ft5_5-Fw5_5)/Gm;
a4_1=asin((D5_1-f*sqrt((1-D5_1.^2+f*f)))/(1+f*f));
a4_2=asin((D5_2-f*sqrt((1-D5_2.^2+f*f)))/(1+f*f));
a4_3=asin((D5_3-f*sqrt((1-D5_3.^2+f*f)))/(1+f*f));
a4_4=asin((D5_4-f*sqrt((1-D5_4.^2+f*f)))/(1+f*f));
a4_5=asin((D5_5-f*sqrt((1-D5_5.^2+f*f)))/(1+f*f));
i4_1=tan(a4_1);
i4_2=tan(a4_2);
i4_3=tan(a4_3);
i4_4=tan(a4_4);
i4_5=tan(a4_5);
figure;
plot(ua5_1,i4_1,'g');
hold on;
plot(ua5_2,i4_2,'m');
hold on;
plot(ua5_3,i4_3,'b');
hold on;
plot(ua5_4,i4_4,'r');
hold on;
plot(ua5_5,i4_5,'k');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%满载时的最大附着率%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[i4_1_max,ua4_1_max]=max(i4_1);%求最大值及其位置
disp('5档变速箱满载时的最大爬坡度的附着率');
disp(i4_1_max);
plot(ua5_1(ua4_1_max),i4_1_max,'.');%画点
hold on;
title('5档变速箱满载时的最大爬坡度');
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
%Iw=2*Iw1+4*Iw2;
Iw=Iw1+Iw2;
If=0.218;%飞轮转动惯量
delta5_1=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig5_1^2*i0^2*nT)/(r*r);
delta5_2=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig5_2^2*i0^2*nT)/(r*r);
delta5_3=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig5_3^2*i0^2*nT)/(r*r);
delta5_4=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig5_4^2*i0^2*nT)/(r*r);
delta5_5=1+(1/Wm)*(Iw/r^2)+(1/Wm)*(If*ig5_5^2*i0^2*nT)/(r*r);
aj5_1=(1/(delta5_1*Wm))*(Ft5_1-(Ffm+Fw5_1));
aj5_2=(1/(delta5_2*Wm))*(Ft5_2-(Ffm+Fw5_2));
aj5_3=(1/(delta5_3*Wm))*(Ft5_3-(Ffm+Fw5_3));
aj5_4=(1/(delta5_4*Wm))*(Ft5_4-(Ffm+Fw5_4));
aj5_5=(1/(delta5_5*Wm))*(Ft5_5-(Ffm+Fw5_5));
figure;
plot(ua5_1,aj5_1,'g');
hold on;
plot(ua5_2,aj5_2,'m');
hold on;
plot(ua5_3,aj5_3,'b');
hold on;
plot(ua5_4,aj5_4,'r');
hold on;
plot(ua5_5,aj5_5,'k');
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
plot(ua5_1,ajd5_1,'g');
hold on;
plot(ua5_2,ajd5_2,'m');
hold on;
plot(ua5_3,ajd5_3,'b');
hold on;
plot(ua5_4,ajd5_4,'r');
hold on;
plot(ua5_5,ajd5_5,'k');
hold on;
ylim([0,20]);%限制Y轴取值，防止图像变形
xlim([0,ua_F_max]);%限制X轴取值，防止图像变形
title('五档变速箱满载时的加速度倒数曲线');
ylabel('1/aj');
xlabel('ua');
legend('1/aj1','1/aj2','1/aj3','1/aj4','1/aj5');
grid on;
%axis([0 100 0 20]);
%%%%%%%%%%%%%%%%%%%%%%%%满载时用2挡起步加速到70km/h的最短加速时间%%%%%%%%%%%%%
t2=trapz(ua5_2,ajd5_2);
for i3=1:length(ua5_3)
    if (ua5_3(i3)>ua5_2(end))
        t3=trapz(ua5_3(i3:end),ajd5_3(i3:end));
        break;
    end
end
for i4=1:length(ua5_4)
    if (ua5_4(i4)>ua5_3(end))
        for j4=1:length(ua5_4)
            if (ua5_4(j4)>=70)
                t4=trapz((ua5_4(i4:j4)),(ajd5_4(i4:j4)));
                break;
            end
        end
        break
    end
end
t=(t2+t3+t4)/3.6;
disp('满载时用2挡起步加速到70km/h的最短加速时间为');
disp(t);
%%%%%%%%%%%%%%%%%%%80km/h分别采用 4档和5档时的百公里油耗%%%%%%%%%%%%%%%%%%%%%
pg=7;
va_80=80;
n_80_4d=va_80.*i0.*ig5_4/r/0.377;
n_80_5d=va_80.*i0.*ig5_5/r/0.377;

for j=1:length(ua5_4)
    if (ua5_4(j)>va_80)
        Ff_80_4d=Ftm5_4(j);
        break;
    end
end
for j=1:length(ua5_5)
    if (ua5_5(j)>va_80)
        Ff_80_5d=Ftm5_5(j);
        break;
    end
end
Pe_4d=Ff_80_4d.*va_80/3.6/nT/1000;
Pe_5d=Ff_80_5d.*va_80/3.6/nT/1000;
disp('n_80_4d,Pe_4d');
disp(n_80_4d);
disp(Pe_4d);
disp('n_80_5d,Pe_5d');
disp(n_80_5d);
disp(Pe_5d);
b_80_4d=260;
b_80_5d=250;
disp('b_80_4d');
disp(b_80_4d);
disp('b_80_5d');
disp(b_80_5d);
Qs_80_4d=Pe_4d*b_80_4d/(1.02*va_80*pg);
Qs_80_5d=Pe_5d*b_80_5d/(1.02*va_80*pg);
disp('Qs_80_4d');
disp(Qs_80_4d);
disp('Qs_80_5d');
disp(Qs_80_5d);
%%%%%%%%%%%%5档时，车速分别为40、50、60、70和80km/h的燃油消耗量%%%%%%%%%%%%%%
va=[40, 50, 60, 70, 80];
n_5d=va.*i0.*ig5_5/r/0.377;
for i=1:length(va)
    for j=1:length(ua5_5)
        if (ua5_5(j)>va(i))
            Ff(i)=Ftm5_5(j);
            break;
        end
    end
end
Pe_5d=Ff.*va/3.6/nT/1000;
disp('n_5d,Pe_5d');
disp(n_5d);
disp(Pe_5d);
b=[290,280,260,250,250];
disp('查图取得燃料消耗率b');
disp(b);
Qt_5d=Pe_5d.*b/(367.1*pg);
disp('5档百公里油耗Qt');
disp(Qt_5d);
