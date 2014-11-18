clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����������������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_min=600;
n_max=4000;
n=n_min:1:n_max;
Tq=-19.313+295.27*(n/1000)-165.44*((n/1000).^2)+40.874*((n/1000).^3)-3.8445*((n/1000).^4);
plot(n,Tq);
hold on;
title('����������������');
ylabel('������ת��(N��m)');
xlabel('������ת��(r/min)');
grid on;

%%%%�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i0=5.83;%��������������
nT=0.85;%����ϵ��еЧ��
r=0.367;%�����뾶
ig5_1=5.56;%f_ǰ׺��ֹ����ʱ��������������
ig5_2=2.769;
ig5_3=1.644;
ig5_4=1.00;
ig5_5=0.793;
Ft5_1=Tq.*ig5_1*i0*nT/r;%������
Ft5_2=Tq.*ig5_2*i0*nT/r;
Ft5_3=Tq.*ig5_3*i0*nT/r;
Ft5_4=Tq.*ig5_4*i0*nT/r;
Ft5_5=Tq.*ig5_5*i0*nT/r;
ua5_1=0.377*r.*n.*(1./ig5_1).*(1/i0);%����
ua5_2=0.377*r.*n.*(1./ig5_2).*(1/i0);
ua5_3=0.377*r.*n.*(1./ig5_3).*(1/i0);
ua5_4=0.377*r.*n.*(1./ig5_4).*(1/i0);
ua5_5=0.377*r.*n.*(1./ig5_5).*(1/i0);
figure;
plot(ua5_1,Ft5_1,'g');
hold on;
plot(ua5_2,Ft5_2,'m');
hold on;
plot(ua5_3,Ft5_3,'b');
hold on;
plot(ua5_4,Ft5_4,'r');
hold on;
plot(ua5_5,Ft5_5,'r');
hold on;
title('�嵵������������-��ʻ����ƽ��ͼ');
ylabel('Ft(N)');
xlabel('ua(km/h)');
legend('Ft1','Ft2','Ft3','Ft4','Ft5');
grid on;



CD_A=2.77;%��������ϵ��*ӭ�����
f=0.013;%��������ϵ��
%Wk=1800;%����װ������
g=9.8;
ua=0:0.1:120;%����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ʱ��������¶�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fw=CD_A*ua.^2/21.15;%��������
Fw5_1=CD_A*ua5_1.^2/21.15;%��������
Fw5_2=CD_A*ua5_2.^2/21.15;
Fw5_3=CD_A*ua5_3.^2/21.15;
Fw5_4=CD_A*ua5_4.^2/21.15;
Fw5_5=CD_A*ua5_5.^2/21.15;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ʱ����߳���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wm=3880;%������
Gm=Wm*g;
Ffm=Gm*f;
Ftm=Ffm+Fw;
%%������߳���
% for i=1:length(ua)*10
%     if (Ft5_5(i)<Ftm(ceil(ua5_5(i))))
%         ua_F_max=ua5_5(i);
%         break;
%     end
% end
ua_F_max=97;%������߳���
figure;
plot(ua5_1,Ft5_1,'g');
hold on;
plot(ua5_2,Ft5_2,'m');
hold on;
plot(ua5_3,Ft5_3,'b');
hold on;
plot(ua5_4,Ft5_4,'r');
hold on;
plot(ua5_5,Ft5_5,'r');
hold on;
plot(ua,Ftm,'k');
hold on;
title('�嵵����������ʱ����߳���');
ylabel('Ft(N)');
xlabel('ua(km/h)');
legend('Ft1','Ft2','Ft3','Ft4','Ft5','Ff+Fw');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ʱ��������¶�%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
plot(ua5_5,i4_5,'r');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ʱ���������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[i4_1_max,ua4_1_max]=max(i4_1);%�����ֵ����λ��
plot(ua5_1(ua4_1_max),i4_1_max,'.');%����
hold on;
title('�嵵����������ʱ��������¶�');
ylabel('i');
xlabel('ua(km/h)');
legend('i1','i2','i3','i4','i5');
grid on;
L=3.2;%���
a=1.947;%������ǰ����루���أ�
hg=0.9;%���ĸߣ����أ�
q=i4_1_max;%���������������ڵĵ�Ч�¶�
CO2=q/(a/L+hg/L*q);
disp('����ʱ�˷�����¶ȵĸ�����');
disp(CO2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����ʱ���ٶȵ�������%%%%%%%%%%%%%%%%%%%%%%%%%%%
Iw1=1.798;%��ǰ��ת������
Iw2=3.598;%�ĺ���ת������
Iw=2*Iw1+4*Iw2;
If=0.218;%����ת������
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
title('�嵵����������ʱ�ļ��ٶ�����');
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
ylim([0,20]);%����Y��ȡֵ����ֹͼ�����
xlim([0,ua_F_max]);%����X��ȡֵ����ֹͼ�����
title('�嵵����������ʱ�ļ��ٶȵ�������');
ylabel('1/aj');
xlabel('ua');
legend('1/aj1','1/aj2','1/aj3','1/aj4','1/aj5');
grid on;
%axis([0 100 0 20]);
%%%%%%%%%%%%%%%%%%%%%%%%����ʱ��2���𲽼��ٵ�70km/h����̼���ʱ��%%%%%%%%%%%%%
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
disp('����ʱ��2���𲽼��ٵ�70km/h����̼���ʱ��Ϊ');
disp(t);
%%%%%%%%%%%%%%%%%%%80km/h�ֱ���� 4����5��ʱ�İٹ����ͺ�%%%%%%%%%%%%%%%%%%%%%
va=80;
n_5d=va*1200/50+400;
disp('5��ʱ������ת��(r/min)');
disp(n_5d);
Tq_5d=-19.313+295.27*(n_5d/1000)-165.44*((n_5d/1000).^2)+40.874*((n_5d/1000).^3)-3.8445*((n_5d/1000).^4);%N*m
Pe_5d=Tq_5d*(2*pi/60*n_5d)/1000;%kW
disp('5��ʱ����������(kW)');
disp(Pe_5d);
b=260;
disp('��ͼȡ��ȼ��������b');
disp(b);
pg=7;
Qs_5d=Pe_5d*b/(1.02*80*pg);
disp('5���ٹ����ͺ�Qs');
disp(Qs_5d);
n_4d=n_5d*ig5_5/ig5_4;
disp('4��ʱ������ת��(r/min)');
disp(n_4d);
Tq_4d=-19.313+295.27*(n_4d/1000)-165.44*((n_4d/1000).^2)+40.874*((n_4d/1000).^3)-3.8445*((n_4d/1000).^4);%N*m
Pe_4d=Tq_4d*(2*pi/60*n_4d)/1000;%kW
disp('4��ʱ����������(kW)');
disp(Pe_4d);
Qs_4d=Pe_4d*b/(1.02*80*pg);
disp('4���ٹ����ͺ�Qs');
disp(Qs_4d);
%%%%%%%%%%%%5��ʱ�����ٷֱ�Ϊ40��50��60��70��80km/h��ȼ��������%%%%%%%%%%%%%%
va=[40, 50, 60, 70, 80];
n_5d=va.*1200/50+400;
Tq_5d=-19.313+295.27*(n_5d./1000)-165.44*((n_5d./1000).^2)+40.874*((n_5d./1000).^3)-3.8445*((n_5d./1000).^4);%N*m
Pe_5d=(Tq_5d).*(2*pi/60*(n_5d))/1000;%kW
b=260;
disp('��ͼȡ��ȼ��������b');
disp(b);
%pg=7;
Qs_5d=Pe_5d.*b/(1.02*80*pg);
disp('5���ٹ����ͺ�Qs');
disp(Qs_5d);