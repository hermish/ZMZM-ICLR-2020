clc;close all;clear all;
k = [7,9,11,15,19,23,27,29,33,35,37,41];
p = [25,50,75,100,125,150,175,200,225,250,275,300];
c_1 = polyfit(k,p,1);
c_2 = polyfit(k,10*sqrt(p),1);
c_3 = polyfit(k,25*log(p),1);

figure;
hold on;
plot(k,p,'r+');
plot(k,c_1(1)*k+c_1(2),'r');
plot(k,10*sqrt(p),'b*');
plot(k,c_2(1)*k+c_2(2),'b');
plot(k,25*log(p),'go');
plot(k,c_3(1)*k+c_3(2),'g');
legend('data','data-regression','sqrt{data}','sqrt{data}-regression',...
    'log{data}','log{data}-regression');

save('c.mat','c_1');