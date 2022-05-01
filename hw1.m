clear;
clc;
N=1000;
sequence=randi([0 1],1,N);
figure;
stem(sequence);
title('图1 随机产生0/1数据');
%tabulate(sequence);
% figure;
% h=histogram(sequence);
n1=length(find(sequence==0))/N;
n2=length(find(sequence==1))/N;
figure;
% plot([0,1],[n1,n2]);
b=bar([0,1],[n1,n2]);
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints; %获取 Bar 对象的 XEndPoints 和 YEndPoints 属性
labels1 = string(b.YData); %获取条形末端的坐标
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
title('图2 0/1概率密度');
qpsk= 2*sequence(1:2:end)+sequence(2:2:end);
const=[1+1i,-1+1i,1-1i,-1-1i];
qpsk=genqammod(qpsk,const);
figure;
subplot(2,1,1);
stem(real(qpsk));
title('图3 QPSK调制后的同向分量');
subplot(2,1,2);
stem(imag(qpsk),'r');
title('图3 QPSK调制后的正交分量');