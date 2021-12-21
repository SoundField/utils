function [CC,k]=cst(f,H,cp,cs,rhob)
close all;clc;pi=3.141592653;eps=10^-4;   %%eps是解c的精度，提高精度更改eps
c1=1500;
%cp=1550;cs=800;
xl=1;CC=zeros(1,1);  %%c1是水中声速，cp是纵波声速，cs是横波声速
rhob1=1000;rhob2=rhob;   %%rhob1是水密度，rhob2是弹性体密度
%  y=tanh(2*pi*f*H/c*sqrt(1-c^2/c1^2))*(rhob1*c^4*sqrt(1-c^2/cp^2))/ ...
% (rhob2*cs^4*sqrt(1-c^2/c1^2))-(2-c^2/cs^2)^2+4*sqrt(1-(c/cs)^2)*sqrt(1-c^2/cp^2);
for c=10:0.1:min(min(c1,cp),cs)-2
    bot=c;top=bot+0.1;cc=(bot+top)/2;    %%bot是二分法单元下限，top是上限，cc是中值
    while (abs(top-bot)>eps)           %%当精度不足eps时进行计算
    aaa=tanh(2*pi*f*H/top*sqrt(1-top^2/c1^2))*(rhob1*top^4*sqrt(1-top^2/cp^2))/ ...
(rhob2*cs^4*sqrt(1-(top/c1)^2))+(2-top^2/cs^2)^2-4*sqrt(1-(top/cs)^2)*sqrt(1-(top/cp)^2);
    bbb=tanh(2*pi*f*H/bot*sqrt(1-bot^2/c1^2))*(rhob1*bot^4*sqrt(1-bot^2/cp^2))/ ...
(rhob2*cs^4*sqrt(1-bot^2/c1^2))+(2-bot^2/cs^2)^2-4*sqrt(1-(bot/cs)^2)*sqrt(1-bot^2/cp^2);
    ccc=tanh(2*pi*f*H/cc*sqrt(1-cc^2/c1^2))*(rhob1*cc^4*sqrt(1-cc^2/cp^2))/ ...
(rhob2*cs^4*sqrt(1-cc^2/c1^2))+(2-cc^2/cs^2)^2-4*sqrt(1-(cc/cs)^2)*sqrt(1-cc^2/cp^2);
   if aaa*bbb>0
       ff=0;    
       break
   else
       ff=1;
         if aaa>0&&ccc>0
               top=cc;cc=(bot+top)/2;
         else
             if aaa>0&&ccc<0
                  bot=cc;cc=(bot+top)/2;
             else
                 if aaa<0&&ccc>0
                       bot=cc;cc=(bot+top)/2;
                 else
                     if aaa<0&&ccc<0
                           top=cc;cc=(bot+top)/2;
                     end
                  end
             end
         end
   end
    end
    if ff==0
        continue
    else 
    format short;   %%输出长整型
    CC(xl)=cc;xl=xl+1;k=2*pi*f/cc;    
%     fprintf('The result:\n方程的根c=%.6f;\n水平波数k=%.6f;\n',cc,k);
    end
end  
% % 深海tanh近似；
% cc=zeros(1500,1);tanhcc=zeros(size(cc));H=2000;
% for c=1:1:1500
%     c1=1500;
%     cc(c)=c;
%     tanhcc(c)=tanh(2*pi*H*sqrt((c/c1)^2)/c);
% end
% plot(cc(1:length(cc)),tanhcc(1:length(cc)),'LineWidth',0.5);
% grid on;set(gca,'FontSize',14);
% xlabel('C(m/s)','fontsize',16);
% ylabel('tanh','fontsize',14);
% title(['H=' num2str(H) '深海tanh近似'],'fontsize',14);