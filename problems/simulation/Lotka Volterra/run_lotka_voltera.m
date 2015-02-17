% run lotka voltera
clear
par = [3 2 1 2 1 1];
for k=1:length(par)
    eval(['c' num2str(k) '= par(' num2str(k) ');']);
end

x_i=[10 15]; %hares
y_i=[5 10]; %lynxes
%IC = [6,3];
% xs=[1:10];
% ys=[1:10];
% wins = ones(length(xs),length(ys));

time_span = 5;
samp_time=.025;
%figure;
options = simset('SrcWorkspace','current');
% for i=1:length(xs)
%     for j=1:length(ys)
%         x_i=xs(i);
%         y_i=ys(j);
%         try
%outstr='c1*u(1)-c2*u(1)*u(2)-c3*u(1)^2';
allx1=[];
allx2=[];
alldx1=[];
alldx2=[];
for i=1:length(x_i)
    for j=1:length(y_i)
        IC = [x_i(i),y_i(j)];
        sim('lotka_voltera.mdl');
        allx1=[allx1;x1];
        allx2=[allx2;x2];
        alldx1=[alldx1;dx1];
        alldx2=[alldx2;dx2];
    end
end

lv1 = [alldx1,allx1,allx2];
lv2 = [alldx2,allx1,allx2];
        figure;
        plot(allx1,'b'); hold on;
        plot(allx2,'r');
        plot(alldx1,'+b'); hold on;
        plot(alldx2,'+r'); hold on;
        legend('x_1','x_2','dx_1','dx_2');
        xlabel('Time step');
        ylabel('Population');
%         catch
%             wins(i,j)=0;
%         end
%     end
% end
% figure;
% surf(xs,ys,wins);
% xlabel('x_i');
% ylabel('y_i');
