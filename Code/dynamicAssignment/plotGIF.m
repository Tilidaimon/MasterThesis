% ---------------mainGIF.m---------
% ------------����gif---------
Delay = 2000 % �ӳٲ���
set(0,'defaultfigurecolor','w');
% ------------��һ��GIF-------û�й켣��---------
pic_num = 1;
steps = length(targets_save.p)
a = 1;
for k = 1:100:steps+Delay
    if k <= steps
        figure(1);
        %     scrsz = get(0,'ScreenSize'); % ��Ļ��С
        %     set(gcf,'Position',scrsz);  %
        %******�Լ���plot
        % �ú��ͷ��ʾ�ػ�

        % ����ͷ��ʾ����
        x = missiles_save.p(k,:,1);
        y = missiles_save.p(k,:,2);
        %         plot(x,y,'bo')
        
        quiver(x',y', cos(missiles_save.angle(k,:))',sin(missiles_save.angle(k,:))',0.3,...
            'LineWidth',1,'AutoScale','off','color','b');
        hold on

        % ���̼�ͷ��ʾĿ��
        x = targets_save.p(k,:,1);
        y = targets_save.p(k,:,2);
        quiver(x',y', cos(targets_save.angle(1,:))',sin(targets_save.angle(1,:))',0.5,...
            'LineWidth',1.5,'AutoScale','off','color','g');
        axis([0 14*10 0 5*10]);
        
        % ����legend
       
        x = missiles_save.p(1,1,1);
        y = missiles_save.p(1,1,2);
        p1 = quiver(x,y, cos(missiles_save.angle(1,1))',sin(missiles_save.angle(1,1))',0.3,...
            'LineWidth',1,'AutoScale','off','color','b');
        hold on;
        x = targets_save.p(1,1,1);
        y = targets_save.p(1,1,2);
        p2 = quiver(x,y, cos(targets_save.angle(1,1))',sin(targets_save.angle(1,1))',0.5,...
            'LineWidth',1.5,'AutoScale','off','color','g');
        hold on;
        legend([p1,p2],["����","Ŀ��"],'Location','northwest');
        xlabel("x/km");ylabel("y/km");
        str = ['t=',num2str(k*0.02-0.02),'s'];
        text(5.5*10,4.5*10,str);
        hold off;
        
        %*************
        drawnow;
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        if pic_num == 1
            imwrite(I,map,'./GIF/move.gif','gif', 'Loopcount',inf,'DelayTime',0);
        else
            imwrite(I,map,'./GIF/move.gif','gif','WriteMode','append','DelayTime',0);
        end
        pic_num = pic_num + 1;
    end
    %*************
    drawnow;
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'./GIF/move.gif','gif', 'Loopcount',inf,'DelayTime',0);
    else
        imwrite(I,map,'./GIF/move.gif','gif','WriteMode','append','DelayTime',0);
    end
    pic_num = pic_num + 1;
    
end

% ------------�ڶ���GIF-------�й켣��---------
pic_num = 1;
steps = length(targets_save.p)
a = 1;
for k = 1:100:steps+Delay
    if k <= steps
        figure(1);
        %     scrsz = get(0,'ScreenSize'); % ��Ļ��С
        %     set(gcf,'Position',scrsz);  %
        %******�Լ���plot
        % ����ͷ��ʾ����
        x = missiles_save.p(k,:,1);
        y = missiles_save.p(k,:,2);
        %         plot(x,y,'bo')
        
        %             plot(x,y,'b.');
        quiver(x',y', cos(missiles_save.angle(k,:))',sin(missiles_save.angle(k,:))',0.3,...
            'LineWidth',0.8,'AutoScale','off','color','b');
        hold on
        % ���̼�ͷ��ʾĿ��
        x = targets_save.p(k,:,1);
        y = targets_save.p(k,:,2);
        plot(x,y,'g.');
        axis([0 14*10 0 5*10]);
        
        % ����legend
        
        x = missiles_save.p(1,1,1);
        y = missiles_save.p(1,1,2);
        p1 = quiver(x,y, cos(missiles_save.angle(1,1))',sin(missiles_save.angle(1,1))',0.3,...
            'LineWidth',0.8,'AutoScale','off','color','b');
        hold on;
        x = targets_save.p(1,1,1);
        y = targets_save.p(1,1,2);
        p2 = plot(x,y,'g.');
        hold on;
        legend([p1,p2],["����","Ŀ��"],'Location','northwest');
        xlabel("x/km");ylabel("y/km");
%         str = ['t=',num2str(k*0.02-0.02),'s'];
%         text(5.5*10000,2.6*10000,str);
        
        hold on;
        
        %*************
        drawnow;
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,256);
        if pic_num == 1
            imwrite(I,map,'./GIF/route.gif','gif', 'Loopcount',inf,'DelayTime',0);
        else
            imwrite(I,map,'./GIF/route.gif','gif','WriteMode','append','DelayTime',0);
        end
        pic_num = pic_num + 1;
    end
    %*************
    drawnow;
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
        imwrite(I,map,'./GIF/route.gif','gif', 'Loopcount',inf,'DelayTime',0);
    else
        imwrite(I,map,'./GIF/route.gif','gif','WriteMode','append','DelayTime',0);
    end
    pic_num = pic_num + 1;
    
end