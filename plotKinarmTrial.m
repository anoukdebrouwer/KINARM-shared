function plotKinarmTrial(fig1,rep,tp,stim,limXY,xyHand,xyGaze,vxyHand,vGaze_angular)
% plotKinarmTrial Plot hand and gaze data of single Kinarm trial
%
% plotKinarmTrial(fig1,rep,tp,stim,limXY,xyHand,xyGaze,vxyHand,vGaze_angular)
% Used by processKinarmData.m

% MIT License
% Copyright (c) 2021 Anouk de Brouwer

% get stimuli
startXY = [stim.startX stim.startY];
startR = [stim.startRadius stim.startRadius];
targetXY = [stim.targetX stim.targetY];
if strcmp(stim.targetShape,'circle')
    targetR = [stim.targetRadius stim.targetRadius];
    curv = [1 1];
elseif strcmp(stim.targetShape,'rect')
    targetR = [stim.targetWidth/2 stim.targetHeight/2];
    curv = [0 0];
end

figure(fig1); clf

% SUBPLOT 1
subplot(2,3,[1 4]); hold on
% plot start position and target
rectangle('Position',[startXY-startR startR*2],'Curvature',[1 1]);
rectangle('Position',[targetXY-targetR targetR*2],'Curvature',curv);
% plot hand and gaze trajectory
pg = plot(xyGaze(:,1),xyGaze(:,2),'.-','LineWidth',1,'MarkerSize',15);
ph = plot(xyHand(:,1),xyHand(:,2),'-','LineWidth',2);
% axes and labels
hold off
axis equal
legend([ph pg],{'hand','gaze'},'location','east');
xlabel('X position (m)')
ylabel('Y position (m)')
title(['Type ' num2str(tp) ' - Rep ' num2str(rep)]);

% SUBPLOT 2: plot hand and gaze position over time
subplot(2,3,[2 5]); hold on
% plot hand and gaze trajectory
plot(xyHand(:,1),':','Color',ph.Color,'LineWidth',2);
plot(xyHand(:,2),'--','Color',ph.Color,'LineWidth',2);
plot(xyGaze(:,1),':','Color',pg.Color,'LineWidth',2);
plot(xyGaze(:,2),'--','Color',pg.Color,'LineWidth',2)
% axes, legend and labels
hold off
xlim([0 length(xyHand)])
legend('x hand','y hand','x gaze','y gaze','location','east');
xlabel('Time since target onset (ms)')
ylabel('Position (m)')
title('Hand and gaze position')

% SUBPLOT 3: plot hand velocity over time
subplot(2,3,3); hold on
% plot hand velocity
plot(vxyHand(:,1),':','Color',ph.Color,'LineWidth',1)
plot(vxyHand(:,2),'--','Color',ph.Color,'LineWidth',1)
vHand = sqrt(vxyHand(:,1).^2+vxyHand(:,2).^2);
plot(vHand,'-','Color',ph.Color,'LineWidth',2);
% axes, legend and labels
hold off
axis([0 length(xyHand) -0.5 1.5])
legend('vx hand','vy hand','v hand','location','northeast');
xlabel('Time since target onset (ms)')
ylabel('Velocity (m/s)')
title('Hand velocity')

% SUBPLOT 4: plot angular gaze velocity over time
subplot(2,3,6); hold on
% plot gaze velocity
plot(vGaze_angular,'-','Color',pg.Color,'LineWidth',2);
% axes and labels
hold off
axis([0 length(vGaze_angular) 0 1500])
xlabel('Time since target onset (ms)')
ylabel('Angular velocity (deg/s)')
title('Gaze velocity')

end