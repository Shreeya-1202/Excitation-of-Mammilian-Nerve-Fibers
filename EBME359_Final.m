clear
close all
clc

%% Question 

sigma = 2e-6; % S/um
Nodes = 41;
n = 1:Nodes;
x = -200; % um
y = 500; % um
current = 2; % mA
del2V = zeros(3,Nodes-2);

D = [2, 12, 20]; % um
for i = 1:length(D)
    L = 100.*D(i);
    axons = (-21+n) .* L;
    r = sqrt((x-axons).^2 + y^2);
    Vext = current./(4*pi*sigma.*r);
    for j = 2:Nodes-1
        del2V(i,j-1) = Vext(j-1) - 2*Vext(j) + Vext(j+1);
    end
end
figure
plot(n(2:end-1),del2V)
xlabel('Node')
ylabel('2nd Difference (mV)')
legend(strcat('D = ',string(D),' um'))

% The max. 2nd difference is at the center node (21) since as D increases,
% L increases and the electrode is closest to node 21. Only for D = 2 um is
% the electrode directly over node 20 instead, which is also reflected in
% the plot.

%% Question 3

x = -6000;

D = [2, 12, 20]; % um
for i = 1:length(D)
    L = 2.501 ./(0.001593549 + 0.03234*exp(-0.42*D(i)))+0.8
    axons = (-21+n) .* L;
    r = sqrt((x-axons).^2 + y^2);
    Vext = current./(4*pi*sigma.*r);
    for j = 2:Nodes-1
        del2V(i,j-1) = Vext(j-1) - 2*Vext(j) + Vext(j+1);
    end
end
figure
plot(n(2:end-1),del2V)
xlabel('Node')
ylabel('2nd Difference (mV)')
legend(strcat('D = ',string(D),' um'))

% The diameter changes the node at which the stimulation occurs for a
% larger magnitude x. At D = 2 um, L = 200 um, the electrode exists at
% 1000/200 = 5 nodes from the center node, therefore the stimulation occurs
% at node 16. At D = 12 um, the stim. occurs at node 20, and at D = 20 um
% the stim. occurs at node 21 by similar logic (closest node to electrode).

%% Question 4

D = 12;

x = [200, 6000, 10000];
for i = 1:length(x)
    L = 100*D;
    axons = (-21+n) .* L;
    r1 = sqrt((-x(i)-axons).^2 + y^2);
    r2 = sqrt((x(i)-axons).^2 + y^2);
    Vext = (current/(4*pi*sigma)).*(1./r1 - 1./r2);
    for j = 2:Nodes-1
        del2V(i,j-1) = Vext(j-1) - 2*Vext(j) + Vext(j+1);
    end
end
figure
plot(n(2:end-1),del2V)
xlabel('Node')
ylabel('2nd Difference (mV)')
legend(strcat('x = ',string(x),' um'))

% b. The maximum occurs at nodes 16 and 26. This is symmetric about the
% center node (21) but not at the center node, since the electrodes are
% located above nodes 16 and 26 and not 21. As the distance between
% electrodes increases, the distance between the stimulations also
% increases.