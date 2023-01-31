%%%%%%%%AMPLITUDE and FREQUENCY of ripples in mm and rad/mm%%%%%
three_d_amp = .0005;
freq = 2*pi; %angular frequency
lambda = 60; %in nm 
%%%%%%%%%%%%%%%%%INITIAL LIGHT RAYS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
three_d_rays = [];
for i = -500:10:500
    for j = -sqrt(500^2-i^2):10:sqrt(500^2-i^2)
        three_d_rays = [three_d_rays; i j];
    end
end
radius = sqrt(three_d_rays(:,1).^2+three_d_rays(:,2).^2);
%%%%%%%%%%%%%%FIND NORMAL OF TANGENT PLANES%%%%%%%%%%%%%%%%%%%%%
n = [];
for i = 1:size(three_d_rays,1)
    n = [n;-1*(1/2000)*three_d_rays(i,1)+three_d_amp*sin(freq*radius(i))*freq*three_d_rays(i,1)/radius(i), -1*(1/2000)*three_d_rays(i,2)+three_d_amp*sin(freq*radius(i))*freq*three_d_rays(i,2)/radius(i), 1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIND RAY Z COMPONENT%%%%%%%%%%%%%%
z = [];
for i = 1:size(three_d_rays,1)
    z = [z; (n(i,1)^2+n(i,2)^2+1)/2];
end
%%%%%%%%%%%%%%%%%%FIND REF VECTOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subtracting normal and ray vector
ref_vec = [];
for i = 1:size(n,1)
    ref_vec = [ref_vec; n(i,1), n(i,2), n(i,3)-z(i)];
end
%%%%%%%%%%%%%%%%%%%%%%%FIND T WHEN Z=1000%%%%%%%%%%%%%%%%%%%%%%%
%3D LINE%
t = [];
for i = 1:size(ref_vec,1)
    t = [t; (1000-((1/4000)*three_d_rays(i,1)^2+(1/4000)*three_d_rays(i,2)^2+three_d_amp*cos(freq*radius(i))))/ref_vec(i,3)];
end
%%%%%%%%%%%%%%%%%%%%%%%%PLOT AT T%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coords = [];
for i = 1:length(t)
    coords = [coords;three_d_rays(i,1)+t(i)*ref_vec(i,1) three_d_rays(i,2)+t(i)*ref_vec(i,2)];
end
figure(4)
plot(coords(:,1),coords(:,2), '*')
xlabel('mm')
ylabel('mm')
title('Focal Plane Scatterplot')

figure(5)
hist3([coords(:,1),coords(:,2)],'CDataMode','auto','FaceColor','interp','nbins',[25 25]);
colorbar
xlabel('mm')
ylabel('mm')
title('Point Spread Function')

percentage = 0;
max_rad = 0;
while percentage < .9
    %YOU CAN CHANGE RADIUS INCREMENT SIZE HERE
    max_rad = max_rad + 1e-2;
    fit_coords = [];
    for i = 1:size(coords,1)
        if sqrt(coords(i,1)^2+coords(i,2)^2)< max_rad
            fit_coords = [fit_coords; coords(i,1),coords(i,2)];
        end
    end
    percentage = size(fit_coords,1)/size(coords,1);
end
percentage = percentage * 100;
fprintf("A radius of %f mm will absorb %f percent of the rays.\n", max_rad, percentage)

figure(6)
plot(fit_coords(:,1),fit_coords(:,2), '*')
xlabel('mm')
ylabel('mm')
title('90% Fit Focal Plane Scatter Plot')

%%%%%%%%%%%%%%%%%%%%%%%%%%RMS OPD and STREHL%%%%%%%%%%%%%%%%%%%%
surf = [];
for i = 1:size(three_d_rays,1)
    surf = [surf; three_d_rays(i,1) three_d_rays(i,2) (three_d_rays(i,1)^2+three_d_rays(i,2)^2)/4000 + three_d_amp*cos(radius(i))];
end

%Perfect Distance
p_dist = [];
for i = 1:size(surf,1)
    p_dist = [p_dist; sqrt((0-surf(i,1))^2+(0-surf(i,2))^2+(1000-surf(i,3))^2)];
end

%Actual Distance
a_dist = [];
for i = 1:size(surf,1)
    a_dist = [a_dist; sqrt((coords(i,1)-surf(i,1))^2+(coords(i,2)-surf(i,2))^2+(1000-surf(i,3))^2)];
end

%path displacement
dis = [];
for i = 1:length(a_dist)
    %convert distance to nm from mm
    dis = [dis; abs(p_dist(i)*1000000-a_dist(i)*1000000)];
end


%RMSE
sub = (dis).^2;
for i = 1:length(sub)
    if isnan(sub(i))
        sub(i) = 0;
    end
end

%in units of nm (average path diff)
rms = sqrt(sum(sub)/length(a_dist));

%average path displacement in wavelengths
rms_wavelengths = rms/lambda;
strehl = exp(-(2*pi*rms_wavelengths)^2);
fprintf("This mirror has a Strehl ratio of %f and an RMS OPD of %f wavelengths.\n", strehl, rms_wavelengths)
