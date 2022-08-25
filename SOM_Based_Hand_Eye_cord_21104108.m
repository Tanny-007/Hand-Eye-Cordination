clc;clear all; close all;

% Hello sir, once I got perfect Circular trajectory but by mistake I saved 
% it as a .png format. Later on I tried to get that trajectory but I didnt get 
% that kind of trajectory. And I guess, to be able to get that Perfect trajectory 
% we need so much iterations which requires much time as I had seen on the
% result section of your paper(Table no. 1). So I am attaching the results
% along with that png file I got. 

%% 3D lattice size =(12x7x4)
for i = 1:12
    for j = 1:7
        for k = 1:4
            for l = 1:3
                for m = 1:3
                    a(i,j,k,l,m) = 0.1*rand;  
                end
            end
        end
    end
end

for i = 1:12
    for j = 1:7
        for k = 1:4
            for l = 1:3
                theta(i,j,k,l) = 0.1*rand;
            end
        end
    end
end

for i = 1:12
    for j = 1:7
        for k = 1:4
            for l = 1:3
                w(i,j,k,l) = 0.1*rand;
            end
        end
    end
end
l1 = 0.254;
l2 = 0.254;
l3 = 0.254;
%eta initialization
etai = 1;
etaf = 0.05;
%epsilon initialization
epsi = 0.9;
epsf = 0.9;
%sigma initialization
sigmai = 2;
sigmaf = 0.1;
%sigma1 initialization
sigma1i = 2.5;
sigma1f = 0.01;
%initial values
t0 = 0.05;
t=0:(2*pi);
v0 = ones(3,1);
vf = ones(3,1);
w_gamma = zeros(3,1);
a_gamma = zeros(3,3);
theta_gamma = zeros(3,1);
del_th = ones(3,1);
del_aws = ones(3,1);
del_a = ones(3,3);

tmax = 50000;
for itrain = 1:tmax
    
    xd = 3*rand-1.5;
    yd = rand;
    zd = rand;
    x = [xd;yd;zd];
    
    for i = 1:12
        for j = 1:7
            for k = 1:4
                d(i,j,k) = sqrt((w(i,j,k,1)-xd)^2+(w(i,j,k,2)-yd)^2+(w(i,j,k,3)-zd)^2);
            end
        end
    end
    m1 = 10000;
    for i = 1:12
        for j = 1:7
            for k = 1:4
                m2 = min(m1,d(i,j,k));
                  if(m2<m1)
                    iw = i;
                    jw = j;
                    kw = k;
                  end
                m1 = m2;
            end
        end
    end
    eta = etai*(etaf/etai)^(itrain/tmax);
    eps = epsi*(epsf/epsi)^(itrain/tmax);
    sigma = sigmai*(sigmaf/sigmai)^(itrain/tmax);
    sigma1 = sigma1i*(sigma1f/sigma1i)^(itrain/tmax);
    h_sum = 0; %Sum Of h gamma function
    th_0_out = zeros(3,1);
    th_1_out = zeros(3,1);
    
    %% Course action
    for i = 1:12
        for j = 1:7
            for k = 1:4
                r = [i;j;k];
                s = [iw;jw;kw];  %LINE 84
                h_gamma = exp(-(norm(r-s))^2/(2*sigma1^2));
                theta_gamma = [theta(i,j,k,1);theta(i,j,k,2);theta(i,j,k,3)];
                a_gamma = [a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                w_gamma = [w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                th_0_out = th_0_out+h_gamma*(theta_gamma+a_gamma*(x-w_gamma)); %{thg, wg, x = 3*1} and {ag = 3*3} 
                h_sum = h_sum + h_gamma;
            end
        end
    end
    th_0_out = th_0_out/h_sum;
        
        R = l2*cos(th_0_out(2))+l3*cos(th_0_out(3))+t0;
        v0(3) = l2*sin(th_0_out(2))+l3*sin(th_0_out(3))+l1;
        v0(1) = R*cos(th_0_out(1));
        v0(2) = R*sin(th_0_out(1));
    delth = zeros(3,1);
    
     %% Fine Action 
        for i = 1:12
            for j = 1:7
                for k = 1:4
                    r = [i;j;k];
                    s = [iw;jw;kw];
                    h_gamma = exp(-(norm(r-s))^2/(2*sigma1^2));
                    a_gamma = [a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                    w_gamma = [w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                    delth = delth + h_gamma*a_gamma*(x-v0);
                end
            end
        end
        th_1_out = th_0_out + delth/h_sum; 
        
        R = l2*cos(th_1_out(2))+l3*cos(th_1_out(3))+t0;
            vf(3) = l2*sin(th_1_out(2))+l3*sin(th_1_out(3))+l1; %z 
            vf(1) = R*cos(th_1_out(1));                         %x 
            vf(2) = R*sin(th_1_out(1));                         %y 
        
      %% weight vector(w)update
        for i = 1:12           
            for j = 1:7
                for k = 1:4
                    r = [i;j;k];
                    s = [iw;jw;kw];
                    h_gamma = exp(-(norm(r-s))^2/(2*sigma1^2));
                    for l = 1:3
                        w(i,j,k,l) = w(i,j,k,l) + eta*h_gamma*(x(l)-w(i,j,k,l));
                    end
                end
            end
        end
      %% theta vector(theta) and Jacobian Matrix(a) update
        
        del_th = zeros(3,1);
        del_aws = zeros(3,1);
        for i = 1:12
            for j = 1:7
                for k = 1:4
                    r = [i;j;k];
                    s = [iw;jw;kw];
                    h_gamma = exp(-(norm(r-s))^2/(2*sigma1^2));
                    theta_gamma = [theta(i,j,k,1);theta(i,j,k,2);theta(i,j,k,3)];
                    a_gamma = [a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                    w_gamma = [w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                    
                    del_th = del_th + h_gamma*(theta_gamma+a_gamma*(v0-w_gamma)); %theta0 out
                    del_aws = del_aws + h_gamma*a_gamma*(vf-v0);     %theta1 out
                end
            end
        end
        del_t = zeros(3,1);
        for i1 = 1:12
            for j1 = 1:7
                for k1 = 1:4
                    r1 = [i1;j1;k1];
                    s = [iw;jw;kw];
                    h_gamma_1 = exp(-(norm(r1-s))^2/(2*sigma1^2));
                    del_t = (th_0_out-(del_th/h_sum))*(h_gamma_1/h_sum); 
                    del_a = ((th_1_out-th_0_out)-(del_aws/h_sum))*(vf-v0)'*(h_gamma_1/(h_sum*norm(vf-v0)^2));
                    for  l = 1:3
                        theta(i1,j1,k1,l) = theta(i1,j1,k1,l)+eps*del_t(l); %theta update
                        for m = 1:3
                            a(i1,j1,k1,l,m) = a(i1,j1,k1,l,m) + eps*del_a(l,m); %ag update
                        end
                    end
                end
            end
        end
        error(itrain)=0.5*((norm(x-vf))^2);
end
figure(1)
plot(error);

%% Testing to track a circle
c=1;
xm=0;
ym=0;
zm=0;
for t=0:0.15:2*pi
    xd1 = 0.1+0.08*cos(t);
    yd1 = 0.2+0.08*sin(t);
    zd1 = 0.12;
    x=[xd1;yd1;zd1];
    for i=1:12
        for j=1:7
            for k=1:4
                d(i,j,k)=sqrt((w(i,j,k,1)-xd1)^2+(w(i,j,k,2)-yd1)^2+(w(i,j,k,3)-zd1)^2);
            end
        end
    end
    m1=50000;
    for i=1:12
        for j=1:7
            for k=1:4
                m2=min(m1,d(i,j,k));
                if m2<m1
                    iw=i;
                    jw=j;
                    kw=k;
                end
                m1=m2;
            end
        end
    end
    eta=etai*(etaf/etai)^(itrain/tmax);
    eps=epsi*(epsf/epsi)^(itrain/tmax);
    sigma=sigmai*(sigmaf/sigmai)^(itrain/tmax);
    sigma1=sigma1i*(sigma1f/sigma1i)^(itrain/tmax);
    h_sum=0;
    th_0_out=zeros(3,1);
    th_1_out=zeros(3,1);
    for i=1:12
        for j=1:7
            for k=1:4
                r=[i;j;k];
                s=[iw;jw;kw];
                h_gamma=exp(-(norm(r-s))^2/(2*sigma1^2));
                th_gamma=[theta(i,j,k,1);theta(i,j,k,2);theta(i,j,k,3)];
                a_gamma=[a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                w_gamma=[w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                th_0_out=th_0_out+h_gamma*(th_gamma+a_gamma*(x-w_gamma));
                h_sum=h_sum+h_gamma;
            end
        end
    end
    th_0_out=th_0_out/h_sum;
    R=l2*cos(th_0_out(2))+l3*cos(th_0_out(3))+t0;
    v0(1)=R*cos(th_0_out(1));
    v0(2)=R*sin(th_0_out(1));
    v0(3)=l2*sin(th_0_out(2))+l3*sin(th_0_out(3))+l1;
    for i=1:12
        for j=1:7
            for k=1:4
                r=[i;j;k];
                s=[iw;jw;kw];
                h_gamma=exp(-(norm(r-s))^2/(2*sigma^2));
                a_gamma=[a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                w_gamma=[w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                delth=delth+h_gamma*a_gamma*(x-v0);
            end
        end
    end
    th_1_out=th_0_out+delth/h_sum;
    R=l2*cos(th_1_out(2))+l3*cos(th_1_out(3))+t0;
    vf(1)=R*cos(th_1_out(1));
    vf(2)=R*sin(th_1_out(1));
    vf(3)=l2*sin(th_1_out(2))+l3*sin(th_1_out(3))+l1;
    xm(1,c)=vf(1);
    ym(1,c)=vf(2);
    zm(1,c)=vf(3);
    c=c+1;
end
figure(2)
plot3(xm,ym,zm);

%% Testing to track a Line 
c=1;
xm1=0;
ym2=0;
zm3=0;
for t=1:2
    xd1 = 0.1+0.2*t;
    yd1 = 0.2+0.3*t;
    zd1 = 0.4;
    x=[xd1;yd1;zd1];
    for i=1:12
        for j=1:7
            for k=1:4
                d(i,j,k)=sqrt((w(i,j,k,1)-xd1)^2+(w(i,j,k,2)-yd1)^2+(w(i,j,k,3)-zd1)^2);
            end
        end
    end
    m1=50000;
    for i=1:12
        for j=1:7
            for k=1:4
                m2=min(m1,d(i,j,k));
                if m2<m1
                    iw=i;
                    jw=j;
                    kw=k;
                end
                m1=m2;
            end
        end
    end
    eta=etai*(etaf/etai)^(itrain/tmax);
    eps=epsi*(epsf/epsi)^(itrain/tmax);
    sigma=sigmai*(sigmaf/sigmai)^(itrain/tmax);
    sigma1=sigma1i*(sigma1f/sigma1i)^(itrain/tmax);
    h_sum=0;
    th_0_out=zeros(3,1);
    th_1_out=zeros(3,1);
    for i=1:12
        for j=1:7
            for k=1:4
                r=[i;j;k];
                s=[iw;jw;kw];
                h_gamma=exp(-(norm(r-s))^2/(2*sigma1^2));
                th_gamma=[theta(i,j,k,1);theta(i,j,k,2);theta(i,j,k,3)];
                a_gamma=[a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                w_gamma=[w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                th_0_out=th_0_out+h_gamma*(th_gamma+a_gamma*(x-w_gamma));
                h_sum=h_sum+h_gamma;
            end
        end
    end
    th_0_out=th_0_out/h_sum;
    R=l2*cos(th_0_out(2))+l3*cos(th_0_out(3))+t0;
    v0(1)=R*cos(th_0_out(1));
    v0(2)=R*sin(th_0_out(1));
    v0(3)=l2*sin(th_0_out(2))+l3*sin(th_0_out(3))+l1;
    for i=1:12
        for j=1:7
            for k=1:4
                r=[i;j;k];
                s=[iw;jw;kw];
                h_gamma=exp(-(norm(r-s))^2/(2*sigma^2));
                a_gamma=[a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                w_gamma=[w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                delth=delth+h_gamma*a_gamma*(x-v0);
            end
        end
    end
    th_1_out=th_0_out+delth/h_sum;
    R=l2*cos(th_1_out(2))+l3*cos(th_1_out(3))+t0;
    vf(1)=R*cos(th_1_out(1));
    vf(2)=R*sin(th_1_out(1));
    vf(3)=l2*sin(th_1_out(2))+l3*sin(th_1_out(3))+l1;
    xm1(1,c)=vf(1);
    ym1(1,c)=vf(2);
    zm1(1,c)=vf(3);
    c=c+1;
end
figure(3)
plot3(xm1,ym1,zm1);

%% Testing to track Five Random points 
c=1;
xm=0;
ym=0;
zm=0;
for t=1:5
    xd1 = rand;
    yd1 = rand;
    zd1 = rand;
    x=[xd1;yd1;zd1];
    for i=1:12
        for j=1:7
            for k=1:4
                d(i,j,k)=sqrt((w(i,j,k,1)-xd1)^2+(w(i,j,k,2)-yd1)^2+(w(i,j,k,3)-zd1)^2);
            end
        end
    end
    m1=50000;
    for i=1:12
        for j=1:7
            for k=1:4
                m2=min(m1,d(i,j,k));
                if m2<m1
                    iw=i;
                    jw=j;
                    kw=k;
                end
                m1=m2;
            end
        end
    end
    eta=etai*(etaf/etai)^(itrain/tmax);
    eps=epsi*(epsf/epsi)^(itrain/tmax);
    sigma=sigmai*(sigmaf/sigmai)^(itrain/tmax);
    sigma1=sigma1i*(sigma1f/sigma1i)^(itrain/tmax);
    h_sum=0;
    th_0_out=zeros(3,1);
    th_1_out=zeros(3,1);
    for i=1:12
        for j=1:7
            for k=1:4
                r=[i;j;k];
                s=[iw;jw;kw];
                h_gamma=exp(-(norm(r-s))^2/(2*sigma1^2));
                th_gamma=[theta(i,j,k,1);theta(i,j,k,2);theta(i,j,k,3)];
                a_gamma=[a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                w_gamma=[w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                th_0_out=th_0_out+h_gamma*(th_gamma+a_gamma*(x-w_gamma));
                h_sum=h_sum+h_gamma;
            end
        end
    end
    th_0_out=th_0_out/h_sum;
    R=l2*cos(th_0_out(2))+l3*cos(th_0_out(3))+t0;
    v0(1)=R*cos(th_0_out(1));
    v0(2)=R*sin(th_0_out(1));
    v0(3)=l2*sin(th_0_out(2))+l3*sin(th_0_out(3))+l1;
    for i=1:12
        for j=1:7
            for k=1:4
                r=[i;j;k];
                s=[iw;jw;kw];
                h_gamma=exp(-(norm(r-s))^2/(2*sigma^2));
                a_gamma=[a(i,j,k,1,1),a(i,j,k,1,2),a(i,j,k,1,3);a(i,j,k,2,1),a(i,j,k,2,2),a(i,j,k,2,3);a(i,j,k,3,1),a(i,j,k,3,2),a(i,j,k,3,3)];
                w_gamma=[w(i,j,k,1);w(i,j,k,2);w(i,j,k,3)];
                delth=delth+h_gamma*a_gamma*(x-v0);
            end
        end
    end
    th_1_out=th_0_out+delth/h_sum;
    R=l2*cos(th_1_out(2))+l3*cos(th_1_out(3))+t0;
    vf(1)=R*cos(th_1_out(1));
    vf(2)=R*sin(th_1_out(1));
    vf(3)=l2*sin(th_1_out(2))+l3*sin(th_1_out(3))+l1;
    xm(1,c)=vf(1);
    ym(1,c)=vf(2);
    zm(1,c)=vf(3);
    c=c+1;
end
figure(4)
plot3(xm,ym,zm);















