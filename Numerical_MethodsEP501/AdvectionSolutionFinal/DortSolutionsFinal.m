%% Evan Dort EP501 Final Exam
clc
clear
close all
%% Problem 3 Part a and b


AdvectAllTheThings(100,9000)
AdvectAllTheThingsArty(100,338)
AdvectAllTheThingsImplicit(100,200)
AdvectAllTheThingsImplicitBig(100,50)
AdvectAllTheThingsImplicitTooBig(100,20)





%% Called Functions

function AdvectAllTheThings(xstep,tstep)

v = 20;
lambda = 0.25;
x = linspace(0,1,xstep);
t = linspace(0,0.05,tstep);
dx = x(2)-x(1);
dt = t(2)-t(1);
lengx = numel(x);
lengt = numel(t);
coeffi = v*(dt/dx);
fx0 = exp(-((x-0.5).^8)/(2*(0.1).^8));

advected = zeros(lengt,lengx);
advected(1,:) = fx0;



for i = 1:lengt-1
    for j = 2:lengx-1
        advected(i+1,j) = (advected(i,j)-coeffi*(advected(i,j)-advected(i,j-1))+((lambda*dt)/(dx^2))*(advected(i,j+1)-2*advected(i,j)+advected(i,j-1)));
    end %for
end %for

figure

imagesc(t,x,advected)
axis xy
shading flat
title('Explicit')
end

function AdvectAllTheThingsArty(xstep,tstep)

v = 20;
lambda = 0.25;
x = linspace(0,1,xstep);
t = linspace(0,0.05,tstep);
dx = x(2)-x(1);
dt = t(2)-t(1);
lengx = numel(x);
lengt = numel(t);
coeffi = v*(dt/dx);
fx0 = exp(-((x-0.5).^8)/(2*(0.1).^8));

advected = zeros(lengt,lengx);
advected(1,:) = fx0;



for i = 1:lengt-1
    for j = 2:lengx-1
        advected(i+1,j) = (advected(i,j)-coeffi*(advected(i,j)-advected(i,j-1))+((lambda*dt)/(dx^2))*(advected(i,j+1)-2*advected(i,j)+advected(i,j-1)));
    end %for
end %for

figure

imagesc(t,x,advected)
axis xy
shading flat
title('Explicit stability edge')
end

function AdvectAllTheThingsImplicit(xstep,tstep)

v = 20;
lambda = 0.25;
x = linspace(0,1,xstep);
t = linspace(0,0.05,tstep);
dx = x(2)-x(1);
dt = t(2)-t(1);
lengx = numel(x);
lengt = numel(t);
fx0 = exp(-((x-0.5).^8)/(2*(0.1).^8));

advected = zeros(lengt,lengx);
advected(1,:) = fx0;


k=2;
for i = 1:lengt-1
    for j = 2:lengx-1
        MyMatSol(j,1) = -v/dx*advected(i,j-1)+(v/dx-1/dt)*advected(i,j);
        MyMat(k,j+1) = lambda/dx^2;
        MyMat(k,j) = -(1/dt)-(2*lambda)/(dx^2);
        MyMat(k,j-1) = lambda/dx^2;
        k=k+1;
    end %for
    k=2;
    advected(i+1,:) = MyMat\MyMatSol;
end %for

figure('Name','Implicit')

imagesc(t,x,advected)
axis xy
shading flat
title('Implicit Method')
end

function AdvectAllTheThingsImplicitBig(xstep,tstep)

v = 20;
lambda = 0.25;
x = linspace(0,1,xstep);
t = linspace(0,0.05,tstep);
dx = x(2)-x(1);
dt = t(2)-t(1);
lengx = numel(x);
lengt = numel(t);
fx0 = exp(-((x-0.5).^8)/(2*(0.1).^8));

advected = zeros(lengt,lengx);
advected(1,:) = fx0;


k=2;
for i = 1:lengt-1
    for j = 2:lengx-1
        MyMatSol(j,1) = -v/dx*advected(i,j-1)+(v/dx-1/dt)*advected(i,j);
        MyMat(k,j+1) = lambda/dx^2;
        MyMat(k,j) = -(1/dt)-(2*lambda)/(dx^2);
        MyMat(k,j-1) = lambda/dx^2;
        k=k+1;
    end %for
    k=2;
    advected(i+1,:) = MyMat\MyMatSol;
end %for

figure('Name','Implicit Big')

imagesc(t,x,advected)
axis xy
shading flat
title('Yup even with a larger time step it is stable')
end

function AdvectAllTheThingsImplicitTooBig(xstep,tstep)

v = 20;
lambda = 0.25;
x = linspace(0,1,xstep);
t = linspace(0,0.05,tstep);
dx = x(2)-x(1);
dt = t(2)-t(1);
lengx = numel(x);
lengt = numel(t);
fx0 = exp(-((x-0.5).^8)/(2*(0.1).^8));

advected = zeros(lengt,lengx);
advected(1,:) = fx0;


k=2;
for i = 1:lengt-1
    for j = 2:lengx-1
        MyMatSol(j,1) = -v/dx*advected(i,j-1)+(v/dx-1/dt)*advected(i,j);
        MyMat(k,j+1) = lambda/dx^2;
        MyMat(k,j) = -(1/dt)-(2*lambda)/(dx^2);
        MyMat(k,j-1) = lambda/dx^2;
        k=k+1;
    end %for
    k=2;
    advected(i+1,:) = MyMat\MyMatSol;
end %for


figure('Name','Implicit Too Big')

imagesc(t,x,advected)
axis xy
shading flat
title('Yup when the time step is too large there are induced artifacts')
end



