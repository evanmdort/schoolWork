%% Evan Dort EP501 HW 2 Dr Zettergren
clc
clear

%% Problem 1
%% Part a
funct = @(d) d.^3-9*d;
eps = 0.1;
x0=2;
z=newton_approx(funct,eps,x0);
zed=round(z)
disp('Derivative for the function x^3-9x at point 2')
disp(zed)
% for some reason I am unable to get my answer to publish here. The code
% works but even leaving my answer variable 'zed' unsupressed, MATLAB
% refuses to print the answer when I publish the document it works for me
% though


%% part b
%I placed the plotting for the Bessel function at the end
%Since the plotting is meant to just estimate the root, I took the code for
%plotting the Bessel function directly from the MATLAB help page for
%besselj and modified it slightly
funct = @(d) besselj(0,d);
eps = 0.1;
x0=1;
bess=newton_approx(funct,eps,x0);
disp('First root of zero order Bessel function')
disp(bess)
%% part c
 funky= @(d) besselj(0,d);
 eps = 0.1;
 x0=1;
 bessie=zeros(6,1);
 bessieans=zeros(6,1);
 bessieans(1)=newton_approx(funky,eps,x0);
 x0=x0+bessieans(1);
 sols=2;
 for i=2:7
    bessie(i)= newton_approx(funky,eps,x0);
    difftol=bessie(i)-bessieans(i-1);
    if difftol>10^-6
        bessieans(sols)=bessie(i);
        sols=sols+1;
    end
    x0=x0+bessie(i);
 end
 
 for i=1:6
    Sid=sprintf('%d root of zero order Bessel Function',i);
    disp(Sid);
    disp(bessieans(i));
 end

%% Problem 2

%% Part a
power=5;
counter=1;
esoteric = @(d) d.^5-15*d.^4+85*d^3-225*d^2+274*d-120;
esotericp = @(d) 5*d^4-60*d^3+255*d^2-450*d+274;
rooties=0;
rootsans=zeros(power,1);
genesis=-10;
sols=2;
rootsans(1) = newton_exact(esoteric,esotericp,genesis);
genesis=-9;
while counter<20

rooties = newton_exact(esoteric,esotericp,genesis);
difftol=rootsans(sols-1)-rooties;
difftol=abs(difftol);
if abs(difftol)>10^-4
    rootsans(sols,1)=rooties;
    sols=sols+1;
end
counter=counter+1;
genesis=genesis+1;
end
disp('Roots of the function x^5-15x^4+85x^3-225x^2+274x-120')
disp(rootsans)

%% Part b
power=3;
counter=1;
imag=sqrt(-1);
entelechy = @(d) d^3-3*d^2+4*d-2;
entelechyp = @(d) 3*d^2-6*d+4;
genesis=-10;
sols=2;
elucidate=zeros(power,1);
elucidate(1) = newton_exact(entelechy,entelechyp,genesis);
genesis=-9;
while counter<20

rooties = newton_exact(entelechy,entelechyp,genesis);
difftol=elucidate(sols-1)-rooties;
difftol=abs(difftol);
if abs(difftol)>10^-4
    elucidate(sols,1)=rooties;
    sols=sols+1;
end
counter=counter+1;
genesis=genesis+1;
end
genesis=-10;
counter=1;
if sols<4
    while counter<20
        guess=genesis+genesis*imag;
        rooties = newton_exact(entelechy,entelechyp,guess);
        difftol=elucidate(sols-1)-rooties;
        difftol2=elucidate(1)-rooties;
        difftol=abs(difftol);
            if abs(difftol)>10^-4 && genesis ~= 0 && abs(difftol2)>10^-4
                 elucidate(sols,1)=rooties;
                 sols=sols+1;
            end
counter=counter+1;
genesis=genesis+1;
    end
end

disp('Real and complex roots of the function x^3-3x^2+4x-2');
disp(elucidate);


%% Problem 3
%% Part a
coefs=[2;-6;4];%define coeffs
a=coefs(1);
b=coefs(2);
c=coefs(3);
numep=-b+sqrt(b^2-4*a*c);%use quadratic formula, I separated the numerator and denominator for simplicity
numen=-b-sqrt(b^2-4*a*c);
denom=2*a;
root1=numep/denom;
root2=numen/denom;
disp('Roots of the function 2x^2-6x+4 found with quadratic formula');
disp(root1);
disp(root2);
%% Part b
pol=[1;-15;85;-225;274;-120];
root=5;
[m,n]=size(pol);
Qn=zeros(m,1);
Qn(1)=pol(1);
for i=2:m
   Qn(i)=(Qn(i-1)*root)+pol(i);
end
disp('Polynomial division of x^5-15x^4+85x^3-225x^2+274x-120 with root 5');
disp('The resultant polynomial is the following');
Sips=sprintf('x^4+(%d)x^3+(%d)x^2+(%d)x+(%d)',Qn(2),Qn(3),Qn(4),Qn(5));
disp(Sips);

%% Part c
    a=1;%define coeffs
    b=-15;
    c=85;
    e=-225;
    g=274;
    h=-120;
funct = @(d) a*d^5+b*d^4+c*d^3+e*d^2+g*d+h;
eps = 0.01;
x0=2.1;

pol=[a; b; c; e; g; h];
    [m,n]=size(pol);
    roots=zeros(m,1);
[m,n]=size(pol);
Cn=zeros(m,1);
Cn(1)=pol(1);
roots(1,1)=newton_approx(funct,eps,x0);
for i=2:m
   Cn(i)=(Cn(i-1)*roots(1,1))+pol(i);
end
Cn(abs(Cn)<1e-5)=0;%
Cink=round(Cn);
rouge=round(roots(1));
Yee=sprintf('Polynomial division of x^5-15x^4+85x^3-225x^2+274x-120 with root %d found via Newton Method',rouge);
disp(Yee);
disp('The resultant polynomial is the following');
Sips=sprintf('x^4+(%.1d)x^3+(%.1d)x^2+(%.1d)x+(%.1d)',Cink(2),Cink(3),Cink(4),Cink(5));
disp(Sips);  
%% Part d and e
times=1;
    a=1;%define coeffs
    b=-15;
    c=85;
    e=-225;
    g=274;
    h=-120;
    pol=[a; b; c; e; g; h];
    [m,n]=size(pol);
    rootsit=zeros(m-1,1);
   
while times<4 
    funct = @(d) a*d^5+b*d^4+c*d^3+e*d^2+g*d+h;
    eps = 0.01;
    x0=4;
    rootsit(times,1)=newton_approx(funct,eps,x0);
    pol=[a; b; c; e; g; h];
    Polys=zeros(m,1);
    Polys(1)=pol(1);
    for i=2:m
        Polys(i)=(Polys(i-1)*rootsit(times,1))+pol(i);
    end
    Polys(abs(Polys)<1e-5)=0;%
    a=Polys(1,1);%define coeffs
    b=Polys(2,1);
    c=Polys(3,1);
    e=Polys(4,1);
    g=Polys(5,1);
    h=Polys(6,1);
    times=times+1;
end
numep=-b+sqrt(b^2-4*a*c);%use quadratic formula, I separated the numerator and denominator for simplicity
numen=-b-sqrt(b^2-4*a*c);
denom=2*a;
rootsit(4,1)=numep/denom;
rootsit(5,1)=numen/denom;
disp('All roots found using Newton Method combined with polynomial division')
disp(rootsit)

%% Problem 4
%% Part a
x0=5;
y0=5;

f = @(x,y) x.^2+y.^2-2.*x-y;

gradf = @(x,y) grad_objfun2Df(x,y);

g = @(x,y) 0.25.*x.^2+y.^2-1;

gradg = @(x,y) grad_objfun2Dg(x,y);

[abstruse1,abstruse2] = newton2D_exact(f,gradf,g,gradg,x0,y0);
abstruse1(abs(abstruse1)<1e-5)=0;
abstruse2(abs(abstruse2)<1e-5)=0;
disp('The roots for the system in problem 4 part a are the following:')
disp(abstruse1)
disp(abstruse2)

%% Part b

x0=imag;
y0=1;
z0=2;

flee = @(x,y,z) x.^2+y.^2+z.^2-6;
gradflee = @(x,y,z) grad_objfun2Dflee(x,y,z);

gee = @(x,y,z) x.^2-y.^2+2*z.^2-2;
gradgee = @(x,y,z) grad_objfun2Dgee(x,y,z);

hee = @(x,y,z) 2*x.^2+y.^2-z.^2;
gradhee = @(x,y,z) grad_objfun2Dhee(x,y,z);


[bituminous1,bituminous2,bituminous3] = newton3D_exact(flee,gradflee,gee,gradgee,hee,gradhee,x0,y0,z0); 
disp('The roots for the system in problem 4 part b are the following:')
disp(bituminous1)
disp(bituminous2)
disp(bituminous3)

%% Work for plotting Bessel function
%Since the plotting is meant to just estimate the root, This code for
%plotting the Bessel function was taken directly from the MATLAB help page for
%besselj and modified slightly
numbs = 0:0.1:20;


J = zeros(1,201);
for i = 0
    J(i+1,:) = besselj(i,numbs);
end

plot(numbs,J)
grid on
title('Bessel Functions','interpreter','latex')
xlabel('z','interpreter','latex')
ylabel('$J_\nu(z)$','interpreter','latex')

%% Here the Code ends and the functions are defined below

function [root,it,success]=newton_approx(f,e,x0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f (input as a handle)
% given a function which computes the derivative

% Error checking of input
narginchk(3,6);   %check for correct number of inputs to function
if (nargin<4)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<5)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<6)
    verbose=false;
end %if

% Make sure we don't start at an inflection point with zero derivative
%{
if (abs(fprime(x0))<tol)
    warning(' Attempting to start Newton iterations near an inflection point, you may wish to restart with a different guess...');
    x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
end %if
%}

% Newton iterations
it=1;
root=x0;
fval=f(root);
converged=false;
while(~converged && it<=maxit)
    derivative=(f(root+e)-f(root))/e;
    if (abs(derivative)<100*tol)    %this (inflection point) will end up kicking the root really far away...
        converged=false;
        warning(' Derivative close to zero, terminating iterations with failed convergence... ');
        break;
    else
        root=root-fval./derivative;    % update root estimate
        fval=f(root);                  % see how far off we are from zero...
        if (verbose)
            fprintf(' iteration: %d; root:  %f + %f i; function value: %f, derivative:  %f \n',it,real(root),imag(root),fval,derivative);
        end %if
        it=it+1;
        converged=abs(fval)<tol;
    end %if
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function

function [root,it,success]=newton_exact(f,fprime,x0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f (input as a handle)
% given a function which computes the derivative

%% Error checking of input
narginchk(3,6);   %check for correct number of inputs to function
if (nargin<4)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<5)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<6)
    verbose=false;
end %if


%% Make sure we don't start at an inflection point with zero derivative
if (abs(fprime(x0))<tol)
    warning(' Attempting to start Newton iterations near an inflection point, you may wish to restart with a different guess...');
    x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
end %if


%% Newton iterations
it=1;
root=x0;
fval=f(root);
converged=false;
while(~converged && it<=maxit)
    derivative=fprime(root);
    if (abs(derivative)<100*tol)    %this (inflection point) will end up kicking the root really far away...
        converged=false;
        warning(' Derivative close to zero, terminating iterations with failed convergence... ');
        break;
    else
        root=root-fval./derivative;    % update root estimate
        fval=f(root);                  % see how far off we are from zero...
        if (verbose)
            fprintf(' iteration: %d; root:  %f + %f i; function value: %f, derivative:  %f \n',it,real(root),imag(root),fval,derivative);
        end %if
        it=it+1;
        converged=abs(fval)<tol;
    end %if
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function

function [rootx,rooty,it,success]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f,g (input as a handle)
% given a function which computes the derivative
% 
% requires:  Gauss_elim.m and backsub.m


%% Need access to linear algebra routines for solves at each iteration
%addpath ../linear_algebra;
%I commented this out because it was throwing a warning since I dont have a
%linear algebra folder I just keep the Gauss elim and backsub in the same
%folder

%% Error checking of input
narginchk(6,9);   %check for correct number of inputs to function
if (nargin<7)
    maxit=100;       %maximum number of iterations allowed
end %if
if (nargin<8)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<9)
    verbose=false;
end %if


%% Make sure we don't start at an inflection point with zero derivative
[gradfx,gradfy]=gradf(x0,y0);
[gradgx,gradgy]=gradg(x0,y0);
if (abs(min([gradfx,gradfy,gradgx,gradgy]))<tol)    %this needs to really check inflection vs. saddle point; will fix later
    warning(' Attempting to start Newton iterations near an inflection/saddle point, you may wish to restart with a different guess...');
    x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
    y0=y0+1;
end %if


%% Newton iterations
it=1;
rootx=x0;
rooty=y0;
fval=f(rootx,rooty);
gval=g(rootx,rooty);
converged=false;
while(~converged && it<=maxit)
    [gradfx,gradfy]=gradf(rootx,rooty);
    [gradgx,gradgy]=gradg(rootx,rooty);
    A=[gradfx,gradfy; ...
        gradgx,gradgy];
    
    fvec=[fval;gval];
    [Amod,ord]=Gauss_elim(A,-1*fvec);
    dxvec=backsub(Amod(ord,:));
    detA=prod(diag(Amod(ord,:)));
    if (abs(detA) < 1e-6)
        error(' Ended up at a point where Newton iteration is singular, try a different starting point')
    end %if
    
    rootx=rootx+dxvec(1);
    rooty=rooty+dxvec(2);
    fval=f(rootx,rooty);                  % see how far off we are from zero...
    gval=g(rootx,rooty);
    if (verbose)
        fprintf(' iteration: %d; x:  %f + %f i; y:  %f + %f i; f: %f, g:  %f\n',it,real(rootx),imag(rootx),real(rooty),imag(rooty),fval,gval);
    end %if
    it=it+1;
    converged=abs(fval)+abs(gval)<tol;
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function

function [rootx,rooty,rootz,it,success]=newton3D_exact(f,gradf,g,gradg,h,gradh,x0,y0,z0,maxit,tol,verbose)

% root=newton_exact(f,fprime)
%
% finds a set of roots corresponding to the function f,g (input as a handle)
% given a function which computes the derivative
% 
% requires:  Gauss_elim.m and backsub.m


%% Need access to linear algebra routines for solves at each iteration
%addpath ../linear_algebra;
%I commented this out because it was throwing a warning since I dont have a
%linear algebra folder I just keep the Gauss elim and backsub in the same
%folder

%% Error checking of input
narginchk(9,12);   %check for correct number of inputs to function
if (nargin<10)
    maxit=10000;       %maximum number of iterations allowed
end %if
if (nargin<11)
    tol=1e-6;        %how close to zero we need to get to cease iterations
end %if
if (nargin<12)
    verbose=false;
end %if


%% Make sure we don't start at an inflection point with zero derivative
[gradfx,gradfy,gradfz]=gradf(x0,y0,z0);
[gradgx,gradgy,gradgz]=gradg(x0,y0,z0);
[gradhx,gradhy,gradhz]=gradh(x0,y0,z0);
if (abs(min([gradfx,gradfy,gradfz,gradgx,gradgy,gradgz,gradhx,gradhy,gradhz]))<tol)    %this needs to really check inflection vs. saddle point; will fix later
    warning(' Attempting to start Newton iterations near an inflection/saddle point, you may wish to restart with a different guess...');
    x0=x0+1;   %bump the guess a ways off of initial value to see if we can get anything sensible
    y0=y0+1;
    z0=z0+1;
end %if


%% Newton iterations
it=1;
rootx=x0;
rooty=y0;
rootz=z0;
fval=f(rootx,rooty,rootz);
gval=g(rootx,rooty,rootz);
hval=h(rootx,rooty,rootz);
converged=false;
while(~converged && it<=maxit)
    [gradfx,gradfy,gradfz]=gradf(rootx,rooty,rootz);
    [gradgx,gradgy,gradgz]=gradg(rootx,rooty,rootz);
    [gradhx,gradhy,gradhz]=gradh(rootx,rooty,rootz);
    A=[gradfx,gradfy,gradfz; ...
        gradgx,gradgy,gradgz; ... 
        gradhx,gradhy,gradhz];
    
    fvec=[fval;gval;hval];
    [Amod,ord]=Gauss_elim(A,-1*fvec);
    dxvec=backsub(Amod(ord,:));
    detA=prod(diag(Amod(ord,:)));
    if (abs(detA) < 1e-6)
        error(' Ended up at a point where Newton iteration is singular, try a different starting point')
    end %if
    
    rootx=rootx+dxvec(1);
    rooty=rooty+dxvec(2);
    rootz=rootz+dxvec(3);
    fval=f(rootx,rooty,rootz);                  % see how far off we are from zero...
    gval=g(rootx,rooty,rootz);
    hval=h(rootx,rooty,rootz);
    if (verbose)
        fprintf(' iteration: %d; x:  %f + %f i; y:  %f + %f i; f: %f, g:  %f\n',it,real(rootx),imag(rootx),real(rooty),imag(rooty),fval,gval);
    end %if
    it=it+1;
    converged=abs(fval)+abs(gval)+abs(hval)<tol;
end %while
it=it-1;

if (~converged)
    warning('Used max number of iterations, or derivative near zero...')
    success=false;
else
    success=true;
end %if

end %function


function [fx,fy]=grad_objfun2Df(x,y)

fx=2.*x-2;
fy=2.*y-1;

end %function

function [gx,gy]=grad_objfun2Dg(x,y)

gx=0.5.*x;
gy=2.*y;

end %function

function [fleex,fleey,fleez]=grad_objfun2Dflee(x,y,z)

fleex=2.*x;
fleey=2.*y;
fleez=2.*z;

end %function

function [geex,geey,geez]=grad_objfun2Dgee(x,y,z)

geex=2.*x;
geey=-2.*y;
geez=4.*z;

end %function

function [heex,heey,heez]=grad_objfun2Dhee(x,y,z)

heex=4.*x;
heey=2.*y;
heez=-2.*z;

end %function
