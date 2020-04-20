%% Evan Dort EP501 HW3 
clc
clear
%% Problem 1
%% a and b
% Perform least squares fit poly of arbitrary order

load('test_lsq.mat')
n = 2;
polymat = zeros(n+1,n+1);
xymat = zeros(n+1,1);
Num=size(x,1);
anseqn = zeros(Num,1);
LINEQN = zeros(Num,1);
MATeqn = zeros(Num,1);

% calculate equation matrix
for i = 1:n+1
    power=i-1;
    for j = 1:n+1
        for sum = 1:Num
            polymat(i,j) = polymat(i,j)+x(sum)^power;
        end %for
        power = power+1;
    end %for  
end %for

%calculate xy matrix
for i = 1:n+1
    power=i-1;
    for sum = 1:Num
            xymat(i) = xymat(i)+ynoisy(sum)*x(sum)^power;
    end %for
end %for

solsmat=inv(polymat)*xymat; %coefficient matrix
power=0;

for point = 1:Num % for loop for arbitrary power
    for coeff = 1:n+1
        anseqn(point,1) = anseqn(point)+x(point)^(power)*solsmat(coeff);
        power = power+1;
    end %for
    power = 0;
end %for

power=0;
for point = 1:Num %same for loop as above but for only a line since the assignment asks for a line
    for coeff = 1:2
        LINEQN(point,1) = LINEQN(point)+x(point)^(power)*solsmat(coeff);
        power = power+1;
    end %for
    power = 0;
end %for

%Built-in MATLAB
MATpoly1 = polyfit(x,ynoisy,n);
MATpoly2 = polyval(x,ynoisy,n);

for point = 1:Num % for loop for arbitrary power from Built-in MATLAB
    for coeff = 1:n+1
        MATeqn(point,1) = MATeqn(point)+x(point)^(power)*MATpoly1(coeff);
        power = power+1;
    end %for
    power = 0;
end %for


%Plot the data
figure(1)
plot(x,ynoisy) %data line
hold on
plot(x,anseqn,'LineWidth', 2) %my poly fit line
hold on
plot (x,LINEQN, '-*','MarkerIndices',1:20:length(x)) %linear equation request as part of HW
hold off
hold on
plot (x,MATeqn, 'o','MarkerIndices',1:30:length(x))
hold off
hold on
plot (x,MATeqn, 'diamond','MarkerIndices',1:30:length(x))
hold off
xlabel('x')
ylabel('ynoisy')
legend({'Noisy Data', 'My Least Square Arbitrary Power', 'My Least Square Linear', 'MATLAB polyfit', 'MATLAB polyval'},'Location','northwest')

disp('Quadratic looks to be a better fit than linear.')

%% Part c and d
%instead of setting up the code to run multiple fits, I ran mine for each
%different power value and then gave my answer as a display at the end

uncert = std(ynoisy);

for i = 1:Num
    
    summation = (ynoisy(i)-anseqn(i))^2/uncert^2;
    
end %for

nu = (Num-1);

chi2 = sqrt((1/nu)*summation);
Final1 = sprintf('For linear, my \x03A7^2 value is 0.0368.');
Final2 = sprintf('For quadratic, my \x03A7^2 value is 0.0066.');
Final3 = sprintf('For cubic, my \x03A7^2 value is 0.0102.');
Syzygy = sprintf('Quadratic has the smallest \x03A7^2 value. Thus, the best fit is quadratic.');
disp(Final1)
disp(Final2)
disp(Final3)
disp(Syzygy)


%% Problem 2

load('test_interp.mat')

%% Part a
for i = 1:512
    xp = xgi(i);
    incx = size(xg,2);

    for j = 1:incx
        if xp > xg(j) && xp < xg(j+1)
            indexx(i) = j;
        end
    end %for

end %for

%% Part b

for i = 1:512
    yp = ygi(i);
    incy = size(yg,2);

    for j = 1:incy
        if yp > yg(j) && yp < yg(j+1)
            indexy(i) = j;
        end
    end %for

end %for

%% Part c

for val = 2:510

        outs(1,1) = f2D((indexx(val)),(indexy(val)));
        outs(2,1) = f2D((indexx(val)+1),(indexy(val)));
        outs(3,1) = f2D((indexx(val)),(indexy(val)+1));
        outs(4,1) = f2D((indexx(val)+1),(indexy(val)+1));
        
        xvec = [xg(indexx(val)) xg(indexx(val)+1) xg(indexx(val)) xg(indexx(val)+1)];
        yvec = [yg(indexy(val)) yg(indexy(val)) yg(indexy(val)+1) yg(indexy(val)+1)];
        
        M = [ones(4,1),xvec(:),yvec(:),xvec(:).*yvec(:)];

        avec = inv(M)*outs;   
        
for i = 1:4  
finterpmanual(val-1) = avec(1)+avec(2)*xvec(i)+avec(3)*yvec(i)+avec(4)*xvec(i)*yvec(i);
end

end


%% Part d

[X,Y] = meshgrid(xg,yg);

finterp = interp2(X,Y,f2D,xgi,ygi);

for counts = 2:510
    xworks(counts-1) = xgi(counts);
end %for

figure(2)
plot(xworks,finterpmanual)
figure(3)
plot(xgi,finterp)
xlim([-6 6])
title('Built-in MATLAB')



