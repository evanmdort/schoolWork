%%Evan Dort EP501 HW 1
clc
clear
%% Problem 1 parts a-b
load testproblem.mat

[m,n]=size(A);
Aorig=zeros(8);%makes a copy of the matrices given in test problem so I can reference when I need to. I know this may be a bit uncessary but I prefer to have this.
borig=zeros(8,1);
b2orig=zeros(8,1);
b3orig=zeros(8,1);
solb=zeros(8,1);
borig=b;
b2orig=b2;
b3orig=b3;
solb=b;
Aorig=A;
i=1;%setup of counters
j=1;
A= [A solb];
    for d=j:m-1
        i=i+(j-1);%This increments it to follow the desired elim element
     for c=i:m-1
        A(i+1,:)=A(i+1,:)-(A(i+1,j))/(A(j,j))*A(j,:);%This is the iteration setup so that it runs through each element for elimination
        i=i+1;
     end
     j=j+1;
     i=1;
    end
    A(abs(A)<1e-14)=0;% I was having rounding error where it would display -0.0000 and in the viewer it showed e-16 so I have this for tolerance
    
 
    
    yee=backsub(A);
    [y1 y2]=Gauss_elim(Aorig,borig);
    z=inv(Aorig)*borig;
    y=y1(y2,:);
    gausans=backsub(y);
    disp('Problem 1')
    disp('My forward elimination')
    disp(yee)
    disp('MATLAB built-in')
    disp(z)
    disp('Gausian Elimination')
    disp(gausans)
    %From here if x,y, and z are viewed it shows that they give the same
    %solution, meaning the forward elim I made works
%% Problem 1 Part c-d    

load lowertriang_testproblem.mat
i=2;
[m n]=size(L);
sums=0;
q=zeros(n,1);
ANS=inv(L)*bL;
q(1)=bL(1,1)/L(1,1);
for i=2:n
    for j=1:i-1
    sums=sums+L(i,j)*q(j);
    end
    q(i)=1/L(i,i)*(bL(i)-sums);
    sums=0;%This portion of code was constructed from general eqn for forward sub in http://webhome.auburn.edu/~tamtiny/lecture%203.pdf
 
end
%This gave the same results as the matlab as shown by the variable ANS
disp('Lower Triangluar MATLAB')
disp(ANS)
disp('My Forward Substitution')
disp(q)

%% Problem 2 Part a-d

[m,n]=size(Aorig);
I=eye(m);
SolA=[Aorig I];
i=1;
j=1;
 for d=j:m-1
        i=i+(j-1);%This increments it to follow the desired elim element
     for c=i:m-1
        SolA(i+1,:)=SolA(i+1,:)-(SolA(i+1,j))/(SolA(j,j))*SolA(j,:);%This is the iteration setup so that it runs through each element for elimination, making it upper triangular
        i=i+1;
     end
     j=j+1;
     i=1;
    end
SolA(abs(SolA)<1e-14)=0;% I was having rounding error where it would display -0.0000 and in the viewer it showed e-16 so I have this for tolerance
i=n;
j=n;
 while j>1 
     %This increments it to follow the desired elim element
     while i>1
        SolA(i-1,:)=SolA(i-1,:)-(SolA(i-1,j))/(SolA(j,j))*SolA(j,:);%This is the iteration setup so that it runs through each element for elimination, making it lower triangular
        i=i-1;
     end
     j=j-1;
     i=j;
 end

SolA(abs(SolA)<1e-14)=0;% I was having rounding error where it would display -0.0000 and in the viewer it showed e-16 so I have this for tolerance

j=1;

while j<=n
SolA(j,:)=SolA(j,:)/SolA(j,j);
j=j+1;
end

Sol=inv(Aorig);
disp('Problem 2')
disp('MATLAB built-in Inverse')
disp(Sol)
disp('My solution with right side of Matrix as Inverse')
disp(SolA)

%% Problem 3

[m,n]=size(Aorig);
Lsol=zeros(m);
Usol=Aorig;
bp=b;
i=1;
j=1;
g=zeros(m,1);
for d=j:m-1
        i=i+(j-1);%This increments it to follow the desired elim element
     for c=i:m-1
        num=(Usol(i+1,j));
        denom=(Usol(j,j));
        Usol(i+1,:)=Usol(i+1,:)-(Usol(i+1,j))/(Usol(j,j))*Usol(j,:);%This is the iteration setup so that it runs through each element for elimination, making it upper triangular
        Lsol(i+1,j)=num/denom;
        i=i+1;
     end
     j=j+1;
     i=1;
end
j=1;
while j<9
    Lsol(j,j)=1;
    j=j+1;
end

g(1)=bp(1,1)/Lsol(1,1);
for i=2:n
    for j=1:i-1
    sums=sums+Lsol(i,j)*g(j);
    end
    g(i)=1/L(i,i)*(bp(i)-sums);
    sums=0;%This portion of code was constructed from general eqn for forward sub in http://webhome.auburn.edu/~tamtiny/lecture%203.pdf
 
end


Usol(abs(Usol)<1e-14)=0;


LUsol=[Usol g];
Ansol=backsub(LUsol);
disp('Problem 3')
disp('Backsub of test problem')
disp(Ansol)

%% Problem 3 c-d

B=[b b2 b3];
[row, column]=size(B);
CombSol=zeros(row, column);
i=1;
while i <= column
CombSol(1,i)=B(1,i)/Lsol(i,i);
i=i+1;
end

for k=1:column
    for i=2:n
     for j=1:i-1
        sums=sums+Lsol(i,j)*CombSol(j,k);
     end
        CombSol(i,k)=1/Lsol(i,i)*(B(i,k)-sums);
        sums=0;%This portion of code was constructed from general eqn for forward sub in http://webhome.auburn.edu/~tamtiny/lecture%203.pdf
 
    end
end


LU1=[Usol CombSol(:,1)];%solution for all 3 matricies
LU2=[Usol CombSol(:,2)];
LU3=[Usol CombSol(:,3)];
LUANS1=backsub(LU1);
LUANS2=backsub(LU2);
LUANS3=backsub(LU3);
disp('Solution of system of equations')
disp(LUANS1)
disp(LUANS2)
disp(LUANS3)

%% Problem 4 I copied and pasted the iterative part of the Jacobi 
%on Git into this so I can keep all my changes on one file
%define variables
load('iterative_testproblem.mat')
tol=1*10^-6;
[m n]=size(Ait);
x0=zeros(m,1);
w=1.4;
%% Setup iterations
maxit=100;    %max number of iterations
n=size(Ait,1);  %system size
residual=10*ones(n,1);
difftot=1e3+tol;   %max sure we enter iterations
x=x0;


%% Perform iterations
it=1;
while(difftot>tol && it<=maxit)
    difftotprev=difftot;
    resprev=residual;
    xprev=x;
    for i=1:n
        residual(i)=bit(i);
        for j=1:i-1
            residual(i)=residual(i)-Ait(i,j)*x(j);
        end %for
        for j=i:n
            residual(i)=residual(i)-Ait(i,j)*xprev(j);
        end %for
        x(i)=xprev(i)+(w)*residual(i)/Ait(i,i);
    end %for
    difftot=sum(abs(residual-resprev));
    
    
    if (difftot>difftotprev)
        error('Solution appears to be diverging, check diagonal dominance...')
    end %if
    it=it+1;
end %while

nit=it-1;
if (n==maxit)
    warning('Solution may not have converged fully...')
end %if

Zet=inv(Ait)*bit;%Comparing this to the variable x, they are the same, so my changes worked.
disp('Problem 4')
disp('MATLAB built-in')
disp(Zet)
disp('My Jacobi Iterative')
disp(x)

%% Problem 5 As in problem 4 I copied and pasted the function
%on Git into this so I can keep all my changes on one file

% [Amod,ord]=elim(A,b,verbose)
%
% This function perform elimination with partial pivoting and scaling as
% described in Section 1.3.2 in the Hoffman textbook (viz. it does Gaussian
% elimination).  Note that the ordering which preserves upper triangularity
% is stored in the ord output variable, such that the upper triangular output
% is given by row-permuted matrix Amod(ord,:).  The verbose flag can be set to
% true or false (or omitted, default=false) in order to print out what the algirthm
% is doing for each elimination step.

%Allocation of space and setup
Amod=cat(2,A,b);          %make a copy of A and modify with RHS of system
n=size(A,1);              %number of unknowns
ord=(1:n)';               %ord is a mapping from input row being operated upon to the actual row that represents in the matrix ordering

%Elimination with scaled, partial pivoting for matrix Amod; note all row
%indices must be screen through ord mapping.

for ir1=1:n-1    
    %check scaled pivot elements to see if reordering should be done
    pivmax=0;
    ipivmax=ir1;      %max pivot element should never be higher than my current position
    for ipiv=ir1:n    %look only below my current position in the matrix
        pivcurr=abs(Amod(ord(ipiv),ir1))/max(abs(Amod(ord(ipiv),:)));      %note that columns never get reordered...
        if (pivcurr>pivmax)
            pivmax=pivcurr;
            ipivmax=ipiv;     %this stores the index into ord for row having largest pivot element
        end %if
    end %for
    
    %reorder if situation calls for it
    if (ipivmax ~= ir1)
        itmp=ord(ir1);
        ord(ir1)=ord(ipivmax);
        ord(ipivmax)=itmp;
        
        if (verbose)
            disp('Interchanging rows:  ');
            disp(itmp);
            disp(' and:  ');
            disp(ord(ir1));
            disp('Current matrix state after interchange:  ');
            disp(Amod(ord,:));
        end %if
    end %if
    
    %perform the elimination for this row, former references to ir1 are now
    %mapped through the ord array
    for ir2=ir1+1:n
        fact=Amod(ord(ir2),ir1);
        Amod(ord(ir2),ir1:n+1)=Amod(ord(ir2),ir1:n+1)-fact/Amod(ord(ir1),ir1).*Amod(ord(ir1),ir1:n+1);    %only need columns ahead of where we are in matrix
    end %for
end
DetA=1;
for i=1:n
    DetA=DetA*Amod(i,i);
end

AnsFive=det(Aorig);
disp('Problem 5')
disp('MATLAB built-in')
disp(AnsFive)
disp('My determinant')
disp(DetA)
