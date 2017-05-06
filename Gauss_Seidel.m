function[Solution,iteration_table,Error_val]=Gauss_Seidel(U,Errorinput,F)
% This function aims to solve any 2D PDE equation in the Gausse Seidel
% method. I will use my Poisson Equation for my project as my input values
% of U and F. This function is unique that it solves for each value in the
% form of matrix multiplication style Ax=b, where A=U inut matrix, b=F of
% input column vector, and x is each iteration of individual element U.
A=U;
b=diag(F);
n=length(b);
% Initial values of X assume set-up to be 0
X=zeros(n,1);
Errorvalue=1;

iteration=0;
while (Errorvalue)>Errorinput 
    iteration=iteration+1; %making as many iterations until it fals within the error given
    Z=X; 
    for i=1:n
        j=1:n;
        j(i)=[]; % This empties out the part where there are no solutions and coefficients are not needed.
        Xtemp=X; % Here Xtemp takes the same values as X. 
        Xtemp(i)=[];% This reduce the Xtemp vector and eliminates spots with no answers.
        X(i)=(b(i)-sum(A(i,j)*Xtemp))/(A(i,i)); % The Gausse Seidel equation. 
    end 
    Xsol(:,iteration)=X;
    % The error for the first unknown is usually larger than the rest of
    % the other errors so it will be used as the error value. 
    Errorvalue=abs((X(2,1)-Z(2,1))/(X(2,1)));
end 
iteration_table=[1:iteration;Xsol]';
Solution=X;
Error_val=Errorvalue;