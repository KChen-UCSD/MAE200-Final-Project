function [X] = DRE(E,A,BC,XP_0,T,R,Q,ricatti)
% Function for solving DRE using RK4
% Credits should be given to Muhan Zhao and Isaac Need for this code

% For X, ricatti = 1 (controller)
% For P, ricatti = 0 (estimator)

X = zeros(size(A,1), size(A,2), size(A,3));
X(:,:,end) = XP_0;
% timestep
h = T/(size(A,3)-1);

% Estimator case set-up
if ricatti == 0
    X(:,:,1) = XP_0;
    X(:,:,end) = X(:,:,2);
end

% For controller case, 
if ricatti == 1
    for i = size(A,3):-1:2
        f1 = RHS_X(X(:,:,i),E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        f2 = RHS_X(X(:,:,i)-f1/2*h,E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        f3 = RHS_X(X(:,:,i)-f2/2*h,E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        f4 = RHS_X(X(:,:,i)-f3*h,E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        X(:,:,i-1) = X(:,:,i) - h*(f1/6 + f2/3 + f3/3 + f4/6);
    end
% For estimator case
elseif ricatti == 0
    for i = 1:size(A,3)-1
        f1 = RHS_P(X(:,:,i),E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        f2 = RHS_P(X(:,:,i)-f1/2*h,E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        f3 = RHS_P(X(:,:,i)-f2/2*h,E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        f4 = RHS_P(X(:,:,i)-f3*h,E(:,:,i),A(:,:,i),BC(:,:,i),R,Q);
        X(:,:,i+1) = X(:,:,i) + h*(f1/6 + f2/3 + f3/3 + f4/6);
    end
end

function dX = RHS_X(X,E,A,B,R,Q)
dX = -inv(E')*A'*X - X*A*inv(E) + X*B*inv(R)*B'*X - inv(E')*Q*inv(E);
end

function dP = RHS_P(P,E,A,C,R,Q)
dP = inv(E')*A'*P + P*A*inv(E) - P*C'*inv(R)*C*P + inv(E')*Q*inv(E);
end

end % End of function DRE


