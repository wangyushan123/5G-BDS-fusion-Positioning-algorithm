function [xAMP, PAMP] = algo_3_ivb_ampKF(x,p,xls,F,H,Q_small,R,G,dt,v1,v2)
 

    soft = @(v,lambda) sign(v).*max(abs(v)-lambda,0);% Soft threshold function
    xAMP = x; 
    PAMP = p;
    alpha=1.36;
    z=xls;
    x_pred = F*xAMP;%(:,k-1)
    r = z - H*x_pred;
    a_raw = r/(0.5*dt^2+eps);       % Reverse push-out acceleration
    sigma_r=std(a_raw);
    lambda = alpha*sigma_r;% Threshold value
    a_hat = soft(a_raw, lambda);    % sparse recovery?
    a_hat=[zeros(3,1);a_hat];

    x_pred = x_pred- G*a_hat ;      % correcting prediction

    P_pred = F*PAMP*F' + Q_small+diag((sqrt(lambda)*0.5*dt^2)^2*ones(6,1));%0.4
    xkk=x_pred;
    xk1k=xkk;
    Pkk=P_pred;
    Pk1k=Pkk;
    nz=size(z,1);

    nx=size(xkk,1);
%     for i=1:N
    %%%%%%%Calculate auxiliary parameter 
    Dk1=(xkk-xk1k)*(xkk-xk1k)'+Pkk;
    
    Dk2=(z-H*xkk)*(z-H*xkk)'+H*Pkk*H';
    
    gama1=trace(Dk1*myinv(Pk1k));

    gama2=trace(Dk2*myinv(R));

    lamda1=(v1+nx)/(v1+gama1);

    lamda2=(v2+nz)/(v2+gama2);
    
    %%%%%%%Update the distribution of state vector
    D_Pk1k=Pk1k/lamda1;

    D_R=R/lamda2;

    zk1k=H*xk1k;
    
    Pzzk1k=H*D_Pk1k*H'+D_R;
    
    Pxzk1k=D_Pk1k*H';

    % update
    K = Pxzk1k*myinv(Pzzk1k);
    xAMP = x_pred + K*(z-H*x_pred);
    PAMP = D_Pk1k-K*H*D_Pk1k;



end
