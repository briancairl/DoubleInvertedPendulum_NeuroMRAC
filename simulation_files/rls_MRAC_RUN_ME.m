clear; clc;

pend = ssPendulum(...
    0.2,    ... % mass 1
    0.2,    ... % mass 2
    0.4,    ... % len. 1
    0.4,    ... % len. 2
    0.5,    ... % motor gain
    0.01,   ... % motor back emf
    0.001,  ... % joint 1 damping friction
    0.001,  ... % joint 2 damping friction
    0.01,  ... % disturbance amplitude
    0.001, ... % time step
    0.1     ...
);


lr              = 0.3;
viz_on          = 0;
end_time        = 4;
end_generation  = 10;
ang_generation  = linspace(pi/16,pi/4,end_generation).*(-1).^(1:end_generation);

D1              = [ zeros(2), eye(2)  ; zeros(2), eye(2)/pend.ts ];      
D2              = [ zeros(2), zeros(2); zeros(2), eye(2)/pend.ts ];
G               = pend.B;
K_G             = norm([pend.B;pend.B]);
G_n             = G/K_G;
X_Lp            = zeros(4,1);
X_NLp           = zeros(4,1);
K_N             = 5;
K_L             =-lqr(pend.A,pend.B,eye(4),2);
M               = D1 - pend.A - pend.B*K_L;
H               = [G.'*M,-G.'*M,-G.'*D2,G.'*D2]/K_G^2;


% RLS
n = 100;
W = zeros(1,8);
z = 0.00001;
X = zeros(8,n);
P = zeros(n,1);


for idx = 1:end_generation
    set(gcf,'NAME',['q(0) =',num2str(ang_generation(idx))])
    pend    = pend.init(-ang_generation(idx),0);
    pend.t  = 0;
    while(pend.t<end_time)

        X_NL            = [pend.q_S ; pend.dq_S];
        X_L             = pend.z;

        meta_z          = [ X_NL; X_NLp; X_L; X_Lp ];
        u_des           = H*meta_z/K_G;

        P(end+1,1)      = u_des;
        P(1)            = [];
        X(1:4,end+1)    = X_NL;
        X(5:8,end  )    = X_L;
        X(:,1)          = [];            
        D_input         = K_L*X_NL - K_N*K_G*W*X(:,end);
        
        L_input         = K_L*X_L;
        pend            = pend.sensorUpdate;
        pend            = pend.linStep( L_input );
        pend            = pend.dynStep( D_input );
        pend.t          = pend.t + pend.ts;

        W               = W + z*(W - transpose( inv(X*X.' + z*eye(8) )*X*P ));

        X_Lp            = X_L;
        X_NLp           = [pend.q_S;pend.dq_S];

        pend            = pend.sample();

        if viz_on
            cla
            plot(pend,'NLM')
            pause(pend.ts)
        end
        
        if(pend.status())
            disp('Failure...')
            break
        end
    end
        pend.plot_states(); 
        pause(1);
        
        name = ['Set2_LS', num2str(idx)];
        pend.dump(name)  
        print('-dpng','-r300',[name,'.png'])
    
end