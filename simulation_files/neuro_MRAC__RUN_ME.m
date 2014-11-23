clear; clc;


nnc  = ssNN('NN_WEIGHTS',[8,80,1],1.0,0);

pend = ssPendulum(...
    0.2,    ... % mass 1
    0.2,    ... % mass 2
    0.4,    ... % len. 1
    0.4,    ... % len. 2
    0.5,    ... % motor gain
    0.01,   ... % motor back emf
    0.001,  ... % joint 1 damping friction
    0.001,  ... % joint 2 damping friction
    0.001,  ... % disturbance amplitude
    0.001, ... % time step
    0.1     ...
);

viz_on          = 1;
lr              = 0.3;
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

for idx = 1:end_generation
    
    for jdx = 0:1:2
        pend    = pend.init(-ang_generation(idx),0);
        pend.t  = 0;
        while(pend.t<end_time)
            
            X_NL            = [pend.q_S ; pend.dq_S];
            X_L             = pend.z;

            meta_z          = [ X_NL; X_NLp; X_L; X_Lp ];
            u_des           = H*meta_z/K_G;

            if  jdx==2      
                D_input     = K_L*X_NL + K_N*K_G*u_des;
            elseif  jdx==1
                nnc         = nnc.activate([X_NL; X_L]/K_G/K_N);
                D_input     = K_L*X_NL + K_N*K_G*nnc.output();
            else
                D_input     = K_L*X_NL;
            end
            
            L_input         = K_L*X_L;
            pend            = pend.sensorUpdate();

            pend            = pend.linStep( L_input );
            pend            = pend.dynStep( D_input );
            pend.t          = pend.t + pend.ts;

            if jdx==1
                nnc         = nnc.reweight(lr,-u_des);
            end
            
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
        
        if  jdx == 2
            name = ['Set2_NN_Train', num2str(idx)];
            pend.dump(name)
            print('-dpng','-r300',[name,'.png'])
        elseif jdx==1
            name = ['Set2_NN_On', num2str(idx)];
            pend.dump(name)
            print('-dpng','-r300',[name,'.png'])
        else
            name = ['Set2_NN_Off', num2str(idx)];
            pend.dump(name)  
            print('-dpng','-r300',[name,'.png'])
        end
    end
    
end