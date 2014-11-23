classdef ssPendulum
    
    properties
        m   = ones(2,1);
        l   = ones(2,1);
        kemf= 0.1;          % motor 1 back-emf
        p1  = 1;            % motor 1 current-torque gain
        n1  = 1;            % joint 1 friction
        n2  = 1;            % joint 2 friction
        M1  = 1;            % motor drive coeffs
        var = 0;            
        
        % Linearized Dynamics
        A       = zeros(4);
        B       = zeros(4,1);
        C       = eye(4,4);
        H       = eye(2);
        E       = eye(2);
        G       = eye(2);
        T       = eye(4);
        
        % State Information
        u_nL    = 0;
        u_L     = 0;
        z       = zeros(4,1);
        dz      = zeros(4,1);
        z0      = zeros(4,1);
        q       = zeros(2,1);
        dq      = zeros(2,1);
        ddq     = zeros(2,1);
        tau     = zeros(2,1);
        
        q_L     = zeros(2,1);
        dq_L    = zeros(2,1);
        q_S     = zeros(2,1);
        dq_S    = zeros(2,1);
        
        dq_S_LP = zeros(2,1);
        LP_T    = 1;
        
        
        % Log/Sim Data
        t   = 0;
        ts  = 1e-3;
        q_S_data = [];
        dq_S_data= [];
        q_L_data = [];
        dq_L_data= [];
        q_data   = [];
        dq_data  = [];
        tau_data = [];
        u_nL_data= [];
        u_L_data = [];
        plot_lims = 2;
    end
    
    
    methods
        
        function obj = ssPendulum(m1,m2,l1,l2,p1,kemf,n1,n2,var,ts,lp_T)
            obj.LP_T    = lp_T;
            obj.var     = var;
            obj.ts      = ts;
            obj.n1      = n1;
            obj.n2      = n2;
            obj.p1      = p1;
            obj.kemf    = kemf;
            
        	obj.m(1)    = m1; 
        	obj.m(2)    = m2; 
        	obj.l(1)    = l1; 
            obj.l(2)    = l2;
                        
            obj.H       = l_H(m1,m2,l1,l2);
            obj.G       = l_G(m1,m2,l1,l2);
            obj.E       = inv(obj.H);
            obj.M1      = (p1*kemf + n1);
            
            obj.T       = [ obj.E, zeros(2); zeros(2), obj.E ];

            
            % Construct the A matrix
            obj.A(1,3)  = 1;
            obj.A(2,4)  = 1;
            
            obj.A(3,1)  = -(obj.G(1,1)*obj.E(1,1) + obj.G(1,2)*obj.E(2,1));
            obj.A(4,1)  = -(obj.G(2,1)*obj.E(1,1) + obj.G(2,2)*obj.E(2,1));
            obj.A(3,2)  = -(obj.G(1,1)*obj.E(1,2) + obj.G(1,2)*obj.E(2,2));
            obj.A(4,2)  = -(obj.G(2,1)*obj.E(1,2) + obj.G(2,2)*obj.E(2,2));    
            obj.A(3,3)  = -obj.M1*obj.E(1,1);
            obj.A(4,3)  = -n2*obj.E(2,1);
            obj.A(3,4)  = -obj.M1*obj.E(1,2);
            obj.A(4,4)  = -n2*obj.E(2,2);
            
            obj.A       = (obj.T)*obj.A/obj.T; 
            
            % Construct the B matrix
            obj.B(3,1)  = p1;
            obj.B       = (obj.T)*obj.B;
            
            
            % Construct the C matrix
            obj.C       = eye(4);
            
            
            % Visualization Setup
            figure(1); 
            clf
            hold on
            obj.plot_lims = 1.5*sum(obj.l(1)+obj.l(2));
            set(gca,'XLIM',obj.plot_lims*[-1,1]/2);
            set(gca,'YLIM',obj.plot_lims*[ 0,1]);
            xlabel('X, m')
            ylabel('Y, m')
            
            axis square
            obj = obj.init(0,0);
            obj.plot('NLM');
        end
        
        
        function figure_set(obj)
            figure(1); 
            clf
            hold on
            obj.plot_lims = 1.5*sum(obj.l(1)+obj.l(2));
            set(gca,'XLIM',obj.plot_lims*[-1,1]/2);
            set(gca,'YLIM',obj.plot_lims*[ 0,1]);
            xlabel('X, m')
            ylabel('Y, m')
            axis square
        end
        
        
        
        function obj = dynStep(obj,u)
            obj.u_nL = u;
            
            % D Matrix
            nlD = nl_D(obj.m(1),obj.m(2),obj.l(1),obj.l(2),obj.q);
            
            % C Matrix
            nlC = nl_C(obj.m(1),obj.m(2),obj.l(1),obj.l(2),obj.q,obj.dq);
            
            % V Vector
            nlV = nl_V(obj.m(1),obj.m(2),obj.l(1),obj.l(2),obj.q);
            
            % Step...
            obj.tau(1)  =   -obj.dq(1)*obj.M1 + u*obj.p1;
            obj.tau(2)  =   -obj.dq(2)*obj.n2;
            
            obj.ddq     =    nlD\(obj.tau-nlC*obj.dq-nlV+ randn(2,1)*obj.var);
            obj.dq      =    obj.dq + obj.ddq*obj.ts;
            obj.q       =    obj.q  + obj.dq*obj.ts;
            
            obj.q(1)    =   piNpi(obj.q(1));
            obj.q(2)    =   piNpi(obj.q(2));
        end
        
        
        
        
        function obj = linStep(obj,u)
            obj.u_L = u;
            obj.dz  = (obj.A*obj.z + obj.B*u);
            obj.z   = obj.z + obj.dz*obj.ts;
            obj.q_L(1:2,1)   = obj.z(1:2,1);
            obj.dq_L(1:2,1)  = obj.z(3:4,1);
        end
        
        
        
        
        function obj = simStep(obj,u)
            obj     = linStep(obj,u);
            obj     = dynStep(obj,u);
            obj     = obj.sensorUpdate();
            obj.t   = obj.t + obj.ts;
        end
        
        
        
        function u_est = inputEstimate(obj,error)

            u_est   =  transpose(error.'*exp(obj.A*obj.t)*obj.B);
            
        end
        
        
        
        function obj = sensorUpdate(obj)
            dq_Sp       =  obj.dq_S;
            obj.q_S     =  obj.q  + randn(2,1)*1e-3;
            obj.dq_S    =  obj.dq + randn(2,1)*1e-6;
            
            obj.dq_S_LP =  obj.dq_S_LP + obj.LP_T*(obj.dq_S-dq_Sp);
           
        end
        
        
        
        
        function e = refError(obj)
            e   = obj.q_L-obj.q_S;
        end
       
        

        function obj = init(obj,q1,q2)
            obj.q_S_data = [];
            obj.dq_S_data= [];
            obj.q_L_data = [];
            obj.dq_L_data= [];
            obj.q_data   = [];
            obj.dq_data  = [];
            obj.tau_data = [];
            obj.u_nL_data= [];
            obj.u_L_data = [];
            
            obj.q_L = [q1;q2];
            obj.dq_L= zeros(2,1);
            obj.q   = [q1;q2];
            obj.dq  = zeros(2,1);
            obj.ddq = zeros(2,1);
            
            obj.z   = zeros(4,1);
            obj.z(1)= q1;
            obj.z(2)= q2;
            obj.z(3)= 0;
            obj.z(4)= 0;
            obj.z0  = obj.z;
        end
        
        
        
        
        function obj = sample(obj)
            obj.q_S_data(:,end+1)   = obj.q_S;
            obj.dq_S_data(:,end+1)  = obj.dq_S;
            obj.q_L_data(:,end+1)   = obj.q_L;
            obj.dq_L_data(:,end+1)  = obj.dq_L;
            obj.q_data(:,end+1)     = obj.q;
            obj.dq_data(:,end+1)    = obj.dq;
            obj.tau_data(:,end+1)   = obj.tau;
            obj.u_nL_data(:,end+1)  = obj.u_nL;
            obj.u_L_data(:,end+1)   = obj.u_L;
        end
        
        
        
        function rpvec = cPBH(obj)
            rpvec       = zeros(4,2);
            rpvec(:,1)  = eig(obj.A);
            for idx = 1:4
                rpvec(:,2)  = rank([obj.A-eye(4)*rpvec(idx,1), obj.B]);
            end
        end
        
        
        
        function rpvec = oPBH(obj)
            rpvec       = zeros(4,2);
            rpvec(:,1)  = eig(obj.A);
            for idx = 1:4
                rpvec(:,2)  = rank([obj.A-eye(4)*rpvec(idx,1); obj.C]);
            end
        end
        
       
        function stat_out = status(obj)
            stat_out = (abs(obj.q(1)) > 3*pi/4);
        end
        
        
        function plot(obj,spec)
            figure(1);
            
            q12 = sum(obj.q(1:2));
            q12L= sum(obj.q_L(1:2));
            

            cla
            hold on
            plot([-obj.plot_lims,obj.plot_lims],[0,0],'k:')
            
            if any(spec=='L')
            
            pL1L= obj.l(1)*[ -sin(obj.q_L(1));  cos(obj.q_L(1)) ];
            pL2L= pL1L+ obj.l(2)*[ -sin(q12L);  cos(q12L)       ];
                
            % The Linear Model Output
            plot([0,        pL1L(1)], [0,         pL1L(2)],'g.-')
            plot([pL1L(1),  pL2L(1)], [pL1L(2),   pL2L(2)],'g.-')
            
            end
            
            
            if any(spec=='N')
            
            pL1 = obj.l(1)*[  -sin(obj.q(1));   cos(obj.q(1))   ];
            pL2 = pL1 + obj.l(2)*[ -sin(q12);   cos(q12)        ];
                
            % Real Output
            plot([0,        pL1(1)],   [0,         pL1(2)],'k.-')
            plot([pL1(1),   pL2(1)],   [pL1(2),    pL2(2)],'k.-')
            
            end
            
            
            if any(spec=='M')           
                
            pM1 = obj.l(1)*[ -sin(obj.q(1));    cos(obj.q(1))   ]/2;
            pM2 = pL1 + obj.l(2)*[ -sin(q12);   cos(q12)        ]/2;
                
            % Masses
            plot(pM1(1),pM1(2),'rx',pM1(1),pM1(2),'ro')
            plot(pM2(1),pM2(2),'rx',pM2(1),pM2(2),'ro')
            
            end
        end
        
        
        
        function plot_states(obj)
            tspan = (1:numel(obj.q_data(1,:)))*obj.ts;
            subplot(4,1,1); plot( tspan,obj.q_data(1,:)   ,'b-', tspan,obj.q_L_data(1,:) ,'r-' ); ylabel('q_1,rad')
            set(gca,'XLIM',[tspan(1),tspan(end)]);
            set(gca,'YLIM',1.75*[-1,1]);
            subplot(4,1,2); plot( tspan,obj.q_data(2,:)   ,'b-', tspan,obj.q_L_data(2,:) ,'r-' ); ylabel('q_2,rad')
            set(gca,'XLIM',[tspan(1),tspan(end)]);
            set(gca,'YLIM',1.75*[-1,1]);
            subplot(4,1,3); plot( tspan,obj.dq_data(1,:)  ,'b-', tspan,obj.dq_L_data(1,:),'r-' ); ylabel('dq_1/dt,rad/s')
            set(gca,'XLIM',[tspan(1),tspan(end)]);
            set(gca,'YLIM',[-12,12]);
            subplot(4,1,4); plot( tspan,obj.dq_data(2,:)  ,'b-', tspan,obj.dq_L_data(2,:),'r-' ); ylabel('dq_2/dt,rad/s')
            set(gca,'XLIM',[tspan(1),tspan(end)]);
            set(gca,'YLIM',[-12,12]);
            xlabel('t,s')
        end
     

        function plot_input(obj)
            clf
            tspan = (1:numel(obj.q_data(1,:)))*obj.ts;
            plot( tspan,obj.u_nL_data(1,:)   ,'b-', tspan,obj.u_L_data(1,:) ,'r-' ); ylabel('q_1,rad')
            set(gca,'XLIM',[tspan(1),tspan(end)]);
            set(gca,'YLIM',20*[-1,1]);
        end

        
        
        
        function playback(obj,step,stop,varargin)
            if nargin == 4
                obj = load(obj,varargin{1});
            end
            figure_set(obj);
            hold on
            for idx = 1:step:numel(obj.q_data(1,:))
                obj.q   = obj.q_data(:,idx);
                obj.q_L = obj.q_L_data(:,idx);
                cla;
                plot(obj,'NLM');
                pause(obj.ts);
                
                if (idx > stop) && ( stop )
                    break;
                end
            end
            hold off
        end
        
        
        
        function dump(obj,savename)
            DATA            = zeros(16,numel(obj.q_S_data(1,:)));
            DATA(1:2  ,:)   = obj.q_S_data;
            DATA(3:4  ,:)   = obj.dq_S_data;
            DATA(5:6  ,:)   = obj.q_L_data;
            DATA(7:8  ,:)   = obj.dq_L_data;
            
            DATA(9:10 ,:)   = obj.q_data;
            DATA(11:12,:)   = obj.dq_data;
            DATA(13:14,:)   = obj.tau_data;
            DATA(15,:)      = obj.u_nL_data;
            DATA(16,:)      = obj.u_L_data;
            
            disp(['ssPendulum : Saved ',savename])
            save(savename,'DATA');
        end
        
        
        function obj = load(obj,savename)
            load(savename);
            disp(['ssPendulum : Loaded ',savename])
            obj.q_S_data    = DATA(1:2  ,:);
            obj.dq_S_data   = DATA(3:4  ,:);
            obj.q_L_data    = DATA(5:6  ,:);
            obj.dq_L_data   = DATA(7:8  ,:);
            
            obj.q_data      = DATA(9:10 ,:);
            obj.dq_data     = DATA(11:12,:);
            obj.tau_data    = DATA(13:14,:);
            obj.u_nL_data   = DATA(15,:);
            obj.u_L_data    = DATA(16,:);
        end
        
        
    end
    
end

function D = nl_D(m1,m2,l1,l2,q)
    b1 = l1/2;
    b2 = l2/2;
    D  = zeros(2,2);
    
    D(1,1) = m1*b1^2 + m2*(l1^2 + b2^2 + 2*l1*b2*cos(q(2)));
    D(1,2) = m2*(b2^2 + l1*b2*cos(q(2)));
    D(2,1) = D(1,2);
    D(2,2) = m2*b2^2;
end

function V  = nl_V(m1,m2,l1,l2,q)
    b1      = l1/2;
    b2      = l2/2;
    V       = zeros(2,1); 
    
    V(1,1)  =-9.81*(m1*b1+m2*l1)*sin(q(1)) - 9.81*m2*b2*sin(sum(q));
    V(2,1)  =-9.81*m2*b2*sin(sum(q));
end

function C  = nl_C(m1,m2,l1,l2,q,dq)

    b2      = l2/2;
    C       = zeros(2,2);  
    C(1,1)  =-l1*b2*m2*sin(q(2))*dq(2);
    C(1,2)  =-l1*b2*m2*sin(q(2))*(dq(1)+dq(2));
    C(2,1)  = l1*b2*m2*sin(q(2))*dq(1);
end

function H  = l_H(m1,m2,l1,l2)
    b1      = l1/2;
    b2      = l2/2;
    
    H       = zeros(2,2);  
    H(1,1)  = m1*b1^2 + m2*(l1^2 + b2^2 + 2*l1*b2);
    H(1,2)  = m2*(b2^2 + l1*b2);
    H(2,1)  = H(1,2);
    H(2,2)  = m2*b2^2;
end

function G  = l_G(m1,m2,l1,l2)
    b1      = l1/2;
    b2      = l2/2;
    
    G       = zeros(2,2);  
    G(1,1)  =-9.81*(m1*b1 + m2*l1 + m2*b1);
    G(1,2)  =-9.81*(m2*b2);
    G(2,1)  =-9.81*(m2*b2);
    G(2,2)  =-9.81*(m2*b2);
end

function x = piNpi(x)

    if x <-pi
        x = x + 2*pi;
    end
    if x > pi
        x = x - 2*pi;
    end
    
end