classdef ssNN_Layer
   
    properties
        num = 0;
        W   = 1;
        dW  = 0;
        Wb  = 1;
        dWb = 0;
        O   = 0;
        I   = 0;
        E   = 0;
        d   = 0;
        s   = 1;
        var = 0;
        nOut= 0;
        nIn = 0;
        RMSE= 100;
    end
    
    
    methods
        
        function nn = ssNN_Layer(nIn,nOut,sen,var,num)
            nn.num  = num;
            nn.var  = var;
            nn.nIn  = nIn;
            nn.nOut = nOut;
            nn.s    = sen;
            nn.W    = randn(nOut,nIn)*1e-3;
            nn.Wb   = randn(nOut,1  )*1e-3;
            nn.dW   = randn(nOut,nIn)*1e-3;
            nn.dWb  = zeros(nOut,1);
            nn.I    = zeros(nIn,1);
            nn.O    = zeros(nOut,1);
            nn.E    = zeros(nOut,1);
            nn.d    = zeros(nOut,1);
        end
        
               
        function nn = activate(nn,varargin)
            if  strcmp(class(varargin{1}),'ssNN_Layer')
                if nn.nIn == varargin{1}.nOut 
                    nn.I = varargin{1}.O;
                else
                    error('Layers are incompatible for activation.')
                end
            else
                nn.I = varargin{1};
            end
            nn.O    = tanh(nn.s*(nn.W*nn.I+nn.Wb));
            
            if any(isnan(nn.O))
                error('Output is NaN')
            end
        end
        
        
        
        function nn = reweight_form(nn,lr,varargin)
            if  strcmp(class(varargin{1}),'ssNN_Layer')
                nn.E = (varargin{1}.W.')*varargin{1}.d;
            else
                nn.E = nn.O-varargin{1};
            end
            m2      = (lr/2)^2;
            nn.d    = nn.s*(1-nn.O.^2).*nn.E;
            nn.dW   = -lr*nn.d*transpose(nn.I) + m2*(nn.dW + randn(nn.nOut,nn.nIn)*nn.var);
            nn.dWb  = -lr*nn.d + m2*(nn.dWb + randn(1,1)*nn.var);
            nn.RMSE = norm(nn.E);
        end        
        
        
        
        function nn = reweight_form_err(nn,lr,error)
            nn.E 	= -error;
            m2      = (lr/2)^2;
            nn.d    = nn.s*(1-nn.O.^2).*nn.E;
            nn.dW   = -lr*nn.d*transpose(nn.I) + m2*(nn.dW + randn(nn.nOut,nn.nIn)*nn.var);
            nn.dWb  = -lr*nn.d + m2*(nn.dWb + randn(1,1)*nn.var);
            nn.RMSE = norm(nn.E);
        end         
        
		
		
        function nn = reweight_set(nn)
            nn.W = nn.W  + nn.dW;
            nn.Wb= nn.Wb + nn.dWb;
        end
        
        
		
        function plot(nn)
            bar3(nn.W)
        end
		
    end
end