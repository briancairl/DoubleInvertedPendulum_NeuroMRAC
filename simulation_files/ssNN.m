classdef ssNN
    properties
        name        = '';
        layers      = {};
        rmse_data   = {};
        N_layers    = 0;
    end
    
    
    methods
        
        function ann = ssNN(name,layerSpec,sen,var)
            for idx = 1:(numel(layerSpec)-1)
                ann.name            = name;
                ann.rmse_data{idx}  = [];
                ann.layers{idx}     = ssNN_Layer(layerSpec(idx),layerSpec(idx+1),sen,var,idx);
                ann.N_layers        = ann.N_layers + 1;
                
                %%%%%%%%%%%%%%%%%%%%
                % Load Old Weights %
                %%%%%%%%%%%%%%%%%%%%
                try 
                    w0 = load([name,'_W_',num2str(idx)]);
                    [wn,wm]= size(w0.w0);
                    if  all([wm,wn] == layerSpec(idx:idx+1))
                        ann.layers{idx}.W = w0.w0; 
                    else
                        warning([ann.name,'_W_',num2str(idx),' exists, but was not dimensionally conistent.'])
                    end
                end 
            end
        end
        
        
        function save(ann)
            for idx = 1:ann.N_layers
                w0 = ann.layers{idx}.W; %#ok<NASGU>
                save([ann.name,'_W_',num2str(idx)],'w0');
            end
        end
        
        
        function ann = activate(ann,input)
            ann.layers{1} = ann.layers{1}.activate(input);
            for idx = 2:ann.N_layers
                ann.layers{idx} = ann.layers{idx}.activate(ann.layers{idx-1});
            end
        end
        
        
        function out = output(ann)
            out = ann.layers{end}.O;
        end
        
        
        function ann = reweight(ann,lr,target)
            ann.layers{ann.N_layers} = ann.layers{ann.N_layers}.reweight_form(lr,target);
            for idx = (ann.N_layers-1):-1:1
                ann.layers{idx} = ann.layers{idx}.reweight_form(lr,ann.layers{idx+1});
            end
            for idx = 1:ann.N_layers
                ann.layers{idx}          = ann.layers{idx}.reweight_set();
                ann.rmse_data{idx}(end+1)= ann.layers{idx}.RMSE;
            end
        end
        
        function ann = reweight_e(ann,lr,error)
            ann.layers{ann.N_layers} = ann.layers{ann.N_layers}.reweight_form_err(lr,error);
            for idx = (ann.N_layers-1):-1:1
                ann.layers{idx} = ann.layers{idx}.reweight_form(lr,ann.layers{idx+1});
            end
            for idx = 1:ann.N_layers
                ann.layers{idx}          = ann.layers{idx}.reweight_set();
                ann.rmse_data{idx}(end+1)= ann.layers{idx}.RMSE;
            end
        end
        
        
        function blackbox(ann,smoothing,gain)
            for idx = 1:ann.N_layers
                subplot(ann.N_layers,1,idx)
                plot(smooth((gain*ann.rmse_data{idx}(1,:)),smoothing));
                ylabel( ['Layer ',num2str(idx)','RMSE'] )
            end
        end
        
        
        function plot_weights(ann)
            for idx = 1:ann.N_layers
                subplot(ann.N_layers,1,idx)
                bar3(ann.layers{idx}.W);
                ylabel( ['Layer ',num2str(idx)','RMSE'] )
            end
        end
        
    end
end