

classdef ANM_FEM_Part2
    methods(Static)
        function CN = Solve(hmax, initial_condition, method, plotting)
            CFL = 0.5;
            T = 1;
            gamma = hmax/2;
            geometry = @circleg ;
            [p ,e , t ] = initmesh( geometry, 'hmax' , hmax );

            if initial_condition == "continuous"
                u = Functions_part2.initial_u_continuous(p(1,:),p(2,:))';
            elseif initial_condition == "discontinuous"
                u = Functions_part2.initial_u_discrete(p(1,:),p(2,:))';
            else
                return
            end

            ui = u;
            M = Functions_part2.Mass(p,t);
            C = Functions_part2.Convection(p,t);
            I = eye( length ( p ));
            % compute the time-stepping
            beta_inf_norm = Functions_part2.Beta_inf_norm(p);
            k = CFL*hmax/beta_inf_norm;
            timing = 0;
            m = round(T/k);
            %m_eps = 0;

            if method == "RV"
                Res = zeros(length(p));
                for i = 1:m   
        
                    timing = timing+ k;
                    %display(timing)    
                    epsilons = Functions_part2.Epsilon(p,t,Res, hmax);
                    %m_eps = mean(epsilons);
                    S = Functions_part2.Stiffness(p,t, epsilons); 
                    A = ( M/k + C/2 + 1/2*S );
                    b = ( ( M/k - C/2 - 1/2*S) * u );
                    A( e(1 ,:) ,:) = I( e(1 ,:) ,:);
                    b ( e (1 ,:)) = 0; 
                    u_prev = u;
                    u = A\b;
    
                    Res = (1/k*(u - u_prev ) + M\C*u) ; 
                    Res = Res/max(u - mean(u));
                    % max(Res)   
                    if plotting == 1
                        figure(1)
                        pdeplot(p,e,t,'XYData',u,'ZData',u)
                        xlabel('x')
                        ylabel('y')
                        zlabel('u')
                        zlim([-0.3,1.2])
                        set(gca,'FontSize',10)
                    %   file_name = 'RV_d_16\image' + string(i) + '.png';
                    %  saveas(gcf,file_name)
                    end
                    
    
                end
           
     
            elseif method == "SUPG"
    
                S = Functions_part2.SUPG_Stiffness(p,t);
                for i = 1:m   
        
                    timing = timing+ k;
                    %display(timing)   
                    A = ( M/k + C/2 + gamma/k*C'  + gamma/2*S );
                    b = ( ( M/k - C/2 + gamma/k*C' - gamma/2*S) * u );
                    A( e(1 ,:) ,:) = I( e(1 ,:) ,:);
                    b ( e (1 ,:)) = 0; 
                    u = A\b;

                    if plotting == 1
                        figure(1)
                        pdeplot(p,e,t,'XYData',u,'ZData',u)
                        xlabel('x')
                        ylabel('y')
                        zlabel('u')
                        zlim([-0.3,1.2])
                        set(gca,'FontSize',10)
                        % file_name = 'SUPG_d_16\image' + string(i) + '.png';
                        % saveas(gcf,file_name)
                    end    
                end
            else
                return;
            end
            e = ui - u;
            error = sqrt(e'*M*e);
            CN = error;
            
            fprintf('For hmax = %0.4f the error =  %f \n\n\n',hmax, error);

        end
        
        function Error = plot_errors(h, initial_conditions, method)
            Error = zeros(1,length(h));
            %eps = zeros(1,length(h));


            for i = 1:length(h)
                disp(i)
                Error(i) = ANM_FEM_Part2.Solve(h(i), initial_conditions,method, 0);
            end


            alpha = polyfit(log(h),log(Error),1);
            disp(alpha)
            %ch = polyval(alpha,h);
            figure(2)
            xticks([1/32, 1/16, 1/8, 1/4])
            loglog(h, Error,'r','LineWidth',2)
            hold on
            loglog(h,h.^alpha(1),'b','LineWidth',2)
            %loglog(h, eps,'r','LineWidth',2)
            grid on

            xlabel('Mesh size h','FontSize',15)
            ylabel('Error','FontSize',15)
            %xlim([0 1/2])
            %legend( {'mean viscosity'} ,'location', 'SouthEast')

            legend( {'Error norm','h^{\alpha}'} ,'location', 'SouthEast')
            set(gca,'FontSize',20)   

        end
    end
end

%#############################################################################
