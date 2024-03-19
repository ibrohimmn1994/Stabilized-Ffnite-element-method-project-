
classdef ANM_FEM_Part1
    methods(Static)
        function CN = GFEM(hmax, command, plotting)
            CFL = 0.5;
            T = 1;

            geometry = @circleg ;
            [p ,e , t ] = initmesh( geometry, 'hmax' , hmax );
            if command == "continuous"
                u = Functions_part1.initial_u_continuous(p(1,:),p(2,:))';
            else
                u = Functions_part1.initial_u_discrete(p(1,:),p(2,:))';
            end

            ui = u;

            M = Functions_part1.Mass(p,t);
            C = Functions_part1.Convection(p,t);

            I = eye( length ( p ));

            % compute the time-stepping
            F_prime_norm = Functions_part1.L_inf_norm(p);
            k = CFL*hmax/F_prime_norm;

            timing = 0;
            m = round(T/k);
            for i = 1:m   
   
                timing = timing+ k;
    
                %display(timing)
    
                A = ( M/k+C/2 );
                b = ( (M/k-C/2) * u );
                A( e(1 ,:) ,:) = I( e(1 ,:) ,:);
                b ( e (1 ,:)) = 0; 
                u = A\b;
   
                if plotting == 1
                    figure(1)
                    pdeplot(p,e,t,'XYData',u,'ZData',u)
                    xlabel('x')
                    ylabel('y')
                    zlabel('u')
                    set(gca,'FontSize',10)
                    %file_name = 'ANM_imagesd16\image' + string(i) + '.png';
                    %saveas(gcf,file_name)
                end
        
            end
            e = ui - u;
            error = sqrt(e'*M*e);
            CN = error;
            fprintf('For hmax = %0.4f the error =  %f \n\n\n',hmax, error);

        end

        function P = plot_errors(h, command)
            Error = zeros(1,length(h));
            for i = 1:length(h)
                Error(i) = ANM_FEM_Part1.GFEM(h(i), command, 0);
            end

            alpha = polyfit(log(h),log(Error),1);
            
            %ch = polyval(alpha,h);
            figure(2)
            xticks([1/32, 1/16, 1/8, 1/4])
            loglog(h, Error,'r','LineWidth',2)
            hold on
            loglog(h,h.^alpha(1),'b','LineWidth',2)

            xlabel('Mesh size h','FontSize',15)
            ylabel('Error','FontSize',15)
            %xlim([0 1/2])
            legend( {'Error norm','h^{\alpha}'} ,'location', 'SouthEast')
            set(gca,'FontSize',20)   

        end
    end
end
%#####################################################################
