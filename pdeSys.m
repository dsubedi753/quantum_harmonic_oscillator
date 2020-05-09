 classdef pdeSys < handle
    % pdeSys Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess = public, SetAccess = private, Hidden = true)
        X_all               % Space Vector plus end noddes
        M                   % The 2nd Derivative matrix
        h                   % Internode Length
        t                   % time
        b_cond = "drc"      % boundary condition type (nul)
        Ua = 0              % End Node 1
        Ub = 0              % End Node 2
        source              % Source Vector
        K = 1i              % K value
        Uf                  %Final Value
    end
    
    properties(Access = public, Hidden = true)
        U                   %Time Stepped Vector
        P                   %Probablity
        Up                  %In momentum-space
    end
    
    properties(Access = public)
        pot_fun =   @(z) zeros(size(z))
        init_cond = @(z) exp(1i*2*pi*z)
        hbar = 1.054571817e-34;
    end
    
    properties(GetAccess = public, SetAccess = private)
        N                   %Number of inner nodes
        a                   %Min end point
        b                   %Max end poiny
        X                   %Space Vector
        pX                  %Momentum Vector
        U0                  %Initial Values
    end
    
    %Maker Methods
    methods
        function obj = pdeSys(a_in,b_in,NumInNodex_in)
            %Contructor method
            
            %Correction for inputs
            narginchk(0,3);
            switch nargin
                case 0
                    a_in = 0;
                    b_in = 1;
                    NumInNodex_in = 100;
                case 1
                    b_in = a_in+1;
                    NumInNodex_in = 100;
                case 2
                    NumInNodex_in = 100;
            end
            
             if (b_in < a_in)
                b_in = a_in + b_in;
             end
            
            % Initialization of the system values
            obj.a = a_in;
            obj.b = b_in;
            obj.N = NumInNodex_in;
            
            obj.Ua = obj.init_cond(obj.a);
            obj.Ub = obj.init_cond(obj.b);
            
            %Make the system
            obj.make_innerNodes()
            obj.make_FDmatrix()
        end
        
        function make_FDmatrix(obj)
            %Makes the finite difference Matrix
            FDM = diag(-2*ones(obj.N,1),0) + diag(ones(obj.N-1,1),1) + diag(ones(obj.N-1,1),-1); 
            obj.M = FDM;
        end
        
        function make_innerNodes(obj)
            obj.X_all = linspace(obj.a,obj.b,obj.N+2).';
            obj.X = obj.X_all(2:end-1);
            obj.h = obj.X(2) - obj.X(1);
            
            obj.U0 = obj.init_cond(obj.X);
            obj.source = obj.pot_fun(obj.X);
        end
        
    end
    
    %Calulating Methods
    methods
        
        function calc_final(obj)
            obj.make_FDmatrix();
            obj.Uf = -obj.M\(obj.h^2*obj.source/obj.K);
        end
        
        function calc_all(obj,end_time)
            obj.make_FDmatrix();
             opts = odeset('Reltol',1e-5,'AbsTol',1e-6, 'Stats','on');
            [t,obj.U] = ode113(@dU,[0,end_time],obj.U0.',opts); % time stepping
            function output = dU(t,U)
                output = obj.K/(obj.h^2)*(obj.M*U) + obj.source;
            end
        end
        
        function calc_prob(obj)
            obj.P = obj.U.*conj(obj.U);
        end
        
        function calc_momentum(obj)
            y = fft(obj.U.');                   % Compute DFT of U
            y = y.';
            m = abs(y);                         % Magnitude
            y(m<1e-6) = 0;                      % Remove Bits 
            obj.pX = obj.hbar*2*pi*(-obj.N/2:obj.N/2-1)'; 
            obj.Up = obj.h/(sqrt(obj.hbar*2*pi)).*fftshift(y);
        end
    end
    
    %Plot Methods
    methods
		function plot_pdf_step(obj)
            figure(1)
            clf
            subplot(2,2,1)
            plot(obj.X,obj.P(1,:).' )
            ylim([min(obj.P, [],'all'), max(obj.P, [],'all')])
            xlabel("x")
            ylabel("P(x)")
            grid on
            
            subplot(2,2,2)
            plot(obj.X, obj.P(round(size(obj.P,1)/3),:).')
            ylim([min(obj.P, [],'all'), max(obj.P, [],'all')])
            xlabel("x")
            ylabel("P(x)")
            grid on

            subplot(2,2,3)
            plot(obj.X, obj.P(round(2*size(obj.P,1)/3),:).')
            ylim([min(obj.P, [],'all'), max(obj.P, [],'all')])
            xlabel("x")
            ylabel("P(x)")
            grid on
            
            subplot(2,2,4)
            plot(obj.X, obj.P(end,:).')
            ylim([min(obj.P, [],'all'), max(obj.P, [],'all')])
            xlabel("x")
            ylabel("P(x)")
            grid on
        end
        
        function plot_frame(obj,frame)
            figure(2)
            clf
            
            [X1,Z1] = meshgrid(-2:2,-2:2);
            Y1 = zeros(size(X1));    
            sxz = surf(X1,Y1,Z1);
            set(sxz, 'FaceColor',[0.1 1 0.5],'FaceAlpha',0.1,'EdgeColor','none')
            hold on
            [X1,Y1] = meshgrid(-2:2,-2:2);
            Z1 = zeros(size(X1));    
            sxy = surf(X1,Y1,Z1);
            set(sxy, 'FaceColor',[0.1 0.5 1],'FaceAlpha',0.1,'EdgeColor','none')
            
            plot3(obj.X, imag(obj.U(frame,:).'),real(obj.U(frame,:).'));
            xlabel("x");
            ylabel("\Psi_i");
            zlabel("\Psi_r");
            
            yline(0)
            grid on
            
            xlim([obj.a obj.b])
            ylim([min(min(imag(obj.U))) max(max(imag(obj.U)))])
            zlim([min(min(real(obj.U))) max(max(real(obj.U)))])
        end
		
		function plot_pdf_frame(obj,frame)
            figure(3)
            clf
            plot(obj.X, obj.P(frame,:));
            xlabel("x");
            ylabel("P(x)");
            
            yline(0)
            grid on
            
            xlim([obj.a obj.b])
            ylim([0 max(max(obj.P))])
        end
        
        function plot_continous(obj, speed, last_frame, name)
            %speed -- speed of plotting
            %last_frame -- last frame to be plotted
            %name -- name of the file if movie is to be created
            
            narginchk(1,4)
            switch nargin
                case 1
                    speed = 5;
                    last_frame = 10000;
                case 2 
                    last_frame = 10000;
            end

        % First Plot
            plth = figure(4);
            clf
            axh = axes(plth);
            Uplot_r = real(obj.U(1,:).');
            Uplot_i = imag(obj.U(1,:).');
            Xplot = obj.X;
            
            [X1,Z1] = meshgrid(-2:2,-2:2);
            Y1 = zeros(size(X1));    
            sxz = surf(axh, X1,Y1,Z1);
            set(sxz, 'FaceColor',[0.1 1 0.5],'FaceAlpha',0.1,'EdgeColor','none')
            hold on
            [X1,Y1] = meshgrid(-2:2,-2:2);
            Z1 = zeros(size(X1));    
            sxy = surf(axh,X1,Y1,Z1);
            set(sxy, 'FaceColor',[0.1 0.5 1],'FaceAlpha',0.1,'EdgeColor','none')
            
            lnh = plot3(axh, Xplot, Uplot_i, Uplot_r);

            xlabel("x");
            ylabel("\Psi_i");
            zlabel("\Psi_r");
            
            yline(0)
            grid on
            
            xlim([obj.a obj.b])
            ylim([min(min(imag(obj.U))) max(max(imag(obj.U)))])
            zlim([min(min(real(obj.U))) max(max(real(obj.U)))])
            
            
            lnh.XDataSource = 'Xplot';
            lnh.YDataSource = 'Uplot_i';
            lnh.ZDataSource = 'Uplot_r';
            
            if nargin > 3
                %Write a Movie
                count = 1;
                movieVector(count) = getframe(gcf);
                
                for ii = 2:speed:min(last_frame,size(obj.U,1))
                    Uplot_i = imag(obj.U(ii,:).');
                    Uplot_r = real(obj.U(ii,:).');
                    refreshdata(lnh, 'caller')
                    drawnow
                movieVector(count) = getframe(gcf);
                count = count +1;
                end
                
                % Movie Writer
                movieObj = VideoWriter(name);
                movieObj.FrameRate = 30;
                open(movieObj);
                writeVideo(movieObj, movieVector);
            else
                % Continue Plotting
                for ii = 2:speed:min(last_frame,size(obj.U,1))
                    Uplot_i = imag(obj.U(ii,:).');
                    Uplot_r = real(obj.U(ii,:).');
                    refreshdata(lnh, 'caller')
                    drawnow
                end
            end
        end
        
        function plot_pdf_continous(obj,speed, last_frame, name)
            %speed -- speed of plotting
            %last_frame -- last frame to be plotted
            %name -- name of the file if movie is to be created
            
            narginchk(1,4)
            switch nargin
                case 1
                    speed = 20;
                    last_frame = 100000;
                case 2 
                    last_frame = 100000;
            end

        % First Plot
			plth = figure(5);
            clf
			axh = axes(plth);
			Pplot = obj.P(1,:).';
			Xplot = obj.X;
			lnh = plot(axh, Xplot, Pplot);


            xlabel("x");
			ylabel("P(x)");
			xlim([obj.a obj.b])
			ylim([0 max(max(obj.P))])

            
            yline(0)
            grid on
            
            
			lnh.XDataSource = 'Xplot';
			lnh.YDataSource = 'Pplot';
            
            if nargin > 3
                %Write a Movie
                count = 1;
                movieVector(count) = getframe(gcf);
                
                for ii = 2:speed:min(last_frame,size(obj.U,1))
					Pplot = obj.P(ii,:).';
					refreshdata(lnh, 'caller')
					drawnow
					movieVector(count) = getframe(gcf);
					count = count +1;
                end
                
                % Movie Writer
                movieObj = VideoWriter(name);
                movieObj.FrameRate = 30;
                open(movieObj);
                writeVideo(movieObj, movieVector);
            else
                % Continue Plotting
                for ii = 2:speed:min(last_frame,size(obj.U,1))
                    Pplot = obj.P(ii,:).';
					refreshdata(lnh, 'caller')
					drawnow
                end
            end
        end
        
        function plot_momentum(obj,frame)
            figure(6)
            clf
            narginchk(1,2)
            if nargin == 1
                frame = 1;
            end
            plot(obj.pX,abs(obj.Up(frame,:)));
            xlabel("Momentum")
            ylabel("\psi(p)")
            grid on
        end
    end
    
    %Set Methods
    methods
        function set.b_cond(obj, cond)
            %Change boundary condition
            if cond == "drc"
                disp('Dirichlet Condition is Set')
                flag = 1;
            elseif cond == "vnm"
                disp('Von-Neumann Condition is Set')
                flag = 1;
            else
                disp('Condition is not valid')
                flag = 0;
            end
            
            if flag == 1  
                obj.b_cond = cond;
                obj.make_FDmatrix()
            end
        end
        
        function set.init_cond(obj,sf_in)
            %Change Initial Condition
            obj.init_cond = sf_in;
            obj.U0 = obj.init_cond(obj.X);
            norm = 1/sqrt(trapz(obj.X, obj.U0.*conj(obj.U0)));
            obj.U0 = obj.U0.*norm;
            obj.Ua = 0;
            obj.Ub = 0;
            % Change the matrices
            obj.make_FDmatrix()
            disp('The initial function has changed')
        end
        
        function set.pot_fun(obj,sf_in)
            %Change source
            obj.pot_fun = sf_in;
            obj.source = obj.pot_fun(obj.X);
        end
    end
    
end

