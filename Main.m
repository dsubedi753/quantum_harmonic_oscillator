% clear, clc
% sys = pdeSys;
option = 1;
failed = true;
arr = @(f) trapz(sys.X, sys.U(f,:).*conj(sys.U(f,:)));
while ~(option == -1)
    fprintf("1- Input a new inital function (Make a new system) \n");
    fprintf("2- Plot a specific frame of wavefunction \n");
    fprintf("3- Continously plot the wavefunction \n");
    fprintf("4- Overview of evolution of Probability Density Function \n");
    fprintf("5- Plot a specific frame of Probability Density Function \n");
    fprintf("6- Continously plot the Probability Density Function \n");
    fprintf("7- Plot wavefunction in momentum domain \n");
    fprintf("-1- Exit \n");
    option = input('Enter an option: ');
    switch option
        case 1
            fprintf('Enter a function handle for the intial contion \n')
            fprintf('For Example @(z) sin(pi*z) sin(2*pi*z) \n');
            fprintf('Make sure the baounday value is zero \n');
            while failed
                sf = input('::>> ');
                if abs(sf(0)) < 1e-10 && abs(sf(1)) < 1e-10 && ~all(sf(linspace(0,1)) == zeros(1,100))
                    failed = false;
                    fprintf("The system is being simulated. This might take a while \n")
                    sys.init_cond = sf ;
                    sys.calc_all(2);
                    fprintf("Calculating Probility Density Function \n")
                    sys.calc_prob;
                    fprintf("Transforming into momentum space \n")
                    sys.calc_momentum;
                else
                    fprintf("The function is not 0 at boundary or is a zero function\n");
                end
            end
            
        case 2
            frame = input("Enter Frame: ");
            sys.plot_frame(frame)
        case 3
            sys.plot_continous(10);
        case 4
            row = input("Enter number of rows");
            column = input("Enter number of column");
            sys.plot_pdf_step(row,column)
        case 5
            frame = input("Enter Frame: ");
            sys.plot_pdf_frame(frame)
        case 6
            sys.plot_pdf_continous(10);
        case 7
            sys.plot_momentum
    end
end