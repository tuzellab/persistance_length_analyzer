%% Fourier analysis with or without decimation
%  Last modified 10/24/2021
%  Written by Dr. Pattipong Wisanpitayakorn
%  Tuzel Lab, Bioengineering Department, Temple University
%  For questions, please email: erkan.tuzel@temple.edu

clear all;
close all;

% In the first example, the files are located at 'Coordinates_for_PLA/Monte_Carlo_Lp3000_pixelsize100nm'
% Use the MATLAB "Current Folder" window to navigate into that folder before running this code.
% The filename format follows: file0001.dat, file0002.dat, etc.
% File numbering starts at 1 and ends at 1000 (Vstart = 1 and Vend = 1000)
% In this sample data set, persistence length is 3 mm (or 3,000 micron).
% Use fit error = 0 in this example, as the data does not have any experimental error

% In the second example, we use the data set located at 'Coordinates_for_PLA/Tuzel_Lab_MT_GaussianTracing_Pixelsize107p5'
% Use the MATLAB "Current Folder" window to navigate into that folder before running this code.
% This dataset contains real experimental data with its associated tracking
% errors. While everything else is the same as the first example, fit_error
% should be set to 1 to demonstrate accuracy of using the modified Fourier
% formula. Number of data points in this data set is 198.


% Bootstrapping.
% To use bootstrapping, change BS from 0 to 1. 
% The code will randomly choose a number of filament equal to BS_size from
% a range of file specify by Vstart and Vend (for example, if BS_size = 50,
% Vstart = 200, and Vend = 400, the code will randomly choose 50 files
% between file0200.dat to file0400.dat for the analysis. The process will
% be repeated for a number of times equal to N_BS. Each iteration will
% produce Lp. The reported Lp at the end is the mean of Lp from all
% iterations.

% Assign range of files to read coordinates from
filenameprefix=sprintf('file');
Vstart = 1;
Vend = 1000;

% Input Parameters
pixelsize = 100; % pixel size in nm
H = 20; % number of Fourier modes to analyze
fit_error = 1; % to fit error function or not

% Bootstrapping parameters
BS = 1; % 1 means performing bootstrapping. 0 means no bootstrapping
N_BS = 50; % number of bootstrapping iterations
BS_size = 200; % size of each bootstrapping bins (number of filaments per iteration)

V = Vend - Vstart +1; %number of filaments

if (BS == 1)
    count = 0;
    fprintf('Progress Percentage (only displays properly when N_BS > 10) \n')
    
    for vvv = 1:N_BS
        
        if(mod(vvv, N_BS/10) == 0)
            count = count+10;
            fprintf('%d %% \n', count) % print progress at every 10%
            % Note that this only display properly when N_BS > 10
        end
        
        clear num_filament filament V
        
        num_filament = Vend-Vstart+1; %total number of filament
        
        % Generate an array of randomized indicies for selected filaments
        % for bootstrapping process
        filament = Vstart+ceil(num_filament*rand(BS_size,1))-1;
        
        for v = 1:BS_size;
            clear theta f_original A B x y xx yy ds an a0 sk;
            
            filename = sprintf('./%s%04d.dat', filenameprefix, filament(v));
            
            % read and store x-y coordinate from the files
            A = load(filename);
            B = size(A);
            for j = 1:B(1)
                x(j) = A(j,1);
                y(j) = A(j,2);
            end
            
            for j = 1:B(1)-1
                dx(j) = x(j+1)-x(j);
                dy(j) = y(j+1)-y(j);
                ds(j) = sqrt((dx(j))^2 + (dy(j))^2); % calculate the distance between bonds
                f_original(j) = atan2(dy(j),dx(j)); % calculate the angle of each bond
                sk(j) = sum(ds) - 0.5*ds(j); % midpoint calculation (sk)
            end
            
            L = sum(ds); %Calculate Length of filament
            storeL(v) = L;
            
            %Fourier calculation
            a0 = 0;
            an=zeros(1,H);
            bn=zeros(1,H);
            for j = 1:B(1)-1
                a0=a0+sqrt(2/L)*f_original(j)*ds(j);
                for n = 1:H
                    an(n) = an(n)+sqrt(2/L)*f_original(j)*cos(n*pi*sk(j)/(L))*ds(j);
                end
                
                f_fourier(j)= sqrt(2/L)*a0/2;
                for n=1:H
                    f_fourier(j) = f_fourier(j) + sqrt(2/L)*an(n)*cos(n*pi*sk(j)/(L));
                end
            end
            
            for i = 1:H
                An(v,i) = an(i)/L; % store a_n values for each mode of each filament
            end
            
        end
        
        mean_an = sum(An)/length(An); % find mean of a_n for each mode
        var_ansq = sum(An.*An)/length(An); % find variance
        variance = var_ansq - mean_an.*mean_an;  % subtract the mean
        
        n = 1:H; %Assign an array with values from 1 to H
        
        % guess values for the fit
        bguess(1)=sqrt(1./(pi^2*n(1)^2.*variance(1)));
        if(fit_error==1)
            bguess(2) = 0;
            bguess(3) = 0;
        else
        end
        
        % specify fit function for Fourier analysis
        if(fit_error == 1)
            Fit_Fourier = @(bguess,n)1./(pi^2*bguess(1)^2*n.^2) + bguess(2)^2*n.^2 +bguess(3)^2;
        else
            Fit_Fourier = @(bguess,n)1./(pi^2*bguess(1)^2*n.^2);
        end
        
        bfit=nlinfit(n,variance,Fit_Fourier,bguess); %Fit
        Lp = bfit(1)^2; %Measured Lp is from the first parameter of the fit
        Store_Lp(vvv) = Lp;
        
        
        % Reconstruct the variance using the values from the fit
        if(fit_error==1)
            variance_fit = 1./(pi^2*bfit(1)^2*n.^2) + bfit(2)^2*n.^2 +bfit(3)^2;
            b2 = bfit(2)^2;
            b3 = bfit(3)^2;
        else
            variance_fit = 1./(pi^2*bfit(1)^2*n.^2);
        end
        
    end
    
    Result_Lp = mean(Store_Lp);
    std_Result_Lp = std(Store_Lp);
    
    Result_Lp = Result_Lp*pixelsize*1E-9*1E6; % Convert Lp units from pixel to micron
    std_Result_Lp = std_Result_Lp*pixelsize*1E-9*1E6; % Convert Lp units from pixel to micron
    fprintf('Lp = %4.f %s %4.f \x3BCm \n', Result_Lp, char(177), std_Result_Lp) % print Lp on screen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif (BS == 0)
    for v = 1:V
        clear theta f_original A B x y xx yy ds an a0 sk;
        
        filename = sprintf('./%s%04d.dat',filenameprefix,Vstart+v-1);
        
        % read and store x-y coordinate from the files
        A = load(filename);
        B = size(A);
        for j = 1:B(1)
            x(j) = A(j,1);
            y(j) = A(j,2);
        end
        
        for j = 1:B(1)-1
            dx(j) = x(j+1)-x(j);
            dy(j) = y(j+1)-y(j);
            ds(j) = sqrt((dx(j))^2 + (dy(j))^2); % calculate the distance between bonds
            f_original(j) = atan2(dy(j),dx(j)); % calculate the angle of each bond
            sk(j) = sum(ds) - 0.5*ds(j); % midpoint calculation (sk)
        end
        
        L = sum(ds); %Calculate Length of filament
        storeL(v) = L;
        
        %Fourier calculation
        a0 = 0;
        an=zeros(1,H);
        bn=zeros(1,H);
        for j = 1:B(1)-1
            a0=a0+sqrt(2/L)*f_original(j)*ds(j);
            for n = 1:H
                an(n) = an(n)+sqrt(2/L)*f_original(j)*cos(n*pi*sk(j)/(L))*ds(j);
            end
            
            f_fourier(j)= sqrt(2/L)*a0/2;
            for n=1:H
                f_fourier(j) = f_fourier(j) + sqrt(2/L)*an(n)*cos(n*pi*sk(j)/(L));
            end
        end
        
        for i = 1:H
            An(v,i) = an(i)/L; % store a_n values for each mode of each filament
        end
    end
    
    mean_an = sum(An)/length(An); % find mean of a_n for each mode
    var_ansq = sum(An.*An)/length(An); % find variance
    variance = var_ansq - mean_an.*mean_an;  % subtract the mean
    
    n = 1:H; %Assign an array with values from 1 to H
    
    % guess values for the fit
    bguess(1)=sqrt(1./(pi^2*n(1)^2.*variance(1)));
    if(fit_error==1)
        bguess(2) = 0;
        bguess(3) = 0;
    else
    end
    
    % specify fit function for Fourier analysis
    if(fit_error == 1)
        Fit_Fourier = @(bguess,n)1./(pi^2*bguess(1)^2*n.^2) + bguess(2)^2*n.^2 +bguess(3)^2;
    else
        Fit_Fourier = @(bguess,n)1./(pi^2*bguess(1)^2*n.^2);
    end
    
    bfit=nlinfit(n,variance,Fit_Fourier,bguess); %Fit
    Result_Lp = bfit(1)^2; %Measured Lp is from the first parameter of the fit
    Result_Lp = Result_Lp*pixelsize*1E-9*1E6; % Convert Lp units from pixel to micron
    
    
    % Reconstruct the variance using the values from the fit
    if(fit_error==1)
        variance_fit = 1./(pi^2*bfit(1)^2*n.^2) + bfit(2)^2*n.^2 +bfit(3)^2;
        b2 = bfit(2)^2;
        b3 = bfit(3)^2;
    else
        variance_fit = 1./(pi^2*bfit(1)^2*n.^2);
    end
    
    fprintf('Lp = %4.f \x3BCm \n', Result_Lp) % print Lp on screen
end

    
    
    % Two options of plot: normal plot, or plot in log scale.
    variance = variance/(pixelsize*1E-9); % convert units from 1/pixel to 1/m
    variance_fit = variance_fit/(pixelsize*1E-9); % convert units from 1/pixel to 1/m
    plot(n,variance,'r',n,variance_fit,'b-')
    xlabel('mode number (n)');
    ylabel('<((a_n-a_0)/L)^2> (1/m)');
    figure
    loglog(n,variance,'r',n,variance_fit,'b-')
    xlabel('mode number (n)');
    ylabel('<((a_n-a_0)/L)^2> (1/m)');