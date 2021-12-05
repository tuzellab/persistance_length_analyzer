%% Curvature analysis with bootstrapping
%  Last modified 11/24/2021
%  Written by Dr. Pattipong Wisanpitayakorn
%  Tuzel Lab,  Bioengineering Department, Temple University
%  For questions, please email: erkan.tuzel@temple.edu

clear all;
close all;

%%% Input Parameters

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

filenameprefix=sprintf('file');
Vstart = 1; % starting file index (1 means file0001.dat)
Vend = 1000; % ending file index
pixelsize = 100; % pixel size in nm

% Bootstrapping parameters
BS = 0; % 1 means performing bootstrapping. 0 means no bootstrapping
N_BS = 10; % number of bootstrapping iterations
BS_size = 50; % size of each bootstrapping bins (number of filaments per iteration)

% From which m (skipping) to which m to measure the curvature of those filaments
% m = 1 >> no skipping
% m = 40 >> skip 40 nodes for each measurement
mstart =1; % first m
mend = 30; % last m

mu_eq_1p5 = 1; % 1 means use fixed µ of 1.5 for Lambda equation

y_int = 0; % 1 means allow y-intercept

%%%

num_filament = Vend-Vstart+1; %total number of filament

count = 0;

if (BS == 1)
    
    for vvv = 1:N_BS % loop over all sets
        fprintf('Progress Percentage (only displays properly when N_BS > 10) \n')
        if(mod(vvv, N_BS/10) == 0)
            count = count+10;
            fprintf('%d %% \n', count) % print progress at every 10%
            % Note that this only display properly when N_BS > 10
        end
        
        clearvars -except pixelsize foldername filenameprefix y_int mu_eq_1p5 Store_Lpfit mstart mend binsize bguess Lp_guess vvv Vstart Vend N_BS BS_size num_filament count;
        close all;
        
        filament = Vstart+ceil(num_filament*rand(BS_size,1))-1;
        
        %   Count m
        countm = 0;
        
        %   Loop over curvature calculation for all skipping indices (m)
        for m = mstart:1:mend  %skipping over m-1 nodes
            
            %       At the start of each loop, deallocate all matrices except the followings
            clearvars -except pixelsize foldername filenameprefix y_int mu_eq_1p5 Store_Lpfit filament N_BS BS_size num_filament vvv  mu_eq_1p5 L y_int h stats bguess mstart mend Vstart Vend Nstart Nend binsize m Numbern mu Lp_guess countm store_all_acdss store_all_bfit write_bfitvsmeands write_curvdist m_write_curvdist exp_mu count;
            
            % Declare L_av for average filament length calculation
            L_av=0;
            
            % kk is an index for storing curvature
            kk=1;
            
            % sizedx is an index for storing the bond spacing on the filament
            sizeds =1;
            
            %       Loop over specified filaments
            for v = 1:BS_size
                
                %           Clear all necessary matrices
                clear theta x y xx yy ds C D E F xskip yskip store state storej;
                
                %           Write coordinates file name to a matrix
                filename = sprintf('./%s%04d.dat', filenameprefix, filament(v));
                
                %           Load coordinates from a file into array A
                A = load(filename);
                
                %           B stores the size of the array A >> example output B = [95 2]
                B = size(A);
                
                %           M stores the number of coordinate points (number of bonds)
                %           pick the first column from B
                M = B(1);
                
                %           Store x and y coordinates from the files
                for j = 1:M
                    x(j) = A(j,1);
                    y(j) = A(j,2);
                end
                
                %           Calculate bond spacing in x- and y-projection
                
                for j = 1:M-1
                    xx(j) = x(j+1)-x(j); % displacement in x
                    yy(j) = y(j+1)-y(j); % displacement in y
                end
                
                %           Calculate the bond spacing using pythagoras
                for j = 1:M-1
                    ds(j) = sqrt((xx(j))^2 + (yy(j))^2);
                end
                
                %           Calculate length of the filament
                L(v) = sum(ds);
                %-------------------------------------------
                
                %           Loop for doing recursion
                %           Recursion is the proccess to get more data out of a filament
                %           If we skipped some nodes and measured the curvature, the
                %           information on those skipped nodes would be wasted
                %           Therefore, after measuring the curvature for the first sets of
                %           coordinates on the filament, starting from the first point, we
                %           now start skipping and measuring the curvature from the second
                %           point on the chian (This point was previously wasted during the
                %           last calculation.
                
                for q = 1:m % recursion  % If no-recursion, use "for q = 1:1"
                    
                    %             Deallocate necessary array
                    clear C D dss E F
                    
                    %               Loop to find out how many big bonds (after skipping) the
                    %               filament can have
                    b = floor ((M-q)/m);
                    
                    % Loop to find the bond spacing of the big bonds
                    for i = 1:b
                        C(i) = x(q+(m)*i)-x(q+(m)*(i-1));  %displacement in x
                        D(i) = y(q+(m)*i)-y(q+(m)*(i-1));  %displacement in y
                        dss(i) = sqrt(C(i).^2 + D(i).^2);  %bond spacing
                        acdss(sizeds) = dss(i); % store bond spacing
                        sizeds = sizeds+1; % ready to store the next one
                    end
                    
                    % Loop to calculate curvature
                    for i = 1:b-1
                        
                        if (C(i+1)>C(i)) && (D(i+1)>D(i)) %if pointing in 1st quadrant
                            dtheta_all(i) = acos((C(i)*C(i+1) + D(i)*D(i+1))/(dss(i)*dss(i+1))); %calculate angle
                        elseif (C(i+1)<C(i)) && (D(i+1)<D(i)) %3rd quadrant
                            dtheta_all(i) = acos((C(i)*C(i+1) + D(i)*D(i+1))/(dss(i)*dss(i+1))); %calculate angle
                        else %2nd and 4th quadrant
                            dtheta_all(i) = -acos((C(i)*C(i+1) + D(i)*D(i+1))/(dss(i)*dss(i+1))); %calculate angle
                        end
                        
                        %Calculate and store curvature from the angle and bond
                        %length
                        curvature(kk)= dtheta_all(i)/((dss(i)+dss(i+1))*0.5);
                        kk=kk+1; %ready to store the next one
                    end
                end % End recursion loop
            end % End loop over filament
            
            % After storing all curvature values for current m, subtract the mean and take
            % the absolute value
            curvature = real(curvature);
            %curvature = abs(curvature-mean(curvature));
            
            [h, stats] = cdfplot(curvature);
            hold on;
            [cdf_y,cdf_x] = ecdf(curvature);
            bfit_guess = [getfield(stats, 'mean') getfield(stats, 'std')];
            cumulative_gaussian = @(bfit,x)0.5*(1+erf((x-bfit(1))/(sqrt(2)*bfit(2))));
            bfit=nlinfit(cdf_x,cdf_y,cumulative_gaussian,bfit_guess);
            cdf_y_fit2 =  normcdf(cdf_x,bfit(1),bfit(2));
            plot(cdf_x,cdf_y_fit2,'r');
            hold off;
            
            store_all_bfit(m-mstart+1) = 1/(bfit(2))^2;
            store_all_acdss(m-mstart+1) = mean(acdss);
            
            %To write measured and fitted curvature distribution functions
        end
        
        % Calculate average length of all filaments
        L_av = sum(L)/length(L);
        
        % Assign a matrix containing all m values
        Numberm=mstart:mend;
        %calculate mu for all m
        
        if(mu_eq_1p5)
            mu_theory = 1.5;
        else
            mu_theory=Numberm./(((2.*Numberm-1).*(1-1./Numberm)/3)+1);
        end
        
        % Estimate gauss persistence length (guess Lp)
        Lp_guess = store_all_bfit(end)/(mu_theory(end)*store_all_acdss(end));
        
        if(y_int)
            Fit_Lp_Function = @(c,x)c(1).*x+c(2);
            Lp_guess = [Lp_guess,0];
            fit=nlinfit(mu_theory.*store_all_acdss,store_all_bfit,Fit_Lp_Function,Lp_guess);
            fit_lambda = fit(1).*mu_theory.*store_all_acdss+fit(2);
            Lp = fit(1);
            y_intercept = fit(2);
            figure;  % plot CDF of curvature distribution
            plot(mu_theory.*store_all_acdss,store_all_bfit,'or',mu_theory.*store_all_acdss,fit_lambda); set(gca,'FontSize',15)
        else
            Fit_Lp_Function = @(c,x)c(1).*x;
            Lpfit=nlinfit(mu_theory.*store_all_acdss,store_all_bfit,Fit_Lp_Function,Lp_guess);
            
            figure;
            mu_ds = mu_theory.*store_all_acdss;
            mu_ds = mu_ds*pixelsize*1E-9*1E6;
            Lambda = store_all_bfit;
            Lambda = Lambda*(pixelsize*1E-9*1E6)^2;
            Lambda_th = mu_theory.*Lpfit.*store_all_acdss;
            Lambda_th = Lambda_th*(pixelsize*1E-9*1E6)^2;
            plot(mu_ds,Lambda,'or',mu_ds,Lambda_th)
            xl=xlabel('\mu \Delta s (\mu m)','FontSize',14);
            yl=ylabel('\Lambda (\mu m^2)','FontSize',14);
            set(gcf,'color','w');
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', 14);
            %yl.Position(1) = -0.55;
            set(findall(gca, 'Type', 'Line'),'LineWidth',4);
        end
        
        Store_Lpfit(vvv,1) = Lpfit;
    end
    
    Result_Lp = mean(Store_Lpfit);
    std_Result_Lp = std(Store_Lpfit);
    
    Result_Lp = Result_Lp*pixelsize*1E-9*1E6; % Convert Lp units from pixel to micron
    std_Result_Lp = std_Result_Lp*pixelsize*1E-9*1E6; % Convert Lp units from pixel to micron
    fprintf('Lp = %4.f %s %4.f \x3BCm \n', Result_Lp, char(177), std_Result_Lp) % print Lp on screen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif(BS == 0)
    %   Count m
    countm = 0;
    
    %   Loop over curvature calculation for all skipping indices (m)
    for m = mstart:1:mend  %skipping over m-1 nodes
        % print m on screen to keep track of progress
        m
        
        %       At the start of each loop, deallocate all matrices except the followings
        clearvars -except pixelsize foldername filenameprefix y_int mu_eq_1p5 Store_Lpfit filament N_BS BS_size num_filament vvv  mu_eq_1p5 L y_int h stats bguess mstart mend Vstart Vend Nstart Nend binsize m Numbern mu Lp_guess countm store_all_acdss store_all_bfit write_bfitvsmeands write_curvdist m_write_curvdist exp_mu;
        
        % Declare L_av for average filament length calculation
        L_av=0;
        
        % kk is an index for storing curvature
        kk=1;
        
        % sizedx is an index for storing the bond spacing on the filament
        sizeds =1;
        
        %       Loop over filaments
        for v = Vstart:Vend
            
            %           Clear all necessary matrices
            clear theta x y xx yy ds C D E F xskip yskip store state storej;
            
            %           Write coordinates file name to a matrix
            filename = sprintf('./%s%04d.dat', filenameprefix, v);
            
            %           Load coordinates from a file into array A
            A = load(filename);
            
            %           B stores the size of the array A >> example output B = [95 2]
            B = size(A);
            
            %           M stores the number of coordinate points (number of bonds)
            %           pick the first column from B
            M = B(1);
            
            %           Store x and y coordinates from the files
            for j = 1:M
                x(j) = A(j,1);
                y(j) = A(j,2);
            end
            
            %           Calculate bond spacing in x- and y-projection
            
            for j = 1:M-1
                xx(j) = x(j+1)-x(j); % displacement in x
                yy(j) = y(j+1)-y(j); % displacement in y
            end
            
            %           Calculate the bond spacing using pythagoras
            for j = 1:M-1
                ds(j) = sqrt((xx(j))^2 + (yy(j))^2);
            end
            
            %           Calculate length of the filament
            L(v) = sum(ds);
            %-------------------------------------------
            
            %           Loop for doing recursion
            %           Recursion is the proccess to get more data out of a filament
            %           If we skipped some nodes and measured the curvature, the
            %           information on those skipped nodes would be wasted
            %           Therefore, after measuring the curvature for the first sets of
            %           coordinates on the filament, starting from the first point, we
            %           now start skipping and measuring the curvature from the second
            %           point on the chian (This point was previously wasted during the
            %           last calculation.
            
            for q = 1:m % recursion  % If no-recursion, use "for q = 1:1"
                
                %             Deallocate necessary array
                clear C D dss E F
                
                %               Loop to find out how many big bonds (after skipping) the
                %               filament can have
                b = floor ((M-q)/m);
                
                % Loop to find the bond spacing of the big bonds
                for i = 1:b
                    C(i) = x(q+(m)*i)-x(q+(m)*(i-1));  %displacement in x
                    D(i) = y(q+(m)*i)-y(q+(m)*(i-1));  %displacement in y
                    dss(i) = sqrt(C(i).^2 + D(i).^2);  %bond spacing
                    acdss(sizeds) = dss(i); % store bond spacing
                    sizeds = sizeds+1; % ready to store the next one
                end
                
                % Loop to calculate curvature
                for i = 1:b-1
                    
                    if (C(i+1)>C(i)) && (D(i+1)>D(i)) %if pointing in 1st quadrant
                        dtheta_all(i) = acos((C(i)*C(i+1) + D(i)*D(i+1))/(dss(i)*dss(i+1))); %calculate angle
                    elseif (C(i+1)<C(i)) && (D(i+1)<D(i)) %3rd quadrant
                        dtheta_all(i) = acos((C(i)*C(i+1) + D(i)*D(i+1))/(dss(i)*dss(i+1))); %calculate angle
                    else %2nd and 4th quadrant
                        dtheta_all(i) = -acos((C(i)*C(i+1) + D(i)*D(i+1))/(dss(i)*dss(i+1))); %calculate angle
                    end
                    
                    %Calculate and store curvature from the angle and bond
                    %length
                    curvature(kk)= dtheta_all(i)/((dss(i)+dss(i+1))*0.5);
                    kk=kk+1; %ready to store the next one
                end
            end % End recursion loop
        end % End loop over filament
        
        % After storing all curvature values for current m, subtract the mean and take
        % the absolute value
        curvature = real(curvature);
        %curvature = abs(curvature-mean(curvature));
        
        [h, stats] = cdfplot(curvature);
        hold on;
        [cdf_y,cdf_x] = ecdf(curvature);
        bfit_guess = [getfield(stats, 'mean') getfield(stats, 'std')];
        cumulative_gaussian = @(bfit,x)0.5*(1+erf((x-bfit(1))/(sqrt(2)*bfit(2))));
        bfit=nlinfit(cdf_x,cdf_y,cumulative_gaussian,bfit_guess);
        cdf_y_fit2 =  normcdf(cdf_x,bfit(1),bfit(2));
        plot(cdf_x,cdf_y_fit2,'r');
        hold off;
        
        store_all_bfit(m-mstart+1) = 1/(bfit(2))^2;
        store_all_acdss(m-mstart+1) = mean(acdss);
        
        %To write measured and fitted curvature distribution functions
    end
    
    % Calculate average length of all filaments
    L_av = sum(L)/length(L);
    
    % Assign a matrix containing all m values
    Numberm=mstart:mend;
    %calculate mu for all m
    
    if(mu_eq_1p5)
        mu_theory = 1.5;
    else
        mu_theory=Numberm./(((2.*Numberm-1).*(1-1./Numberm)/3)+1);
    end
    
    % Estimate gauss persistence length (guess Lp)
    Lp_guess = store_all_bfit(end)/(mu_theory(end)*store_all_acdss(end));
    
    if(y_int)
        Fit_Lp_Function = @(c,x)c(1).*x+c(2);
        Lp_guess = [Lp_guess,0];
        fit=nlinfit(mu_theory.*store_all_acdss,store_all_bfit,Fit_Lp_Function,Lp_guess);
        fit_lambda = fit(1).*mu_theory.*store_all_acdss+fit(2);
        Lp = fit(1);
        y_intercept = fit(2);
        figure;  % plot CDF of curvature distribution
        plot(mu_theory.*store_all_acdss,store_all_bfit,'or',mu_theory.*store_all_acdss,fit_lambda); set(gca,'FontSize',15)
    else
        Fit_Lp_Function = @(c,x)c(1).*x;
        Lpfit=nlinfit(mu_theory.*store_all_acdss,store_all_bfit,Fit_Lp_Function,Lp_guess);
        
        figure;
        mu_ds = mu_theory.*store_all_acdss;
        mu_ds = mu_ds*pixelsize*1E-9*1E6;
        Lambda = store_all_bfit;
        Lambda = Lambda*(pixelsize*1E-9*1E6)^2;
        Lambda_th = mu_theory.*Lpfit.*store_all_acdss;
        Lambda_th = Lambda_th*(pixelsize*1E-9*1E6)^2;
        plot(mu_ds,Lambda,'or',mu_ds,Lambda_th)
        xl=xlabel('\mu \Delta s (\mu m)','FontSize',14);
        yl=ylabel('\Lambda (\mu m^2)','FontSize',14);
        set(gcf,'color','w');
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 14);
        %yl.Position(1) = -0.55;
        set(findall(gca, 'Type', 'Line'),'LineWidth',4);
    end
    
    Result_Lp = Lpfit;
    Result_Lp = Result_Lp*pixelsize*1E-9*1E6; % Convert Lp units from pixel to micron
    fprintf('Lp = %4.f \x3BCm \n', Result_Lp) % print Lp on screen
end

