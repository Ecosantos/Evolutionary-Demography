% Lamvaresti.m: a program to estimate the sampling variance in log stochastic lambda
%  using approximation formulae from Doak et al. 2005 (equation numbers refer to this 
%  paper and its appendix 3). 

% You must have the symbolic math toolbox of Matlab  to use this program. 
% This program uses two functions (secder.m and eigenall.m) from the website of 
%  programs that accompany Morris and Doak (2002): www.sinauer.com/PVA/
% The general form of data entry used here is quite similar to other, simpler
%  programs also on this website, including Vitalsens.m and Stochsens.m; reading 
%  through these programs may help you understand the structures used here if you
%  having trouble.
%  As a working example, this sample program is written to run a model for the desert
%  tortoise life history and data that are also explained in the Morris and Doak book. 

% One warning: the symbolic logic routines and the simulations to estimate correlations
%  in beta variable means and variances are time-consuming, with one to several minutes
%  between different steps. Be patient. 

clear all;
global yrsam kknums mmnums     % global variables used by called functions
randn('state',sum(100*clock)); % seeding random numbers
rand('state',sum(100*clock));  % seeding random numbers
warning off                    % MATLAB:divideByZero

%***************** Parameters that must be input by user ***********************
% First, give symbolic names for each vital rate to be used in this program. For the desert tortoise, 
%  these are: first, six survival rates (for stages 2-7); next, 5 growth rates (stages 2-6); 
%  and, finally, three fecundities (stages 5-7). 
% These symbolic definitions are given below, and then the vector of these names (Svr) is defined.
syms  v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14  % vital rates as symbolic variables
Svr = [v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 ]; % vector of symbolic vital rates 
% Next, give the mean Vital rate values:
realvrmeans= [0.77 0.95		0.79		0.93		0.91		0.87 ...
0.3		0.31		0.19		0.22		0.052	0.42	0.69	0.69]; 
% Then estimated true temporal variances (not standard deviations) of the Vital Rates:
realvrvars=[0.0006		0.0006		0.053		0.0006		0.049 ...
0.02		0.0096 0.017		0.048		0.0006		0.0006	0	0	0];
% Next, you must say what the distribution is for each vital rate: this program only distinguishes 
%  between beta-distributed variables (coded as 1) and all others, assumed to be fecundities 
%  or similarly distributed parameters (coded 2). 
vrtypes= [ones(1,11),2,2,2]; 
% Then, you must give a the full estimated matrix of temporal correlations between the vital rates.
% We do this here by putting the matrix for the desert tortoise directly in the code (see also
%  Table 8.2 in Morris and Doak 2002). You could also load a matlab binary data file that has your
% correlation matrix.  
realcorrmx = ...
 [  1.0000    0.4690    0.3818   -0.5968   -0.0140   -0.9997   -0.7044    0.8978    0.9561   -0.5140   -0.9941         0         0         0;
    0.4690    1.0000    0.9952    0.4289    0.8764   -0.4961   -0.9571    0.8100    0.1891    0.5166   -0.5630         0         0         0;
    0.3818    0.9952    1.0000    0.5139    0.9189   -0.4103   -0.9250    0.7500    0.0937    0.5967   -0.4809         0         0         0;
   -0.5968    0.4289    0.5139    1.0000    0.8108    0.5714   -0.1494   -0.1820   -0.8059    0.9951    0.5050         0         0         0;
   -0.0140    0.8764    0.9189    0.8108    1.0000   -0.0171   -0.7000    0.4281   -0.3069    0.8651   -0.0959         0         0         0;
   -0.9997   -0.4961   -0.4103    0.5714   -0.0171    1.0000    0.7260   -0.9109   -0.9463    0.4870    0.9968         0         0         0;
   -0.7044   -0.9571   -0.9250   -0.1494   -0.7000    0.7260    1.0000   -0.9450   -0.4650   -0.2469    0.7779         0         0         0;
    0.8978    0.8100    0.7500   -0.1820    0.4281   -0.9109   -0.9450    1.0000    0.7289   -0.0834   -0.9407         0         0         0;
    0.9561    0.1891    0.0937   -0.8059   -0.3069   -0.9463   -0.4650    0.7289    1.0000   -0.7430   -0.9179         0         0         0;
   -0.5140    0.5166    0.5967    0.9951    0.8651    0.4870   -0.2469   -0.0834   -0.7430    1.0000    0.4167         0         0         0;
   -0.9941   -0.5630   -0.4809    0.5050   -0.0959    0.9968    0.7779   -0.9407   -0.9179    0.4167    1.0000         0         0         0;
         0         0         0         0         0         0         0         0         0         0         0         1         0         0;
         0         0         0         0         0         0         0         0         0         0         0         0         1         0;
         0         0         0         0         0         0         0         0         0         0         0         0         0         1 ];                       
% Define how the different vital rates combine to make each matrix element, doing this 
%  by defining the entire symbolic matrix:
symmx = [0  0           0           0           0            v12            v13             v14 
        v1  v1*(1-v7)   0           0           0              0            0               0
        0   v1*v7       v1*(1-v7)    0          0              0            0               0
        0   0           v1*v7       v2*(1-v8)   0               0            0               0
        0   0           0           v2*v8        v3*(1-v9)      0           0               0
        0   0           0           0           v3*v9            v4*v10      0                0
        0   0           0           0           0               v4*v10       v5*(1-v11)     0
        0   0           0           0           0               0            v5*v11         v6];
    
% Now, what are the sampling intensities for each vital rate and the durations of sampling that 
%  you want to have run calculations for? insams is a matrix with columns of sampled number of 
%  individuals used to estimate each vital rate (in the same order as for the means and variances
%  above) and rows for different sets of these samples to run. For example, the insams defined 
%  below has one set of sampling of 30 individuals for each vital rate, and one set of sampling 100 
%  individuals for each rate; remember that these sampling patterns can be those used or ones you 
%  might want to consider. 
insams = [ones(1,14)*30;  ones(1,14)*100];
% Then input each sampling duration that you want to consider: each number here is one duration to try: 
yrsams = [3 5 10 20];
outputfilename = 'EstiVarTort.txt'; % The name of the file to save output data to
%*************** End of Parameter inputs: Proceeding to calculations **************************
  
%First Step: Basic calculations and estimation of the deterministic vital rate sensitivities
estiouts=[];                            % The variable to store output data
realmx = subs(symmx,Svr,realvrmeans);   % Making a matrix of the mean numerical values
nmx = length(realmx);                   % Size of pop mx.
nvr = length(realvrvars);               % Number of vital rates

[lambdas,lambda1,W,w,V,v]= eigenall(realmx);    % Use eigenall.m to get eigenvalues
sensmx = v*w'/(v'*w);                           % Get sensitivities of matrix elements
vrsens = zeros(1,nvr);                          % Initialize vital rate sens.
for xx=1:nvr  % A loop to calculate sensitivity for each vital rate
	% First get derivatives of elements with respect to vital rates:
	diffofvr = double(subs(diff(symmx,Svr(xx)),Svr,realvrmeans));
    vrsensbyelements(:,:,xx) = diffofvr; 
    % Then, sum up to get row of total vital rate sensitivities:
 	vrsens(xx) = double(sum(sum(sensmx.*diffofvr))); 
end; % xx

% Second Step: Calculate stochastic lambda and its sensitivities to the matrix element means
mx = realmx; % Set mx equal to the name of stored pop'n matrix 
vrcovmx = realcorrmx.*(sqrt(realvrvars')*sqrt(realvrvars)); % Make a covariance matrix
tau=(vrsens)*vrcovmx*(vrsens'); % tau as in Tuljapurkar (1991), but estimated by vital rates
% Estimate  log(lambda_S), the log of stochastic lambda:
loglamS = log(lambda1) - 0.5*(1/(lambda1^2))*tau; 

squloglamderivs=[]; % Here, we are define the three storage variables for the final calcs:
squVarsums =[];
squCorrsums = [];

for ii=1:nvr % Loop to get the values needed to estimate the derivatives: d(log(lambda_S))/d(vi)
    kkllsum=0;
    for kk=1:nvr
        for ll = 1:nvr
            dSldi =0; 
            dSkdi =0;           
dSldi = sum(sum( sensmx.*double(subs(diff(diff(symmx,Svr(ll)),Svr(ii)),Svr,realvrmeans)) ));
dSkdi = sum(sum( sensmx.*double(subs(diff(diff(symmx,Svr(kk)),Svr(ii)),Svr,realvrmeans)) ));
            for xx = 1:nmx
                for yy =1:nmx
dSldi = dSldi + vrsensbyelements(xx,yy,ii)*sum(sum(secder(mx,xx,yy).*vrsensbyelements(:,:,ll) ));
dSkdi = dSkdi + vrsensbyelements(xx,yy,ii)*sum(sum(secder(mx,xx,yy).*vrsensbyelements(:,:,kk) ));
                end
            end
            kkllsum= kkllsum + vrcovmx(kk,ll)*(dSldi*vrsens(kk) +dSkdi*vrsens(ll));
        end 
    end
    % The derivatives of log(lambda_S) with respect to each vital rate:    
    loglamderivs(ii) =     vrsens(ii)/lambda1 + vrsens(ii)*tau/(lambda1^3) -  kkllsum/(2*lambda1^2);
    % The square of each derivative, which multiples with the variance in each rate in equation 2. 
    squloglamderivs(ii) = (loglamderivs(ii))^2; 
    % The sums that multiple with the variances of the variances terms in equation 2:
    squVarsums(ii) = (1/lambda1^4)*(  sum( vrsens(ii)*vrsens.*sqrt(realvrvars).*realcorrmx(ii,:)) )^2;
    
     disp('The vital rate number and sensitivity of log(lambda_S) to this vital rate');
     disp([ii,loglamderivs(ii)]);
end; %ii
% Finally, the matrix of values that multiple with the variances of correlations in equation 2:
squCorrsums = (1/lambda1^4)*((sqrt(realvrvars')*sqrt(realvrvars)).*(vrsens'*vrsens)).^2;
clear v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11 v12 v13 v14 Svr symmx; %making space in memory

% Third Step: estimate sampling variance in log(lambda_S)for different sampling patterns 
for ii = 1:length(insams(:,1))% Loop through each set of sampling intensities
    SamNs = insams(ii,:); % The vector of within year sample sizes to use 
    for jj=1:nvr            % A loop to use simulation to estimate the correlation of means and standard 
        if vrtypes(jj) ==1  %  deviations in sampled values for beta-distributed variables:
            mn = realvrmeans(jj);
            va = realvrvars(jj);
            vv = mn*((mn.*(1-mn)/(va))-1); % calculate the beta parameters
            ww = (1-mn).*((mn.*(1-mn)/(va))-1);
             aa = betarnd(vv,ww,SamNs(jj),10000); % Draw 10,000 sets of values
            aavars = var(aa);
            aaSD= sqrt(aavars);
            aameans=mean(aa);
            
            aacov=cov([aaSD',aameans']);
            vrvrvarcovs(jj) = aacov(1,2);
        else vrvrvarcovs(jj)=0;
        end;
        betacorrcontribut(jj) = 2*vrvrvarcovs(jj).*loglamderivs(jj).*(1/lambda1^2).*(sum( ...
            vrsens(jj)*vrsens.*sqrt(realvrvars(jj)).*realcorrmx(jj,:)) );
        disp('The vital rate number and, next line, beta-value correlation contribution to variance');
        disp(jj); disp(betacorrcontribut(jj));

    end;   
    clear aa aavars aaSD aameans aacov; % making space in memory
    for yy=1:length(yrsams); % Loop through the sampling durations
        yrs = yrsams(yy); % number of years of data
        
        for xx=1:nvr %loop to estimate  within-year sampling variances of each vital rate:
            if vrtypes(xx) == 1;   inyrvar(xx) = realvrmeans(xx)*(1-realvrmeans(xx));  end; % binomials
            if vrtypes(xx) == 2;   inyrvar(xx) = realvrmeans(xx);   end; % using Poisson variance for fecundities
        end
        % Next, estimate the total sampling variance for mean values (equation A6):
        meanvars = (1/yrs).*(realvrvars + inyrvar./SamNs);   
        % Then, the variances for the corrected variance estimates (equation A9):
        correctedvarvars = (2*yrs/(yrs-1)^2)*realvrvars.*(realvrvars + 2*(inyrvar./SamNs));                        
        SDvars = (correctedvarvars./(4.*realvrvars)); % Transform correctedvarvars to get variances of SDs
        SDvars(isnan(SDvars)) = 0;
        correlvars = (yrs/(yrs-1)^2)*(realcorrmx.^2 -1).^2; % The variances of the correlations
        
        % At Last, get the outputs: 
        % 1. The sampling variance in the estimate of deterministic log(lambda): this is also the 
        %     sampling variance in log(lambda_S) generated by sampling variance of the mean vital rates: 
        DeterLogLamVar = sum(squloglamderivs.*meanvars);
        % 2. Sampling variance of log(lambda_S) from just variance in means and variances of vital rates: 
        VarLogLamVar   = sum(squloglamderivs.*meanvars + squVarsums.*SDvars );
        % 3. Sampling variance of log(lambda_S) from variances of means, variances, and correlations, 
        %      but without the effects of beta variable correlations
        FullLogLamVar = sum(squloglamderivs.*meanvars + squVarsums.*SDvars + 0.5*sum(squCorrsums.*correlvars) );
        % 4. The best of sampling variance of log(lambda_S) with the effects of beta variable correlations
        FullLogLamVarADDED = FullLogLamVar+sum(betacorrcontribut);
        
        % Save the data: as now written, the outputs are one row for each combination of sampling duration and 
        %  intensity. The columns of data are: sampling intensity for the first vital rate; sampling
        %  duration; sampling variance (SV) for deterministic log(lambda); SV for log(lambda_S) from SV in 
        %  vital rate means and variances; SV for log(lambda_S) from SV in means, variances, and correlations; 
        %  SV for log(lambda_S) from all sources; estimated log(lambda_S) for the input parameters; and, 
        %  estimated log(deterministic lambda). 
        estiouts = [estiouts;[SamNs(1) yrs  DeterLogLamVar VarLogLamVar FullLogLamVar FullLogLamVarADDED  ...
                    loglamS  log(lambda1)]];
        disp('The sampling intensity set, sampling duration set, and sampling variance in log(lambda_S)');
        disp([ii,yy, FullLogLamVarADDED]);           
    end; %yy
end; %ii
save(outputfilename, 'estiouts','-ASCII'); % This saves a file with the data in estiouts
disp('DONE!');

