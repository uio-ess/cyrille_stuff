% analysis of the data OCL
% HV10
% spectra, images and current

clear

%% calibration of the camera : photons / count
% dr = '/Users/cyrillethomas/Documents/ESS Instrumentation Projects/A2T/OCL experiments/lab/2018-02-27-camera-calibration/'
% 
% calib= load([dr 'camera_OCL_calibration_27022018'])
% ph_count = 1./calib.pp(1)
% 
ph_count = 5.7817 % ph/count

%%
% elem charge (C)
qe = 1.6e-19 ;


%%
% range of file numbers 
aa= [4173 4175:4183]

% init loop 
N = length(aa)
cam_exp = zeros(1,N);
light_counts = zeros(1,N);
Charge_i = zeros(1,N);
Proton_N = zeros(1,N);
BK =  zeros(1,N) ;

%% loop files and extraction of the data 

for ii=1:10
    
    % generate the file number
    data_nam = strcat(sprintf('%16.16d-HV10.h5', aa(ii)) )
    
    
    % read the image data 
    im = h5read(data_nam,'/data/images/CAM1/data') ;
    % get the background of the image from one corner
    im_bk = im(1:100,1:100) ;
    BK(ii) = mean(im_bk) ;
    
    % read the current and time vectors
    cur =  h5read(data_nam,'/data/wavefront/ps4264py/y_data') ;
    t_cur =  h5read(data_nam,'/data/wavefront/ps4264py/x_data') ;
    
    % read exposure of the camera
    t_exp = hdf5read(data_nam,'/data/images/CAM1/','CAM1:det1:AcquireTime_RBV') ;
    cam_exp(ii) = t_exp ;

    % total counts less background 
    light_counts(ii) = sum(im(:))-BK(ii) ; 
    
    % find time index from trigger to cam expos length
    id_c = find(t_cur>0 & t_cur<t_exp) ;
    
    % sample time 
    dt = mean(diff(t_cur)) ;
    
    % charge and protons per image exposure
    Chrg  = sum(-cur(id_c) * dt );
    Charge_i(ii) = Chrg ; 
    Proton_N(ii) = Chrg  / qe ; 
    
    %%
    figure(1)
    imagesc(im)
    drawnow
    
    %%
    
    figure(2)
    plot(t_cur(:)*1e3,-(cur(:)/1e-9))
    drawnow
    % ylim([0 10])
    
    % xlim([0 1])
    
    
end

%%

fts = 22
%%
figure(4)

plot(Proton_N,light_counts,'s','markersize',15)
xlabel('Protons' , 'fontsize',fts)
ylabel('Counts' , 'fontsize',fts)

set(gca,'FontSize',fts)
%%
N_counts = [mean(light_counts) std(light_counts) min(light_counts) max(light_counts)]

N_P = [mean(Proton_N) std(Proton_N) min(Proton_N) max(Proton_N)]


%%
N_counts./N_counts(1)

N_P./N_P(1)
%%
T_len = 0.95 
D = 50 /1.4 
L = 1185
Omega = ( D/2 / L )^2 / 4 

ph_Pr = N_counts(1) * ph_count  ./N_P(1) / Omega / T_len

%%
figure(4)

plot(Charge_i ./ 1e-9 ,light_counts * ph_count  ./N_P(1) / Omega / T_len,'s','markersize',15)
xlabel('Charge (nC)' , 'fontsize',fts)
ylabel('Photons' , 'fontsize',fts)

set(gca,'FontSize',fts)


