% These parameters controls all aspects of GPS_SW_RCX.

% Hardware properties
FS = 8.192e6;                        %sampling frequency [Hz]
FC = 2.048e6;                        %incoming center frequency [Hz]
IF_BANDWIDTH = 2.5e6;                %IF filter bandwidth, used to calculate CNo
R_INS_APC = [-2e-2,0e-2,0e-2]';      %Body frame vector from acc. to antenna phase center

% Detection control
USE_AIDING = 1;                      %directive to use external data to aid detection and tracking
N_CODES_AQU = 60;                    %number of codes to us in acquisition
DB_DETECTION = 26;                   %CNo detection threshold [dB]
FREQ_STEP = 1;                       %freq. steps used in searching doppler shifts [Hz]
FD_SIZE = 5;                         %magnitude of max. doppler shifts to search through [Hz]

% Correlator settings
N_CORR = 2*14+1;                     %number of correlators, must be odd
MAX_SPACE = 1.0;                     %maximum correlator spacing, in chips
SEARCH_WIN = 0.5;
CF_TYPE = 'tent';                    %type of curve fitting - 'tent' or 'poly'
COH_INT_TIME = 10e-3;                %coherent integration time [seconds]
COH_INT_SAM = round(COH_INT_TIME*FS);%number of samples during the integration time
CORR_LOSS = -0.0;                    %correlation losses, used to calculate CNo [dB]

% Tracking controls
D1FDT1 = 20;                         %d1f/dt1, used to optimize 2nd-order PLL [Hz/s]
D2FDT2 = 2e3;                        %d2f/dt2, used to optimize 3rd-order PLL [Hz/s^2]
HNUM = 1e-3;                         %DLL convergence rate, [0 1]
MAGNITUDE_SMOOTHING = 5e-2;          %coefficent for moving average magnitude estimate smoothing, [0 1]
DB_LIMIT_SMT = 6;                    %if mean SNR is below this, we are not tracking [dB]
DB_LIMIT_INS = 3;                    %if instananeus SNR is below this, we are not tracking [dB]
LOL_M = 6;                           %lock at last M points to determine track/dont-track
LOL_N = 0;                           %if less than N were valid, ignore mesurement
LOL_X = 1;                           %if less than X were valid, reset doppler frequency (0 to disable)

% --PLL controls
USE_PLL = 1;                         % To use only the FLL, set USE_PLL = 0.  If you want to use the PLL, set this to 1
PLL_LOOP_ORDER = 2;                  % This is the loop order to use for the PLL
PLL_SWITCH_TIME = 0.0;               % This is the time at which to switch over to the PLL [seconds]

% --FLL controls
FLL_BANDWIDTH = 20;                  %bandwidth in Hz
A_FLL = (1.89*FLL_BANDWIDTH).^2;     %DE = w_nF^2 = (1.89*B_LF)^2
B_FLL = (sqrt(2)*1.89*FLL_BANDWIDTH);%EE = sqrt(2)*w_nF = sqrt(2)*1.89*B_LF

% Data files
work_dir        = 'D:\Desktop\data\20120308_surveyed_point\';   %local var
RAW_FILE        = [work_dir 'raw\data.log']; %raw data file name
MAIN_FILE       = ['']; %main memory file name
TRACK_DIRECTORY = [work_dir '5msec\'];   %directory that stores tracks
IS_COMPLEX      = 0;                      %directive to read complex data

% Peripheral data files
BIT_START_TIME_FILE = [work_dir 'bst.mat'];     %name of bit-start-time file
RECEIVER_FILE   = [work_dir 'rec.mat'];         %name of receiver information file
DATA_PLOTTER    = '..\RAIN_Electrical\Data Plotter'; %addres to Data Plotter by Joakim, to read accels.

% Position solution
POSITION_DB_LIMIT = 10;              %minimum SNR needed to be included in position solution [dB]
POSITION_EXL_RNG = 10;               %exclution range around bad samples

% Velocity solution
DOPPLER_DB_LIMIT  = 10;              %minimum SNR needed to be included in doppler [dB]
DOPPLER_MAX_PHI_ERR = 60;            %maximum phase tracking error [deg]
DOPPLER_EXL_RNG = 15;                %exclution range around bad samples

% Fixed parameters, do not change
W_FC = 2*pi*FC;                      %angular center frequency
TP = 1/FS;                           %sample spacing
T = 1e-3;                            %duration of a C/A code period [seconds]
FSAMP_MSEC = FS*T;                   %number of samples in 1 millisecond
ONE_MSEC_SAM = round(FSAMP_MSEC);    %~number of samples in one millisecond of data
CHIPS_PER_CODE = 1023;               %number of chips per code period
CHIP_RATE = 1023e3;                  %code chipping rate
L1 = 10.23e6*154;                    %L1 frequency [Hz]
C = 299792458;                       %Speed of light [m/s]
INITIALIZE = 1;                      %directive to initialize persistent variables