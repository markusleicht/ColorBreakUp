function cbu_loop

% AUTHOR: Markus Leicht
% EMAIL: markusleicht84@gmx.de
% DATA: 2022-02-25  
% VERSION: 1.0.0.0
% DOWNLOAD: https://github.com/markusleicht/ColorBreakUp 
%
% BACKGROUND:
% The function CBU_LOOP is a supplement to the main function CBU_MODEL (also available via GitHub). Both functions were written during the author's doctorate. The associated doctoral thesis refers to the the basic model and the supplementary loop and provides further information on the subject of color break-up (CBU): Leicht, M. (2022). Perception of Color Break-Up [Doctoral thesis]. Ilmenau University of Technology.
%
% PURPOSE:
% The CBU_LOOP allows to calculate the model output from the basic function CBU_MODEL (i.a. model indices: cbu score, non-cbu score, reference score) under variation of one input parameter (while keeping all other parameters constant) in order to determine the effect of the variable parameter on cbu perception.
% The following parameters can be varied in the corresponding code segments: Segment A - frame rate, Segment B - eye movement velocity, Segment C - duty cycle, Segment D - content width, Segment E - subframe number.
%
% BENEFIT:
% Every single Segment A to E of the CBU_LOOP is executed separately, accessing the main model function CBU_MODEL (see Segments II to XI; directly implemented within the loop function; line 125-1964) and reading in the constants and the variable parameter (range/steps) that are specified in the code lines of the CBU_LOOP, i.e., the CBU_LOOP allows to automatically run the CBU_MODEL multiple times for the defined conditions. Otherwise the variable parameter would have to be adjusted manually within the CBU_MODEL and the code of the CBU_MODEL would have to be executed manually a multiple times.
%
% STRUCTURE:
% Segment A - Variation of Frame Rate
% Segment B - Variation of Horizontal Eye Movement Velocity (unfinished)
% Segment C - Variation of Duty Cycle (unfinished)
% Segment D - Variation of Content Width (unfinished)
% Segment E - Variation of Subframe Number (unfinished)
%
% INPUT (EXEMPLARY FOR SEGMENT A):
% Step 1 - Determination of the constant parameter (line 69-73: subframe number, duty cycle, content width, eye movement velocity, content movement velocity). 
% Step 2 - Definition of the min-max range and the steps between the discrete values of frame rate as the variable parameter (line 78-80).
% Step 3 - Determination of the basic parameter within the CBU_MODEL (Segment II: Data Entry Mask; line 139-303).
% Step 4 - The CBU_LOOP function has to be executed by running single segments exclusively, not the whole function. For a variation of the frame rate, Segment A has to be executed (current model input refers to the example of choice that is discussed in Chapter 5 of the associated doctoral thesis).
%
% DATA OUTPUT (EXEMPLARY FOR SEGMENT A):
% After running the segment, it's data input/output is summarized within an overview table in the command window. The table lists the loop input, i.e., subframe number [-], duty cycle [-], content width [px], frame rate [Hz], and eye/content movement velocity [px/fr], and the loop output, i.e., cbu score - option 1 [deg2], non-cbu score - option 1 [deg2], reference score - option 1 [deg2], cbu score - option 2 [-], non-cbu score - option 2 [-], and reference score - option 2 [-], for every single loop (loop number).
%
% GRAPHIC OUTPUT (EXEMPLARY FOR SEGMENT A):
% Figure 1: Model Indices in Dependency of Frame Rate (Option 1) >> refers to listed loop output in overview table (output before standardization and filtering)
% Figure 2: Model Indices in Dependency of Frame Rate (Option 2) >> refers to listed loop output in overview table (output before standardization and filtering)
% Figure 5: Empirical vs. Theoretical CBU (Option 1, w/o Stevens)
% Figure 6: Empirical vs. Theoretical CBU (Option 2, w/o Stevens)
% Figure 7: Empirical vs. Theoretical CBU (Option 2, Stevens)
% Figure 12: Scatter Plot CBU Comparison (Sample, Option 1, Standardized, w/o Stevens)
% Figure 13: Scatter Plot CBU Comparison (Sample, Option 2, Standardized, w/o Stevens)
% Figure 14: Scatter Plot CBU Comparison (Sample, Option 2, Standardized, Stevens)
% Figure 15: Boxplot CBU Comparison (Sample, Option 1, Standardized, w/o Stevens)
% Figure 16: Boxplot CBU Comparison (Sample, Option 2, Standardized, w/o Stevens)
% Figure 17: Boxplot CBU Comparison (Sample, Option 2, Standardized, Stevens)
%
% DESCRIPTION:
% A basic code description is embedded in the MATLAB function itself. However, the fundamental structure of the loop and its build-up are described in the corresponding doctoral thesis (see Chapter 5).
%
% COMPATIBILITY:
% The function was written with MATLAB Version R2018b (OS: Windows 10, 64bit).      
%
% MATLAB SEARCH PATH
% For the graphic output, the external functions >>cprintf<< (https://de.mathworks.com/matlabcentral/fileexchange/24093-cprintf-display-formatted-colored-text-in-command-window), >>blandaltmanplot<< (https://de.mathworks.com/matlabcentral/fileexchange/71052-blandaltmanplot), and >>boundedline<< (https://de.mathworks.com/matlabcentral/fileexchange/27485-boundedline-m) must be available by adding the folders to the MATLAB search path.
%
% Copyright 2022 Markus Leicht

%% SEGMENT A - VARIATION OF FRAME RATE [HZ]
% first step - determination of the constants (subframe number, duty cycle, content width, eye movement velocity, content movement velocity) 
% second step - definition of the min-max range (min./max.) and the steps between the discrete values of frame rate as the variable
% third step - determination of the basic parameter within the CBU_MODEL (Segment II: Data Entry Mask)
% fourth step - run segment exclusively, not the whole function (!)

clearvars % deletes all parameter in order to set the status to zero

% CONSTANTS >> values that are not changed during the execution of the loop
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subframe_number = 3.0;
duty_cycle = 0.30;
content_width_px = 40.0;
eye_velo_hor_pxframe_fix = 24.0; % defined eye movement velocity value at FRAME_RATE_HZ_MIN as the lowest frame rate applied (see VARIABLE), eye movement velocity in [PX/FR] depends on the applied frame rate, therefore, the eye movement velocity in [PX/FR] can only be determined for a specific frame rate, during the execution of the loop for various frame rates (see VARIABLE) the eye movement velocity in [PX/FR] is changed (see LOOP_TABLE) in order to achieve a stable eye movement velocity in [DEG/SEC] for various frame rates since the eye movement velocity should be a constant and not a variable in this segement 
content_velo_hor_pxframe_fix = 24.0; % defined content movement velocity value at FRAME_RATE_HZ_MIN as the lowest frame rate applied (see VARIABLE), content movement velocity in [PX/FR] depends on the applied frame rate, therefore, the content movement velocity in [PX/FR] can only be determined for a specific frame rate, during the execution of the loop for various frame rates (see VARIABLE) the content movement velocity in [PX/FR] is changed (see LOOP_TABLE) in order to achieve a stable content movement velocity in [DEG/SEC] for various frame rates since the content movement velocity should be a constant and not a variable in this segement
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% VARIABLE >> value that is changed during the execution of the loop (frame rate needs to be a positive integer equal or greater than 20Hz, no upper limit) 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
frame_rate_hz_min = 30.0; 
frame_rate_hz_max = 420.0; 
frame_rate_hz_step = 5.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% LOOP TABLE >> includes constants and variables as well as the resulting model indices
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
loop_table = zeros((frame_rate_hz_max-frame_rate_hz_min)/frame_rate_hz_step+1,12); 
loop_table(1:end,1) = subframe_number; 
loop_table(1:end,2) = duty_cycle; 
loop_table(1:end,3) = content_width_px;

aa = 0; % adaption of eye and content movement velocity in [PX/FR] for changing frame rate in order to guarantee a stable eye and content movement velocity in [DEG/SEC]

for bb = frame_rate_hz_min:frame_rate_hz_step:frame_rate_hz_max

aa = aa + 1;    
    
loop_table(aa,4) = bb;

loop_table(aa,5) = eye_velo_hor_pxframe_fix/(bb/frame_rate_hz_min);

loop_table(aa,6) = content_velo_hor_pxframe_fix/(bb/frame_rate_hz_min);

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CONDITION (1) >> for calculation of theoretial cbu scores 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ((frame_rate_hz_max-frame_rate_hz_min)/frame_rate_hz_step) ~= fix((frame_rate_hz_max-frame_rate_hz_min)/frame_rate_hz_step)
    
    cprintf('err','ATTENTION - lower and upper borders of frame rate range do not match with frame rate steps (see VARIABLE) \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% BASIC CODE FROM CBU_MODEL
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% within Segment A of CBU_LOOP the core elements of the main CBU_MODEL are embedded (Segment II to XI, see below) for calculation of model indices (cbu index, non-cbu index, reference index) under variable conditions which are altered by the loop itself
% it must be noted that there are a view code line alterations within the embedded segments below in comparison to the fundamental code lines within CBU_MODEL, see following explanation:
%   Segment II: all newly defined variables within the CBU_LOOP are deactivated within the CBU_MODEL based code lines (frame_rate_hz, subframe_number, duty_cycle, content_width_px, content_velo_hor_pxframe, eye_velo_hor_pxframe)
%   Segment III: frame_number is floored / rounded down since (1) the entered model parameter can lead to decimal numbers for the calculated frame_number in some cases and (2) the decimal numbered frame_number needs to be rounded down since it must be clear how many full frame cycles to consider during the following steps of calculation (within CBU_MODEL only integer numbers for frame_number are possile because of the limitation of data entry for the underlying parameter like content or eye movement velocity for example)
%   Segment IV: data entry rules allow the input parameters content_velo_hor_pxframe and eye_velo_hor_pxframe to be decimal numbers since the code lines for the control of the entry of only positive integers including zero are deactivated (within CBU_MODEL the input of content_velo_hor_pxframe and eye_velo_hor_pxframe needs to be a positive integer including zero)
%   Segment XI: frame_number_not_rounded is activated in order to multiply the single frame cycle based model output (option 2) with the exact decimal number of frames cycles presented within the defined sequence; this alteration is necessary to make different frame rates comparable, otherwise different ON times of light stimulation within one single frame cycle would lead to clear falsification of the model output for different frame rates to be compared (within CBU_MODEL the single frame cycle based model output is multiplied with frame_number which is always a positive integer as already mentioned earlier)
% the listed alterations of the embedded code lines from CBU_MODEL do not lead to different model outputs when comparing CBU_LOOP and CBU_MODEL by determining the same input parameter, alterations are only made to allow running the CBU_LOOP for a wider range of input parameter

for loop_count = 1:1:size(loop_table,1)
     
frame_rate_hz = loop_table(loop_count,4); % specific frame rate for the execution of actual loop

eye_velo_hor_pxframe = loop_table(loop_count,5); % specific eye movement velocity in [PX/FR] for the execution of actual loop (adapted in order to guarantee a stable eye movement velocity in [DEG/SEC] for various frame rates)

content_velo_hor_pxframe = loop_table(loop_count,6); % specific content movement velocity in [PX/FR] for the execution of actual loop (adapted in order to guarantee a stable content movement velocity in [DEG/SEC] for various frame rates)

clearvars -except subframe_number duty_cycle content_width_px frame_rate_hz frame_rate_hz_min frame_rate_hz_max frame_rate_hz_step eye_velo_hor_pxframe content_velo_hor_pxframe output_table loop_table loop_count % deletes all parameter (except the mentioned) in order to not take previous runs into account

%% SEGMENT II - DATA ENTRY MASK

% ANATOMICAL PARAMETER OF STANDARD OBSERVER EYE REFERING TO GULLSTRAND-EMSLEY SCHEMATIC HUMAN EYE MODEL

CV_N_relax_m = 0.00706; % corneal vertex (CV) to first nodal point (N) - relaxed
CV_N2_relax_m = 0.00736; % corneal vertex (CV) to second nodal point (N´) - relaxed
CV_F2_relax_m = 0.02389; % corneal vertex (CV) to second principal focus (F´) - relaxed

CV_N_acc_m = 0.00656; % corneal vertex (CV) to first nodal point (N) - accommodated (10.9D)
CV_N2_acc_m = 0.00691; % corneal vertex (CV) to second nodal point (N´) - accommodated (10.9D)
CV_F2_acc_m = 0.02125; % corneal vertex (CV) to second principal focus (F´) - accommodated (10.9D)

CV_Z_m = 0.01350; % cornela vertex (CV) to eyes center of rotation (Z), not based on Gullstrand-Emsley Eye Model

RR_m = 0.01200; % radius of retina (RR), simplified assumption that schematic retina has spherical shape of r = 12mm (Smith & Atchison, 1997, S.678), new findings show that emmetropic eyes tend strongly to oblate retinal shape (Atchison et al., 2005) 

% COGNITIVE PARAMETER REGARDING TEMPORAL AND SPATIAL SUMMATION OF VISUAL STIMULI

bloch_time_cones_sec = 0.0;
bloch_time_rods_sec = 0.0;

ricco_space_mm = 0.0;

% VIEWING PARAMETER REGARDING OBSERVER POSITION TO DISPLAY UNIT

viewing_distance_m = 0.571690; % distance from display level to corneal vertex of observer eye

% HARDWARE PARAMETER REGARDING DISPLAY UNIT

% frame_rate_hz = 90.0; % number of frames displayed per second
frame_duration_sec = 1./frame_rate_hz; % duration of one single frame cycle

% subframe_number = 3; % number of subframes displayed within on single frame cycle, specification of subframes (light intensity & color, see below) is designed for up to 6 subframes (could be extended if necessary)
subframe_rate_hz = frame_rate_hz.*subframe_number; % duration of one single subframe cycle, x-times (depending on number of subframes) the frame rate of the display
subframe_duration_sec = 1./subframe_rate_hz; % duration of one single subframe (on and off time included)

% duty_cycle = 0.30; % on-off-ratio of one single subframe cycle (on-time divided by complete subframe-time)

subframe_on_sec = (1./subframe_rate_hz).*duty_cycle; % time subframe is switched on
subframe_off_sec = (1./subframe_rate_hz).*(1-duty_cycle); % time subframe is switched off

resolution_hor_px = 912; resolution_ver_px = 1140; % count of pixels in horizontal / vertical direction

display_aspect_ratio_hor = 16; display_aspect_ratio_ver = 10; % ratio between horizontal and vertical dimension of display

display_dia_inch = 21.386496; % diagonal dimension of display
    display_dia_m = display_dia_inch.*(25.4./1000);

display_hor_inch = (cos(atan(display_aspect_ratio_ver./display_aspect_ratio_hor))).*display_dia_inch; % horizontal dimension of display
    display_hor_m = display_hor_inch.*(25.4./1000);

display_ver_inch = (sin(atan(display_aspect_ratio_ver./display_aspect_ratio_hor))).*display_dia_inch; % vertical dimension of display
    display_ver_m = display_ver_inch.*(25.4./1000);
       
pixel_pitch_hor_m = ((cos(atan(display_aspect_ratio_ver/display_aspect_ratio_hor))).*display_dia_m)./resolution_hor_px; % distance of two horizontal neighbour pixels measured from both centers of the pixel
    pixel_pitch_hor_deg = abs(((atand((pixel_pitch_hor_m./2)./(viewing_distance_m+CV_N_relax_m))))-((atand((-pixel_pitch_hor_m./2)./(viewing_distance_m+CV_N_relax_m)))));
pixel_pitch_ver_m = ((sin(atan(display_aspect_ratio_ver/display_aspect_ratio_hor))).*display_dia_m)./resolution_ver_px; % distance of two vertical neighbour pixels measured from both centers of the pixel
    pixel_pitch_ver_deg = abs(((atand((pixel_pitch_ver_m./2)./(viewing_distance_m+CV_N_relax_m))))-((atand((-pixel_pitch_ver_m./2)./(viewing_distance_m+CV_N_relax_m)))));

viewing_angle_dia_deg = abs(((atand((display_dia_m./2)./(viewing_distance_m+CV_N_relax_m))))-((atand((-display_dia_m./2)./(viewing_distance_m+CV_N_relax_m))))); % angle over which the diagonal dimension of the display is seen by the observer
viewing_angle_hor_deg =  abs(((atand((display_hor_m./2)./(viewing_distance_m+CV_N_relax_m))))-((atand((-display_hor_m./2)./(viewing_distance_m+CV_N_relax_m))))); % angle over which the horizontal dimension of the display is seen by the observer 
viewing_angle_ver_deg =  abs(((atand((display_ver_m./2)./(viewing_distance_m+CV_N_relax_m))))-((atand((-display_ver_m./2)./(viewing_distance_m+CV_N_relax_m))))); % angle over which the vertical dimension of the display is seen by the observer
    
% LIGHTING PARAMETER (INTENSITY/COLOR) REGARDING THE DISPLAYING UNIT (SUBFRAME NUMBER LIMITED TO MAX. 6SF)

% subframe radiance (Le) in [W/(SR.*M2)] as the physical quantity related to light intensity as the power of light

sf1_radiance = 0.122407; % red subframe, see chromaticity coordinates (value determined from led_current_values_112018 measurement w/o correction of 360HZ measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file)
sf2_radiance = 0.245612; % green subframe, see chromaticity coordinates (value determined from led_current_values_112018 measurement w/o correction of 360HZ measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file)
sf3_radiance = 0.157638; % blue subframe, see chromaticity coordinates (value determined from led_current_values_112018 measurement w/o correction of 360HZ measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file)
sf4_radiance = 0.000000; % blank subframe (reserved for extended subframe calculation)
sf5_radiance = 0.000000; % blank subframe (reserved for extended subframe calculation)
sf6_radiance = 0.000000; % blank subframe (reserved for extended subframe calculation)
 
sf_radiance = [sf1_radiance; sf2_radiance; sf3_radiance; sf4_radiance; sf5_radiance; sf6_radiance]; % summarizing table
control = sf_radiance~=0; a = sum(control,2); b = zeros(6,1); b(1:subframe_number,1) = 1; if ~(sum(a==b)==6); error (['ATTENTION - control radiance specification of all ' num2str(subframe_number) ' subframes to be calculated']); end % control if the the defined subframe number (see above) is specified
sf_radiance( all(~sf_radiance,2), : ) = []; % strikes the zeros at the end of the summarizing table in case of not applying all six possible SF

% subframe luminance (Lv) in [CD/M2] as the photometric quantity related to the perceived light intensity by the human eye

sf1_luminance = 27.3882; % red subframe, see chromaticity coordinates (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file) 
sf2_luminance = 126.0303; % green subframe, see chromaticity coordinates below (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file) 
sf3_luminance = 5.1380; % blue subframe, see chromaticity coordinates below (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file)  
sf4_luminance = 0.0000; % blank subframe (reserved for extended subframe calculation)
sf5_luminance = 0.0000; % blank subframe (reserved for extended subframe calculation)
sf6_luminance = 0.0000; % blank subframe (reserved for extended subframe calculation)

sf_luminance = [sf1_luminance; sf2_luminance; sf3_luminance; sf4_luminance; sf5_luminance; sf6_luminance]; % summarizing table
control = sf_luminance~=0; a = sum(control,2); b = zeros(6,1); b(1:subframe_number,1) = 1; if ~(sum(a==b)==6); error (['ATTENTION - control luminance specification of all ' num2str(subframe_number) ' subframes to be calculated']); end % control if the the defined subframe number (see above) is specified
sf_luminance( all(~sf_luminance,2), : ) = []; % strikes the zeros at the end of the summarizing table in case of not applying all six possible SF

% subframe chromaticity coordinate (x) for light color specification

sf1_xcoordinate = 0.6912; % red subframe (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file)
sf2_xcoordinate = 0.3018; % green subframe (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file) 
sf3_xcoordinate = 0.1527; % blue subframe (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file) 
sf4_xcoordinate = 0.0000; % blank subframe (reserved for extended subframe calculation)
sf5_xcoordinate = 0.0000; % blank subframe (reserved for extended subframe calculation)
sf6_xcoordinate = 0.0000; % blank subframe (reserved for extended subframe calculation)

sf_xcoordinate = [sf1_xcoordinate; sf2_xcoordinate; sf3_xcoordinate; sf4_xcoordinate; sf5_xcoordinate; sf6_xcoordinate]; % summarizing table
control = sf_xcoordinate~=0; a = sum(control,2); b = zeros(6,1); b(1:subframe_number,1) = 1; if ~(sum(a==b)==6); error (['control chromaticity coordinates (x) of all ' num2str(subframe_number) ' subframes to be calculated']); end % control if the the defined subframe number (see above) is specified
sf_xcoordinate( all(~sf_xcoordinate,2), : ) = []; % strikes the zeros at the end of the summarizing table in case of not applying all six possible SF

% subframe chromaticity coordinate (y) for light color specification

sf1_ycoordinate = 0.3078; % red subframe (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file) 
sf2_ycoordinate = 0.6276; % green subframe (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file) 
sf3_ycoordinate = 0.0240; % blue subframe (value determined from led_current_values_11&122018 measurement, mean value of various FR and BFR at maximum light intensity of projector, settings of measurement described in xls.file)
sf4_ycoordinate = 0.0000; % blank subframe (reserved for extended subframe calculation)
sf5_ycoordinate = 0.0000; % blank subframe (reserved for extended subframe calculation)
sf6_ycoordinate = 0.0000; % blank subframe (reserved for extended subframe calculation)

sf_ycoordinate = [sf1_ycoordinate; sf2_ycoordinate; sf3_ycoordinate; sf4_ycoordinate; sf5_ycoordinate; sf6_ycoordinate]; % summarizing table
control = sf_ycoordinate~=0; a = sum(control,2); b = zeros(6,1); b(1:subframe_number,1) = 1; if ~(sum(a==b)==6); error (['ATTENTION - control chromaticity coordinates (y) of all ' num2str(subframe_number) ' subframes to be calculated']); end % control if the the defined subframe number (see above) is specified
sf_ycoordinate( all(~sf_ycoordinate,2), : ) = []; % strikes the zeros at the end of the summarizing table in case of not applying all six possible SF

% determination of MIXED COLOR specification (from blending of defined single SF specification, see specification of single SF above) 
% for determination of the mixed color a transition of chromaticity coordinates (xyY) to tristimulus values (XYZ) back to chromaticity coordinates (xyY) is necessary (see calculation below) 
% xyY for single SFs (defined above) >> XYZ for single SFs >> XYZ for mixed color (via simple addition of single XYZ values) >> xyY of mixed color 

sf_Xtristimulus = sf_xcoordinate.*(sf_luminance./sf_ycoordinate); % calculation of tristimulus values (X,Y,Z) of single SFs from chromaticity coordinates (x,y) and luminance (Lv or Y) of single SFs
sf_Ytristimulus = sf_luminance;
sf_Ztristimulus = (sf_luminance./sf_ycoordinate).*(1-sf_xcoordinate-sf_ycoordinate);

mix_Xtristimulus = sum(sf_Xtristimulus); % calculation of tristimulus values (X,Y,Z) of mixed color from tristimulus values (X,Y,Z) of single SFs (by simply adding the values of the single SF)
mix_Ytristimulus = sum(sf_Ytristimulus);
mix_Ztristimulus = sum(sf_Ztristimulus);

sf_xcoordinate_mix = mix_Xtristimulus./(mix_Xtristimulus+mix_Ytristimulus+mix_Ztristimulus); % calculation of chromaticity coordinates (x,y,Y) of the resulting mixed color from the tristimulus values (X,Y,Z)
sf_ycoordinate_mix = mix_Ytristimulus./(mix_Xtristimulus+mix_Ytristimulus+mix_Ztristimulus); 
sf_luminance_mix = mix_Ytristimulus;

% CONTENT PARAMETER REGARDING SIZE, POSITIONING AND VELOCITY

content_start_hor_px = -236.0; % position of left edge of the content at time 0, point of origin is the horizontal center of the display (left - negative, right - positive), start position is defined in [PX] because it reflects the process of content creation better than defining the content position in [DEG], furthermore non-linearity of positioning in [DEG] from center to periphery is bypassed  
    content_start_hor_deg = (atand((content_start_hor_px.*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m))); % start position of horizontal content movement in [DEG] (refering to first horizontal content edge)	

content_stop_hor_px = 196.0; % position of left edge of the content after movement, point of origin is the horizontal center of the display (left - negative, right - positive), start position is defined in [PX] because it reflects the process of content creation better than defining the content position in [DEG], furthermore non-linearity of positioning in [DEG] from center to periphery is bypassed
    content_stop_hor_deg = (atand((content_stop_hor_px.*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m))); % start position of horizontal content movement in [DEG] (refering to first horizontal content edge)

content_distance_hor_px = content_stop_hor_px - content_start_hor_px; % distance of content movement in horizontal direction
    content_distance_hor_deg = content_stop_hor_deg - content_start_hor_deg; % distance of content movement path [DEG] in horizontal direction  

% content_width_px = 40.0; % horizontal dimension of content, content is assumed to have a simple quadratic or rectangular shape (e.g. vertical bar)

% content_velo_hor_pxframe = 8.0; % content movement velocity [PX/FR] in horizontal direction 

content_position_ver_px = 40.0; % stable position of upper edge [PX] of the content during execution of the movement pattern, point of origin is the vertical center of the display (down - negative, up - positive), vertical content position is defined in [PX] because it reflects the process of content creation better than defining the content position in [DEG], furthermore non-linearity of positioning in [DEG] from center to periphery is bypassed
    content_position_ver_deg = (atand((content_position_ver_px.*pixel_pitch_ver_m)./(viewing_distance_m+CV_N_relax_m))); % position of upper vertical content edge in [DEG] (stable over time since no eye movement in vertical direction)

content_height_px = 80.0; % vertical dimension of content, content is assumed to have a simple quadratic or rectangular shape (e.g. vertical bar)     
    content_height_deg =  abs((atand((content_position_ver_px.*pixel_pitch_ver_m)./(viewing_distance_m+CV_N_relax_m)))-(atand(((content_position_ver_px-content_height_px).*pixel_pitch_ver_m)./(viewing_distance_m+CV_N_relax_m)))); % stable content height in [DEG] (refering to distance between upper and lower content edge)
   
% EYE MOVEMENT PARAMETER REGARDING POSITIONING AND VELOCITY

eye_start_hor_px = -216.0; % start point of horizontal eye movement, point of origin is the horizontal center of the display (left - negative, right - positive), definition of the used unit [PX] for start position of eye movement adapted from start position of content movement in [PX] in order to keep both definitions comparable (see reasons above)
    eye_start_hor_deg = atand((eye_start_hor_px.*pixel_pitch_hor_m)/(viewing_distance_m+CV_N_relax_m));	% start position of horizontal eye movement in [DEG]

eye_stop_hor_px = 216.0; % stop point of horizontal eye movement, point of origin is the horizontal center of the display (left - negative, right - positive), definition of the used unit [PX] for start position of eye movement adapted from start position of content movement in [PX] in order to keep both definitions comparable (see reasons above)
    eye_stop_hor_deg = atand((eye_stop_hor_px.*pixel_pitch_hor_m)/(viewing_distance_m+CV_N_relax_m)); % stop position of horizontal eye movement in [DEG]

eye_distance_hor_px = eye_stop_hor_px - eye_start_hor_px; % distance eye moves in horizontal direction
    eye_distance_hor_deg = eye_stop_hor_deg - eye_start_hor_deg; % distance of eye movement path [DEG] in horizontal direction 

% eye_velo_hor_pxframe = 8.0; % eye movement velocity [PX/FR] in horizontal direction    

eye_position_ver_px = 0.0; % stable vertical position of eye during movement pattern (in horizontal direction, see below), point of origin is the vertical center of the display (down - negative, up - positive), definition of the used unit [PX] for vertical position of eye movement adapted from vertical position of content in [PX] in order to keep both definitions comparable (see reasons above) 
    eye_position_ver_deg = atand((eye_position_ver_px.*pixel_pitch_ver_m)/(viewing_distance_m+CV_N_relax_m)); % stable vertical eye position in [DEG]

%% SEGMENT III - EXPERIMENTAL CATEGORY (AND NUMBER OF FRAME RATES THAT NEED TO BE CONSIDERED FOR CBU CALCULATION)
% Experimental Category 0 (CAT0) = no content movement and no eye movement existent
% Experimental Category 1 (CAT1) = content movement but no eye movement existent
% Experimental Category 2 (CAT2) = eye movement but no content movement existent
% Experimental Category 3 (CAT3) = content movement as well as eye movement existent
% 
% frame_number = number of frames that have to be considered for CBU calculation for the defined movement pattern of content and eye

if content_velo_hor_pxframe == 0 && eye_velo_hor_pxframe == 0
    cat = 0;
    error('ATTENTION - experimental category 0 (eye and content movement velocity equals 0.0) is not considered for CBU calculation since category 0 necessarily leads to a classification within CBU phase 0 in which no CBU perception can be expected');

elseif content_velo_hor_pxframe > 0 && eye_velo_hor_pxframe == 0
    cat = 1;
    error('ATTENTION - experimental category 1 (eye movement velocity equals 0.0) is not considered for CBU calculation since category 1 necessarily leads to a classification within CBU phase 0 in which no CBU perception can be expected');

elseif content_velo_hor_pxframe == 0 && eye_velo_hor_pxframe > 0
    cat = 2;
    frame_number = floor(eye_distance_hor_px./eye_velo_hor_pxframe+1); % frame number within CAT2 is calculated on basis of the time eye movement is executed (since no content movement is executed at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included  

elseif content_velo_hor_pxframe > 0 && eye_velo_hor_pxframe > 0
    cat = 3;
    
    if (content_distance_hor_px/content_velo_hor_pxframe) >= (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is equal or longer than the time of eye movement (eye_distance_px/eye_velo_pxframe), than the frame number to be considered is calculated on basis of the time eye movement is executed, even when content movement is executed for a longer period of time eye movement time is still basis for calculation since CBU only occurs when eye movement is in progress, otherwise it will be CAT1 (eye fix, content variable) that does not provoke CBU at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included, statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning) 
        frame_number = floor(eye_distance_hor_px./eye_velo_hor_pxframe+1);
    elseif (content_distance_hor_px/content_velo_hor_pxframe) < (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is shorter than the time of eye movement (eye_distance_px/eye_velo_pxframe) than the frame number is calculated on basis of the time content movement is executed since the content disappears directly after reaching the stop point of its movement path (without presentation of content no CBU can occur); calculation of time in [FRAME] by the term (content_distance_px./content_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning)
        frame_number = floor(content_distance_hor_px./content_velo_hor_pxframe+1);
    end

end

%% SEGMENT IV - DATA ENTRY RULES 
% guarantees robustness against data entry error by exclusion of wrongly entered data

if ~isscalar(frame_number) || frame_number < 1 || frame_number ~= fix(frame_number)
    error('ATTENTION - frame_number needs to be a positive integer');    

elseif ~isscalar(subframe_number) || subframe_number < 2 || subframe_number > 6 || subframe_number ~= fix(subframe_number)
    error('ATTENTION - subframe_number needs to be a positive integer from 2 to 6');

elseif ~isscalar(frame_rate_hz) || frame_rate_hz < 20 || frame_rate_hz ~= fix(frame_rate_hz)
    error('ATTENTION - frame_rate_hz needs to be a positive integer equal or greater than 20HZ >> undershooting this frame rate threshold will lead to negative side effects like flickering since the temporal summation of all subframes within one frame can not be guaranteed');

elseif ~isscalar(duty_cycle) || duty_cycle < 0 || duty_cycle > 1
    error('ATTENTION - duty_cycle needs to be a decimal number ranging from 0 to 1');    
    
elseif ~isscalar(resolution_hor_px) || resolution_hor_px < 1 || resolution_hor_px ~= fix(resolution_hor_px)
    error('ATTENTION - resolution_hor_px needs to be a positive integer');

elseif ~isscalar(resolution_ver_px) || resolution_ver_px < 1 || resolution_ver_px ~= fix(resolution_ver_px)
    error('ATTENTION - resolution_ver_px needs to be a positive integer');    

elseif ~isscalar(display_aspect_ratio_hor) || display_aspect_ratio_hor < 1 || display_aspect_ratio_hor ~= fix(display_aspect_ratio_hor)
    error('ATTENTION - display_format_hor needs to be a positive integer');    
    
elseif ~isscalar(display_aspect_ratio_ver) || display_aspect_ratio_ver < 1 || display_aspect_ratio_ver ~= fix(display_aspect_ratio_ver)
    error('ATTENTION - display_format_ver needs to be a positive integer');

elseif ~isscalar(display_dia_inch) || display_dia_inch < 0
    error('ATTENTION - display_dia_inch needs to be a positive decimal number');
    
elseif ~isscalar(viewing_distance_m) || viewing_distance_m < 0
    error('ATTENTION - viewing_distance_m needs to be a positive decimal number');

elseif ~isscalar(content_height_px) || content_height_px < 1 || content_height_px ~= fix(content_height_px) || content_height_px > resolution_ver_px
    error('ATTENTION - content_height_px needs to be a positive integer and can not exceed display dimensions');    

elseif ~isscalar(content_position_ver_px) || content_position_ver_px ~= fix(content_position_ver_px) || abs(content_position_ver_px-content_height_px) > resolution_ver_px./2 || abs(content_position_ver_px) > resolution_ver_px./2
    error('ATTENTION - content_position_ver_px needs to be a integer and can not exceed display dimensions (check content_position_ver_px and content_height_px)');    
    
elseif ~isscalar(content_width_px) || content_width_px < 1 || content_width_px ~= fix(content_width_px) || content_width_px > resolution_hor_px
    error('ATTENTION - content_width_px needs to be a positive integer and can not exceed display dimensions');

elseif ~isscalar(content_start_hor_px) || content_start_hor_px ~= fix(content_start_hor_px) || abs(content_start_hor_px+content_width_px) > resolution_hor_px./2  || abs(content_start_hor_px) > resolution_hor_px./2
    error('ATTENTION - content_start_hor_px needs to be a integer and can not exceed display dimensions (check content_start_hor_px and content_width_px)');

elseif ~isscalar(content_stop_hor_px) || content_stop_hor_px ~= fix(content_stop_hor_px) || abs(content_stop_hor_px+content_width_px) > resolution_hor_px./2 || abs(content_stop_hor_px) > resolution_hor_px./2
    error('ATTENTION - content_stop_hor_px needs to be a integer and can not exceed display dimensions (check content_stop_hor_px and content_width_px)');    
    
% elseif ~isscalar(content_velo_hor_pxframe) || content_velo_hor_pxframe < 0 || content_velo_hor_pxframe ~= fix(content_velo_hor_pxframe)
%     error('ATTENTION - content_velo_hor_pxframe needs to be a positive integer (including zero)');

elseif ~isscalar(eye_position_ver_px) || eye_position_ver_px ~= fix(eye_position_ver_px) || abs(eye_position_ver_px) > resolution_ver_px./2
    error('ATTENTION - eye_position_ver_px needs to be a integer and can not exceed display dimensions');    
    
elseif ~isscalar(eye_start_hor_px) || eye_start_hor_px ~= fix(eye_start_hor_px) || abs(eye_start_hor_px) > resolution_hor_px./2
    error('ATTENTION - eye_start_hor_px needs to be a integer and can not exceed display dimensions');

elseif ~isscalar(eye_stop_hor_px) || eye_stop_hor_px ~= fix(eye_stop_hor_px) || abs(eye_stop_hor_px) > resolution_hor_px./2
    error('ATTENTION - eye_stop_hor_px needs to be a integer and can not exceed display dimensions');    
    
% elseif ~isscalar(eye_velo_hor_pxframe) || eye_velo_hor_pxframe < 0 || eye_velo_hor_pxframe ~= fix(eye_velo_hor_pxframe)
%    error('ATTENTION - eye_velo_hor_pxframe needs to be a positive integer (including zero)');
       
elseif ~isscalar(sf1_radiance) || sf1_radiance < 0
    error('ATTENTION - sf1_radiance needs to be a positive decimal number');
    
elseif ~isscalar(sf2_radiance) || sf2_radiance < 0
    error('ATTENTION - sf2_radiance needs to be a positive decimal number');
    
elseif ~isscalar(sf3_radiance) || sf3_radiance < 0
    error('ATTENTION - sf3_radiance needs to be a positive decimal number');
    
elseif ~isscalar(sf4_radiance) || sf4_radiance < 0
    error('ATTENTION - sf4_radiance needs to be a positive decimal number');
    
elseif ~isscalar(sf5_radiance) || sf5_radiance < 0
    error('ATTENTION - sf5_radiance needs to be a positive decimal number');
    
elseif ~isscalar(sf6_radiance) || sf6_radiance < 0
    error('ATTENTION - sf6_radiance needs to be a positive decimal number');
    
elseif ~isscalar(sf1_luminance) || sf1_luminance < 0
    error('ATTENTION - sf1_luminance needs to be a positive decimal number');
    
elseif ~isscalar(sf2_luminance) || sf2_luminance < 0
    error('ATTENTION - sf2_luminance needs to be a positive decimal number');
    
elseif ~isscalar(sf3_luminance) || sf3_luminance < 0
    error('ATTENTION - sf3_luminance needs to be a positive decimal number');
    
elseif ~isscalar(sf4_luminance) || sf4_luminance < 0
    error('ATTENTION - sf4_luminance needs to be a positive decimal number');
    
elseif ~isscalar(sf5_luminance) || sf5_luminance < 0
    error('ATTENTION - sf5_luminance needs to be a positive decimal number');
    
elseif ~isscalar(sf6_luminance) || sf6_luminance < 0
    error('ATTENTION - sf6_luminance needs to be a positive decimal number');
            
elseif ~isscalar(sf1_xcoordinate) || sf1_xcoordinate < 0
    error('ATTENTION - sf1_xcoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf2_xcoordinate) || sf2_xcoordinate < 0
    error('ATTENTION - sf2_xcoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf3_xcoordinate) || sf3_xcoordinate < 0
    error('ATTENTION - sf3_xcoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf4_xcoordinate) || sf4_xcoordinate < 0
    error('ATTENTION - sf4_xcoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf5_xcoordinate) || sf5_xcoordinate < 0
    error('ATTENTION - sf5_xcoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf6_xcoordinate) || sf6_xcoordinate < 0
    error('ATTENTION - sf6_xcoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf1_ycoordinate) || sf1_ycoordinate < 0
    error('ATTENTION - sf1_ycoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf2_ycoordinate) || sf2_ycoordinate < 0
    error('ATTENTION - sf2_ycoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf3_ycoordinate) || sf3_ycoordinate < 0
    error('ATTENTION - sf3_ycoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf4_ycoordinate) || sf4_ycoordinate < 0
    error('ATTENTION - sf4_ycoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf5_ycoordinate) || sf5_ycoordinate < 0
    error('ATTENTION - sf5_ycoordinate needs to be a positive decimal number');
    
elseif ~isscalar(sf6_ycoordinate) || sf6_ycoordinate < 0
    error('ATTENTION - sf6_ycoordinate needs to be a positive decimal number');    
    
elseif ~isscalar(CV_N_relax_m) || CV_N_relax_m < 0
    error('ATTENTION - CV_N_relax_m needs to be a positive decimal number');    

elseif ~isscalar(CV_N2_relax_m) || CV_N2_relax_m < 0
    error('ATTENTION - CV_N2_relax_m needs to be a positive decimal number');        

elseif ~isscalar(CV_F2_relax_m) || CV_F2_relax_m < 0
    error('ATTENTION - CV_F2_relax_m needs to be a positive decimal number');        

elseif ~isscalar(CV_N_acc_m) || CV_N_acc_m < 0
    error('ATTENTION - CV_N_acc_m needs to be a positive decimal number');        
 
elseif ~isscalar(CV_N2_acc_m) || CV_N2_acc_m < 0
    error('ATTENTION - CV_N2_acc_m needs to be a positive decimal number');        

elseif ~isscalar(CV_F2_acc_m) || CV_F2_acc_m < 0
    error('ATTENTION - CV_F2_acc_m needs to be a positive decimal number');        
    
elseif ~isscalar(CV_Z_m) || CV_Z_m < 0
    error('ATTENTION - CV_Z_m needs to be a positive decimal number');        
    
elseif ~isscalar(RR_m) || RR_m < 0
    error('ATTENTION - RR_m needs to be a positive decimal number');
        
elseif ~isscalar(bloch_time_cones_sec) || bloch_time_cones_sec < 0
    error('ATTENTION - bloch_time_cones_sec needs to be a positive decimal number');
    
elseif ~isscalar(bloch_time_rods_sec) || bloch_time_rods_sec < 0
    error('ATTENTION - bloch_time_rods_sec needs to be a positive decimal number');    
    
elseif ~isscalar(ricco_space_mm) || ricco_space_mm < 0
    error('ATTENTION - ricco_space_mm needs to be a positive decimal number');
    
end

%% SEGMENT V - SUBFRAME CHARACTERISTICS (CALCULATION OF POSITION AND TIMING AS WELL AS INTENSITY AND COLOR)
%  calculation of retinal positioning (horizontal, vertical) and presentation timing as well as intensity of luminance and color of all presented subframes (f1_sf1, f1_sf2, f1_sf3, f2_sf1, f2_sf2, f2_sf3, f3_sf1, f3_sf2, f3_sf3 ...)

% SUBFRAME_LOCATION_TABLE_HOR_X1 shows horizontal retinal position of subframes (four basic points of rhomboid --> A´s, C´s, E´s, G´s) in [DEG] and is arranged as follows ...
% column 1 = frame 1, column 2 = frame 2, column 3 = frame 3, ...
% row 1 - 4 = subframe 1 (row 1 = point A, row 2 = point C, row 3 = point E, row 4 = point G), row 5 - 8 = subframe 2, row 9 - 12 = subframe 3, ...    

    subframe_location_table_hor_x1 = zeros(subframe_number.*4,frame_number);
    
    % transfer stable [PX]/[PX/FR] values regarding width, movement and velocity of content respectively movement and velocity of eye to the corresponding variable display position dependent [DEG]/[DEG/SEC] values  (before starting to calculate basic points A, C, E and G)
    % ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    content_width_variable_deg_table = zeros(frame_number,1); % table of all content widths in [DEG] for all frame cycles during the specified sequence (variable since display unit is assumed to be a flat screen)
    content_velo_hor_variable_degs_table = zeros(frame_number,1); % table of all content movement velocities in [DEG/SEC] for all frame cycles during the specified sequence (variable since display unit is assumed to be a flat screen)
    eye_velo_hor_variable_degs_table = zeros(frame_number,1); % table of all eye movement velocities in [DEG/SEC] for all frame cycles during the specified sequence (variable since display unit is assumed to be a flat screen)
    
    for ii = 1:frame_number
                              
        % content width
        
        eye_position_hor_variable_deg = atand(((eye_start_hor_px+(ii-1).*eye_velo_hor_pxframe).*pixel_pitch_hor_m)/(viewing_distance_m+CV_N_relax_m)); % variable horizontal eye position [DEG] in dependency of point in time during the specified sequence
        
        content_position_hor_edge1_variable_deg = (atand(((content_start_hor_px+(ii-1).*content_velo_hor_pxframe).*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m))); % variable position of first horizontal content edge [DEG] in dependency of point in time during the specified sequence
        content_position_hor_edge2_variable_deg = (atand(((content_start_hor_px+(ii-1).*content_velo_hor_pxframe+content_width_px).*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m))); % variable position of second horizontal content edge [DEG] in dependency of point in time during the specified sequence
        
        content_width_variable_deg = abs((content_position_hor_edge1_variable_deg-eye_position_hor_variable_deg)-(content_position_hor_edge2_variable_deg-eye_position_hor_variable_deg)); % variable content width [DEG] in dependency of point in time (refering to distance between first and second content edge) during the specified sequence
        
        content_width_variable_deg_table(ii,1) = content_width_variable_deg; % table of all content widths in [DEG] for all frame cycles during the specified sequence (variable since display unit is assumed to be a flat screen)
        
        % content movement velocity
        
        content_velo_hor_variable_degs = abs((atand(((content_start_hor_px+(ii-1).*content_velo_hor_pxframe+content_width_px./2).*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m)))-(atand(((content_start_hor_px+(ii).*content_velo_hor_pxframe+content_width_px./2).*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m))))/frame_duration_sec; % variable horizontal content movement velocity [DEG/SEC] in dependency of point in time (only positive values for left to right movement, see general model assumptions) during the specified sequence
        
        content_velo_hor_variable_degs_table(ii,1) = content_velo_hor_variable_degs; % table of all content movement velocities in [DEG/SEC] for all frame cycles during the specified sequence (variable since display unit is assumed to be a flat screen)
        
        % content movement distance
        
        content_move_hor_frame_complete_variable_deg = content_velo_hor_variable_degs.*frame_duration_sec; % variable content movement distance [DEG] in horizontal direction while one complete frame cycle is presented (depending on point in time during the specified sequence)
        
        % eye movement velocity
        
        eye_velo_hor_variable_degs = abs((atand(((eye_start_hor_px+(ii-1).*eye_velo_hor_pxframe).*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m)))-(atand(((eye_start_hor_px+(ii).*eye_velo_hor_pxframe).*pixel_pitch_hor_m)./(viewing_distance_m+CV_N_relax_m))))/frame_duration_sec; % variable horizontal eye movement velocity [DEG/SEC] in dependency of point in time during the specified sequence (only positive values for left to right movement, see general model assumptions)
        
        eye_velo_hor_variable_degs_table(ii,1) = eye_velo_hor_variable_degs; % table of all eye movement velocities in [DEG/SEC] for all frame cycles during the specified sequence (variable since display unit is assumed to be a flat screen)
        
        % eye movement distance
        
        eye_move_hor_frame_complete_variable_deg = eye_velo_hor_variable_degs.*frame_duration_sec; % eye movement distance [DEG] in horizontal direction while one complete frame cycle is presented (depending on point in time during the specified sequence)
        eye_move_hor_subframe_on_variable_deg = eye_velo_hor_variable_degs.*subframe_on_sec; % eye movement distance [DEG] in horizontal direction while subframe is on within one subframe cycle (depending on point in time during the specified sequence)
        eye_move_hor_subframe_off_variable_deg = eye_velo_hor_variable_degs.*subframe_off_sec; % eye movement distance [DEG] in horizontal direction while subframe is off within one subframe cycle (depending on point in time during the specified sequence)
        % ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------             

        for jj = 1:4:(subframe_number.*4) % calculation for all A´s
            if ii == 1 && jj == 1
                subframe_location_table_hor_x1(jj,ii) = content_start_hor_deg + content_width_variable_deg - eye_start_hor_deg;
            elseif ii == 1 && jj > 1
                subframe_location_table_hor_x1(jj,ii) = subframe_location_table_hor_x1((jj-4),ii) - eye_move_hor_subframe_on_variable_deg - eye_move_hor_subframe_off_variable_deg;
            elseif ii > 1 && jj == 1
                subframe_location_table_hor_x1(jj,ii) = subframe_location_table_hor_x1(jj,ii-1) - eye_move_hor_frame_complete_variable_deg + content_move_hor_frame_complete_variable_deg;
            else 
                subframe_location_table_hor_x1(jj,ii) = subframe_location_table_hor_x1(jj-4,ii-1) - eye_move_hor_subframe_on_variable_deg - eye_move_hor_subframe_off_variable_deg - eye_move_hor_frame_complete_variable_deg + content_move_hor_frame_complete_variable_deg;
            end
        end
        
        for jj = 2:4:(subframe_number.*4) % calculation for all C´s
            subframe_location_table_hor_x1(jj,ii) = subframe_location_table_hor_x1(jj-1,ii) - content_width_variable_deg;
        end
        
        for jj = 3:4:(subframe_number.*4) % calculation for all E´s
            subframe_location_table_hor_x1(jj,ii) = subframe_location_table_hor_x1(jj-2,ii) - content_width_variable_deg - eye_move_hor_subframe_on_variable_deg;
        end
        
        for jj = 4:4:(subframe_number.*4) % calculation for all G´s
            subframe_location_table_hor_x1(jj,ii) = subframe_location_table_hor_x1(jj-3,ii) - eye_move_hor_subframe_on_variable_deg;
        end
    end
    
% SUBFRAME_LOCATION_TABLE_HOR_MM >> transition of retinal position specification from angle value in [DEG] (see subframe_location_table_hor_x1 in the section above) to distance value in [MM] (see subframe_location_table_hor_mm below)
% the calculations for retinal positions of the displayed content above are made under the assumption that the optical path goes through the first (N) and second nodal point (N´) of the viewers eye which has the benefit that the angle of the incident light ray is similar to the angle of the outgoing light ray - no further calculation is necessary (see subframe_location_table_hor_x1), this works perfectly well as long as relevant retinal positions can be expressed in [DEG]
% however, if the retinal position has to be stated in [MM] a conversion from the angle in [DEG] to the distance in [MM] for the retinal position is not trivial, the problem is that the retinal position in [DEG] (see subframe_location_table_hor_x1) can not be transfered into retinal position in [MM] via simple calculation of the circular arc by using the the parameter (1) retinal radius (see RR_m) and (2) angle x1 (see subframe_location_table_hor_x1) starting from the second nodal point (N´), the reason is that the midpoint of the circle (Z) representing the retinal shape (approximation!) is not equal to the second nodal point (N´), which is relevant for the transition of the optical path from the display unit through the optics of the eye to the retina, instead angle x1 has to be converted into angle x2 (see code for subframe_location_table_hor_x2 below) at first, afterwards the calculation of the circular arc can be operated by using the parameter (1) retinal radius (see RR_m) and (2) angle x2 (see subframe_location_table_hor_x2) starting from the midpoint of the circle representing the retinal shape, the result is the retinal position as distance to the point of origin of the coordinate system (fovea centralis) in [MM] (see subframe_location_table_hor_mm) >> there is also a graphical illustration of the geometrical relations >> see confluence >> theoretical cbu model >> display-to-retina correction (helpful for a better understanding)
    
    for ii = 1:1:size(subframe_location_table_hor_x1,2)
        
        for jj = 1:1:size(subframe_location_table_hor_x1,1)
            
            if subframe_location_table_hor_x1(jj,ii) < 0
                
                subframe_location_table_hor_x2(jj,ii) = (90-(180-(180-(180-abs(subframe_location_table_hor_x1(jj,ii))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(subframe_location_table_hor_x1(jj,ii))-90))))).*(tan((deg2rad(abs(subframe_location_table_hor_x1(jj,ii))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))))).*(-1);
                subframe_location_table_hor_mm(jj,ii) = ((subframe_location_table_hor_x2(jj,ii)./360).*2.*pi.*RR_m).*1000;
                
            else
                
                subframe_location_table_hor_x2(jj,ii) = 90-(180-(180-(180-abs(subframe_location_table_hor_x1(jj,ii))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(subframe_location_table_hor_x1(jj,ii))-90))))).*(tan((deg2rad(abs(subframe_location_table_hor_x1(jj,ii))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))));
                subframe_location_table_hor_mm(jj,ii) = ((subframe_location_table_hor_x2(jj,ii)./360).*2.*pi.*RR_m).*1000;
            
            end
            
        end
        
    end

% SUBFRAME_LOCATION_TABLE_VER_X1 shows vertical retinal position of subframes in [DEG]
% the calculation of the subframe positions in vertical direction is trivial (in comparison to the horizontal positions of the subframes, see above) since there is no content or eye movement in vertical direction
% therefore, there is only an upper and lower content edge which is identical for all subframes
% the two values are listed in the SUBFRAME_LOCATION_TABLE_VER_X1 (2row vector, row1 = upper content edge = mathematically higher value, row2 = lower content edge = mathematically lower value)

subframe_location_table_ver_x1(1,1) = content_position_ver_deg - eye_position_ver_deg; % stable upper vertical content edge
subframe_location_table_ver_x1(2,1) = content_position_ver_deg - content_height_deg - eye_position_ver_deg; % lower content edge

% SUBFRAME_LOCATION_TABLE_VER_MM >> transition of retinal position specification from angle value in [DEG] (see subframe_location_table_ver_x1 in the section above) to distance value in [MM] (see subframe_location_table_ver_mm below)
% the calculations for retinal positions of the displayed content above are made under the assumption that the optical path goes through the first (N) and second nodal point (N´) of the viewers eye which has the benefit that the angle of the incident light ray is similar to the angle of the outgoing light ray - no further calculation is necessary (see subframe_location_table_ver_x1), this works perfectly well as long as relevant retinal positions can be expressed in [DEG]
% however, if the retinal position has to be stated in [MM] a conversion from the angle in [DEG] to the distance in [MM] for the retinal position is not trivial, the problem is that the retinal position in [DEG] (see subframe_location_table_ver_x1) can not be transfered into retinal position in [MM] via simple calculation of the circular arc by using the the parameter (1) retinal radius (see RR_m) and (2) angle x1 (see subframe_location_table_ver_x1) starting from the second nodal point (N´), the reason is that the midpoint of the circle (Z) representing the retinal shape (approximation!) is not equal to the second nodal point (N´), which is relevant for the transition of the optical path from the display unit through the optics of the eye to the retina, instead angle x1 has to be converted into angle x2 (see code for subframe_location_table_ver_x2 below) at first, afterwards the calculation of the circular arc can be operated by using the parameter (1) retinal radius (see RR_m) and (2) angle x2 (see subframe_location_table_ver_x2) starting from the midpoint of the circle representing the retinal shape, the result is the retinal position as distance to the point of origin of the coordinate system (fovea centralis) in [MM] (see subframe_location_table_ver_mm) >> there is also a graphical illustration of the geometrical relations >> see confluence >> theoretical cbu model >> display-to-retina correction (helpful for a better understanding)

    for ii = 1:1:size(subframe_location_table_ver_x1,2)
        
        for jj = 1:1:size(subframe_location_table_ver_x1,1)
            
            if subframe_location_table_ver_x1(jj,ii) < 0
                
                subframe_location_table_ver_x2(jj,ii) = (90-(180-(180-(180-abs(subframe_location_table_ver_x1(jj,ii))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(subframe_location_table_ver_x1(jj,ii))-90))))).*(tan((deg2rad(abs(subframe_location_table_ver_x1(jj,ii))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))))).*(-1);
                subframe_location_table_ver_mm(jj,ii) = ((subframe_location_table_ver_x2(jj,ii)./360).*2.*pi.*RR_m).*1000;
                
            else
                
                subframe_location_table_ver_x2(jj,ii) = 90-(180-(180-(180-abs(subframe_location_table_ver_x1(jj,ii))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(subframe_location_table_ver_x1(jj,ii))-90))))).*(tan((deg2rad(abs(subframe_location_table_ver_x1(jj,ii))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))));
                subframe_location_table_ver_mm(jj,ii) = ((subframe_location_table_ver_x2(jj,ii)./360).*2.*pi.*RR_m).*1000;
            
            end
            
        end
        
    end

% SUBFRAME_TIME_TABLE shows time of retinal presentation for subframes (four basic points of rhomboid --> A´s, C´s, E´s, G´s) and is arranged as follows ...
% column 1 = frame 1, column 2 = frame 2, column 3 = frame 3, ...
% row 1 - 4 = subframe 1 (row 1 = point A, row 2 = point C, row 3 = point E, row 4 = point G), row 5 - 8 = subframe 2, row 9 - 12 = subframe 3, ...
    
    subframe_time_table = zeros(subframe_number.*4,frame_number);
    
    for ii = 1:frame_number
        
        for jj = 1:4:(subframe_number.*4)  % calculation for all A´s            
            if ii == 1 && jj <= 2
                subframe_time_table(jj,ii) = 0.0;
            elseif ii == 1 && jj > 2
                subframe_time_table(jj,ii) = subframe_time_table((jj-4),ii) + subframe_duration_sec ;
            elseif ii > 1 && jj == 1
                subframe_time_table(jj,ii) = subframe_time_table(jj,ii-1) + frame_duration_sec;
            else 
                subframe_time_table(jj,ii) = subframe_time_table(jj-4,ii-1) + subframe_duration_sec + frame_duration_sec;
            end
        end
       
        for jj = 2:4:(subframe_number.*4) % calculation for all C´s
            subframe_time_table(jj,ii) = subframe_time_table(jj-1,ii);
        end
        
        for jj = 3:4:(subframe_number.*4) % calculation for all E´s
            subframe_time_table(jj,ii) = subframe_time_table(jj-2,ii) + subframe_on_sec;
        end
        
        for jj = 4:4:(subframe_number.*4) % calculation for all G´s
            subframe_time_table(jj,ii) = subframe_time_table(jj-3,ii) + subframe_on_sec;
        end    
    end
   
% SUBFRAME_INTENSITY_TABLE is available in two versions with two different measurands for light intensity of the relevant subframes
% subframe_radiance_intensity_table shows light intensity as physical quantity in radiance [W/(SR.*M2)]
% subframe_luminance_intensity_table shows light intensity as photometric quantity in luminance [CD/M2]
% both tables are arranged identically as follows ...
% column 1 = frame 1, column 2 = frame 2, column 3 = frame 3, ...
% row 1 - 4 = subframe 1 (row 1 = point A, row 2 = point C, row 3 = point E, row 4 = point G), row 5 - 8 = subframe 2, row 9 - 12 = subframe 3, ...    
    
for aa = 1:1:subframe_number    
    subframe_radiance_intensity_table((aa.*4-3):(aa.*4),1:frame_number) = sf_radiance(aa,1); % for the moment radiance is identical for within subframe points A, C, E and G (as well as for same subframes in different frames) 
    subframe_luminance_intensity_table((aa.*4-3):(aa.*4),1:frame_number) = sf_luminance(aa,1); % for the moment luminance is identical for within subframe points A, C, E and G (as well as for same subframes in different frames)
end
    
% SUBFRAME_COLOR_TABLE is available in two versions with two different measurands that are necessary to specify the light color of the relevant subframes
% subframe_xcoordinate_color_table shows the chromaticity coordinate x
% subframe_ycoordinate_color_table shows the chromaticity coordinate y
% both tables are arranged identically as follows ...
% column 1 = frame 1, column 2 = frame 2, column 3 = frame 3, ...
% row 1 - 4 = subframe 1 (row 1 = point A, row 2 = point C, row 3 = point E, row 4 = point G), row 5 - 8 = subframe 2, row 9 - 12 = subframe 3, ...

for aa = 1:1:subframe_number    
    subframe_xcoordinate_color_table((aa.*4-3):(aa.*4),1:frame_number) = sf_xcoordinate(aa,1); % for the moment chromaticity coordinate x is identical for within subframe points A, C, E and G (as well as for same subframes in different frames)
    subframe_ycoordinate_color_table((aa.*4-3):(aa.*4),1:frame_number) = sf_ycoordinate(aa,1); % for the moment chromaticity coordinate y is identical for within subframe points A, C, E and G (as well as for same subframes in different frames)
end

%% SEGMENT VI - COLOR BREAK-UP PHASE    
% Phase 0 (PH0) = no CBU occurs 
% Phase 1 (PH1) = CBU occurs with following color pattern: B/C/W/Y/R
% Phase 2 (PH2) = CBU occurs with following color pattern: B/C/G/Y/R
% Phase 3 (PH3) = CBU occurs with following color pattern: B/GAP/G/GAP/R

eye_velo_hor_variable_degs_table_index = find(eye_velo_hor_variable_degs_table == max(eye_velo_hor_variable_degs_table), 1); % determination of CBU phase explicitly refers to the frame cycle with the highest eye movement velocity in [DEG/SEC] (variable angular values because of the flat screen characteristic of the display unit, frame cycle within which the corresponding content position is the most central regarding its display unit position leads to the higest angular values)

if subframe_location_table_hor_x1(1,eye_velo_hor_variable_degs_table_index) == subframe_location_table_hor_x1(5,eye_velo_hor_variable_degs_table_index) % proof assumption that A(SF1) is equal to A(SF2) for verification of PH0; since a movement pattern with constant retinal velocity in [PX/FR] (velocity values are constant for complete sequence) respectively [DEG/SEC] (at least constant velocity values within a frame cycle) is presumed (see Segment I), it is sufficient to check the assumption for the first two subframes of the chosen frame cycle only; it should generally be noted that (at present) it is not possible to confirm this if loop statement that would result in a classification within PH0 since the classification of the experimental category in the code further above (see Segment III) stops the code from running when CAT0 or CAT1 are determined which are the only categories that result in PH0, for reasons of completeness the if loop refering to PH0 is not excluded from the code;
    cbu_phase = 0;  
    
elseif subframe_location_table_hor_x1(3,eye_velo_hor_variable_degs_table_index) <= subframe_location_table_hor_x1((subframe_number-1)*4+1,eye_velo_hor_variable_degs_table_index) && subframe_location_table_hor_x1(3,eye_velo_hor_variable_degs_table_index) > subframe_location_table_hor_x1((subframe_number-1)*4+3,eye_velo_hor_variable_degs_table_index) % proof assumption that E(SF1) is equal to or smaller than A(last SF) and larger than E(last SF) for verification of PH1 (calculation for the chosen frame cycle);
    cbu_phase = 1;
    
elseif subframe_location_table_hor_x1(3,eye_velo_hor_variable_degs_table_index) > subframe_location_table_hor_x1((subframe_number-1)*4+1,eye_velo_hor_variable_degs_table_index) && subframe_location_table_hor_x1(3,eye_velo_hor_variable_degs_table_index) <= subframe_location_table_hor_x1(5,eye_velo_hor_variable_degs_table_index) % proof assumption that E(SF1) is greater than A(last SF) and equal to or smaller than A(SF2) for verification of PH2 (calculation for the chosen frame cycle);
    cbu_phase = 2;
    
elseif subframe_location_table_hor_x1(3,eye_velo_hor_variable_degs_table_index) > subframe_location_table_hor_x1(5,eye_velo_hor_variable_degs_table_index) % proof assumption that E(SF1) is greater than A(SF2) for verification of PH3 (calculation for the chosen frame cycle);
    cbu_phase = 3;
    
end

%% SEGMENT VII - PHASE TRANSITION THRESHOLDS
% Phase 0 (PH0) = no CBU occurs 
% Phase 1 (PH1) = CBU occurs with following color pattern: B/C/W/Y/R
% Phase 2 (PH2) = CBU occurs with following color pattern: B/C/G/Y/R
% Phase 3 (PH3) = CBU occurs with following color pattern: B/GAP/G/GAP/R

transition_ph1_ph2_pxs = (content_width_px.*subframe_number.*frame_rate_hz)/(2-duty_cycle); % eye movement velocity threshold in [PX/SEC] for transition from PH1 to PH2
transition_ph1_ph2_pxframe = (content_width_px.*subframe_number)/(2-duty_cycle); % eye movement velocity threshold in [PX/FR] for transition from PH1 to PH2
transition_ph1_ph2_degs = (content_width_variable_deg_table(eye_velo_hor_variable_degs_table_index,1).*subframe_number.*frame_rate_hz)/(2-duty_cycle); % eye movement velocity threshold in [DEG/SEC] for transition from PH1 to PH2 (content width in [DEG] within formula refers to the frame cycle with the highest eye movement velocity in [DEG/SEC], variable angular values because of the flat screen characteristic of the display unit)
  
transition_ph2_ph3_pxs = (content_width_px.*subframe_number.*frame_rate_hz)/(1-duty_cycle); % eye movement velocity threshold in [PX/SEC] for transition from PH2 to PH3
transition_ph2_ph3_pxframe = (content_width_px.*subframe_number)/(1-duty_cycle); % eye movement velocity threshold in [PX/FR] for transition from PH2 to PH3
transition_ph2_ph3_degs = (content_width_variable_deg_table(eye_velo_hor_variable_degs_table_index,1).*subframe_number.*frame_rate_hz)/(1-duty_cycle); % eye movement velocity threshold in [DEG/SEC] for transition from PH2 to PH3 (content width in [DEG] within formula refers to the frame cycle with the highest eye movement velocity in [DEG/SEC], variable angular values because of the flat screen characteristic of the display unit) 

%% SEGMENT VIII - INTERSECTIONS AND UNIONS

location_table_short = subframe_location_table_hor_x1(1:2:end,1:1:end); % location_table_short consists only of A´s (start of subframe) and E´s (stop of subframe), structure of table as follows: column 1 = frame 1, column 2 = frame 2, column 3 = frame 3, ... row 1 - 2 = subframe 1 (row 1 = point A, row 2 = point E), row 3 - 4 = subframe 2, row 5 - 6 = subframe 3, ... 

% CBU_TABLE_INTERSECT shows intersection sets (set theory!) for different subframe combinations, the existance of other SF within the calculated range is not excluded (see cbu_assi_one and cbu_assi_tow for separated SF combinations) 
% cbu_table_intersect is filled with the correct data by reading out the location_table_short (via loops, see code) in order to generate an specific arrangement of the cbu_table_interscect structure
%
% general build-up of the cbu_table_intersect (row count, column count):
% -----------------------------------------------------------------------------------------------------------------------
% row count = 2*[((subframe_number*subframe_number)+subframe_number)/2] --> gaussian sum formula (example for 5SF)
% row 1/2 = SF1 (row1 = start point A, row2 = stop point E), row 3/4 = SF2, row 5/6 = SF3, row7/8 = SF4, row9/10 = SF5
% row 11/12 = SF1 i SF2, row 13/14 = SF2 i SF3, row 15/16 = SF3 i SF4, row 17/18 = SF4 i SF5 
% row 19/20 = SF1 i SF2 i SF3, row 21/22 = SF2 i SF3 i SF4, row 23/24 = SF3 i SF4 i SF5  
% row 25/26 = SF1 i SF2 i SF3 i SF4, row 27/28 = SF2 i SF3 i SF4 i SF5
% row 29/30 = SF1 i SF2 i SF3 i SF4 i SF5
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% column_count = frame_number
% column 1 = F1, column 2 = F2, column 3 = F3, column 4 = F4, ...
% -----------------------------------------------------------------------------------------------------------------------

cbu_table_intersect = zeros(2*((subframe_number.*subframe_number+subframe_number)/2), frame_number); % blank cbu_table_intersect filled with zeros at first
cbu_table_intersect(1:1:subframe_number*2,1:1:frame_number) = location_table_short(1:1:subframe_number*2,1:1:frame_number); % upper part of cbu_table_intersect consisting of location_table_short  

% arrangement of the cbu_table_intersect structure (rows) is described with an example of a frame with 4SF included (classic signs of set theory do not work in MATLAB so the following abbreviations are used: i = intersection of sets, u = union of sets, \ = difference of sets):
%
% SF combination                C 1
% -----------------------------------> section 1 (number of subdivisions = subframe_number)
% SF1                           [1]
% SF2                           [2]
% SF3                           [3]
% SF4                           [4]
% -----------------------------------> section 2 (number of subdivisions = subframe_number - 1)
% SF1 i SF2                     [5]
% SF2 i SF3                     [6]
% SF3 i SF4                     [7]
% -----------------------------------> section 3 (number of subdivisions = subframe_number - 2)
% SF1 i SF2 i SF3               [8]
% SF2 i SF3 i SF4               [9]
% -----------------------------------> section 4 (number of subdivisions = subframe_number - 3)
% SF1 i SF2 i SF3 i SF4        [10]
% -----------------------------------
%
% example of cbu_table_intersect structure (rows) for a frame with 5SF included:
%
% SF combination                C 1
% -----------------------------------> section 1 (number of subdivisions = subframe_number)
% SF1                           [1]
% SF2                           [2]
% SF3                           [3]
% SF4                           [4]
% SF5                           [5]
% -----------------------------------> section 2 (number of subdivisions = subframe_number - 1)
% SF1 i SF2                     [6]
% SF2 i SF3                     [7]
% SF3 i SF4                     [8]
% SF4 i SF5                     [9]
% -----------------------------------> section 3 (number of subdivisions = subframe_number - 2)
% SF1 i SF2 i SF3              [10]
% SF2 i SF3 i SF4              [11)
% SF3 i SF4 i SF5              [12]
% -----------------------------------> section 4 (number of subdivisions = subframe_number - 3)
% SF1 i SF2 i SF3 i SF4        [13]
% SF2 i SF3 i SF4 i SF5        [14]
% -----------------------------------> section 5 (number of subdivisions = subframe_number - 4)
% SF1 i SF2 i SF3 i SF4 i SF5  [15]
% -----------------------------------

    for ss = 1:1:subframe_number-1 % for loops (ss-loop, xx-loop, ii-loop ,jj-loop) to fill lower part of cbu_table_intersect by running through all SF scenarios. Resulting order: SF1&SF2 // SF1&SF2&SF3 // SF1&SF2&SF3&SF4 // SF2&SF3 // SF2&SF3&SF4 // SF3&SF4 (exemplary for 4SF)
        
        for xx = 1:1:subframe_number-ss
    
            if xx == 1 % for correct row rearrangement of cbu_table_intersect (see discrepancy between description under CBU_TABLE_INTERSECT and LOOPS, like ss-loop, xx-loop ...) 
                
                yy = 0;
                
            elseif xx == 2
                
                yy = subframe_number-xx;
                
            elseif xx > 2 && (subframe_number - xx) == 1
                
                yy = ((((subframe_number-2).*(subframe_number-2))+(subframe_number-2))/2);
       
            elseif xx > 2 && (subframe_number - xx) ~= 1
                
                yy = yy + (subframe_number - xx);
                
            end
            
            for ii = 1:1:frame_number
    
                for jj = ss:1:xx+ss-1
                    
                    second = [location_table_short(jj*2+2,ii), location_table_short(jj*2+1,ii)]; % second closed vector [E = SF stop, A = SF start] to intersect with first vector
        
                    if jj == 1
                       
                       first = [location_table_short(jj+1,ii), location_table_short(jj,ii)]; % first closed vector [E = SF stop, A = SF start] to intersect with second vector
        
                    elseif ss ~= 1 && ss == jj
                    
                       first = [location_table_short(jj*2,ii), location_table_short(jj*2-1,ii)];
            
                    else
                        
                       first = out;
                       
                    end
                    
                    % Start Author: Xavier Beudaert
                    % calculation of intersection set of two closed vectors [E = SF stop, A = SF start], first and second (see above), repetition via loops above
                    % https://de.mathworks.com/matlabcentral/fileexchange/31753-range-intersection

                    out1(1:(numel(second)+(numel(first)-2)))=0;

                    k=1;

                    while isempty(first)==0 && isempty(second)==0

                        if first(1)>second(1)        
                            temp=second;
                            second=first;
                            first=temp;
                        end

                        if first(2)<second(1)
                            first=first(3:end);
                            continue;

                        elseif first(2)==second(1)
                            out1(k)=second(1);
                            out1(k+1)=second(1);
                            k=k+2;

                            first=first(3:end);
                            continue;

                        else        
                            if first(2)==second(2)        
                                out1(k)=second(1);
                                out1(k+1)=second(2);
                                k=k+2; 

                                first=first(3:end);
                                second=second(3:end);

                            elseif first(2)<second(2)
                                out1(k)=second(1);
                                out1(k+1)=first(2);
                                k=k+2;

                                first=first(3:end);

                            else
                                out1(k)=second(1);
                                out1(k+1)=second(2);
                                k=k+2;

                                second=second(3:end);
                            end
                        end
                    end

                    out=out1(1:k-1);
                    out=out'; % result of vector intersection
                    
                    % End Author: Xavier Beudaert
                
                end
                
                if isempty(out) % if "out" is empty because there is no vector intersection from the investigated vectors than the output is [0,0]
                    out = [0,0];
                else % in all other cases "out" gives the resulting vector intersection out
                    out = out;
                end
                
                    cbu_table_intersect([(2*(xx+ss+yy+subframe_number-1)),(2*(xx+ss+yy+subframe_number-1))-1], ii) = out; % positioning of vector intersection result on right place in cbu_table_intersect in order to assure the arrangement described under CBU_TABLE_INTERSECT
                
            end
        end
    end

  
% CBU_TABLE_UNION shows union sets (set theory!) for different subframe combinations, the existance of other SF within the calculated range is not excluded (see cbu_assi_one and cbu_assi_tow for separated SF combinations)
% cbu_table_union is filled with the correct data by reading out the location_table_short (via loops, see code) in order to generate an specific arrangement of the cbu_table_union structure 
%
% general build-up of the cbu_table_union (row count, column count):
% -----------------------------------------------------------------------------------------------------------------------
% row count = 2*[((subframe_number*subframe_number)+subframe_number)/2] --> gaussian sum formula (example for 5SF)
% row 1/2 = SF1 (row1 = start point A, row2 = stop point E), row 3/4 = SF2, row 5/6 = SF3, row7/8 = SF4, row9/10 = SF5
% row 11/12 = SF1 u SF2, row 13/14 = SF2 u SF3, row 15/16 = SF3 u SF4, row 17/18 = SF4 u SF5 
% row 19/20 = SF1 u SF2 u SF3, row 21/22 = SF2 u SF3 u SF4, row 23/24 = SF3 u SF4 u SF5  
% row 25/26 = SF1 u SF2 u SF3 u SF4, row 27/28 = SF2 u SF3 u SF4 u SF5
% row 29/30 = SF1 u SF2 u SF3 u SF4 u SF5
% column_count = frame_number
% column 1 = F1, column 2 = F2, column 3 = F3, column 4 = F4, ...
% -----------------------------------------------------------------------------------------------------------------------

cbu_table_union = zeros(2*((subframe_number.*subframe_number+subframe_number)/2), frame_number); % blank cbu_table_union filled with zeros at first
cbu_table_union(1:1:subframe_number*2,1:1:frame_number) = location_table_short(1:1:subframe_number*2,1:1:frame_number); % upper part of cbu_table_union consisting of location_table_short

% arrangement of the cbu_table_union structure (rows) is described with an example of a frame with 4SF included (classic signs of set theory do not work in MATLAB so the following abbreviations are used: i = intersection of sets, u = union of sets, \ = difference of sets):
%
% SF combination                C 1
% -----------------------------------> section 1 (number of subdivisions = subframe_number)
% SF1                           [1]
% SF2                           [2]
% SF3                           [3]
% SF4                           [4]
% -----------------------------------> section 2 (number of subdivisions = subframe_number - 1)
% SF1 u SF2                     [5]
% SF2 u SF3                     [6] 
% SF3 u SF4                     [7]
% -----------------------------------> section 3 (number of subdivisions = subframe_number - 2)
% SF1 u SF2 u SF3               [8]
% SF2 u SF3 u SF4               [9]
% -----------------------------------> section 4 (number of subdivisions = subframe_number - 3)
% SF1 u SF2 u SF3 u SF4        [10]
% -----------------------------------
%
% example of cbu_table_union structure (rows) for a frame with 5SF included:
%
% SF combination                C 1
% -----------------------------------> section 1 (number of subdivisions = subframe_number)
% SF1                           [1]
% SF2                           [2]
% SF3                           [3]
% SF4                           [4]
% SF5                           [5]
% -----------------------------------> section 2 (number of subdivisions = subframe_number - 1)
% SF1 u SF2                     [6]
% SF2 u SF3                     [7]
% SF3 u SF4                     [8]
% SF4 u SF5                     [9]
% -----------------------------------> section 3 (number of subdivisions = subframe_number - 2)
% SF1 u SF2 u SF3              [10]
% SF2 u SF3 u SF4              [11]
% SF3 u SF4 u SF5              [12]
% -----------------------------------> section 4 (number of subdivisions = subframe_number - 3)
% SF1 u SF2 u SF3 u SF4        [13]
% SF2 u SF3 u SF4 u SF5        [14]
% -----------------------------------> section 5 (number of subdivisions = subframe_number - 4)
% SF1 u SF2 u SF3 u SF4 u SF5  [15] 
% -----------------------------------

    for ss = 1:1:subframe_number-1 % for loops (ss-loop, xx-loop, ii-loop ,jj-loop) to fill lower part of cbu_table_union by running through all SF scenarios. Resulting order: SF1orSF2 // SF1orSF2orSF3 // SF1orSF2orSF3orSF4 // SF2orSF3 // SF2orSF3orSF4 // SF3orSF4 (exemplary for 4SF)

        for xx = 1:1:subframe_number-ss

            if xx == 1 % for correct row rearrangement of cbu_table_union (see discrepancy between description under CBU_TABLE_UNION and LOOPS, like ss-loop, xx-loop ...) 
                
                yy = 0;
                
            elseif xx == 2
                
                yy = subframe_number-xx;
                
            elseif xx > 2 && (subframe_number - xx) == 1
                
                yy = ((((subframe_number-2).*(subframe_number-2))+(subframe_number-2))/2);
       
            elseif xx > 2 && (subframe_number - xx) ~= 1
                
                yy = yy + (subframe_number - xx);
                
            end
            
            for ii = 1:1:frame_number
    
                for jj = ss:1:xx+ss-1
                    
                    second = [location_table_short(jj*2+2,ii), location_table_short(jj*2+1,ii)]; % second closed vector [E = SF stop, A = SF start] to intersect with first vector
        
                    if jj == 1
                       
                       first = [location_table_short(jj+1,ii), location_table_short(jj,ii)]; % first closed vector [E = SF stop, A = SF start] to intersect with second vector
        
                    elseif ss ~= 1 && ss == jj
                    
                       first = [location_table_short(jj*2,ii), location_table_short(jj*2-1,ii)];
            
                    else
                        
                       first = y;
                       
                    end
                    
                    % Start Author: David Goodmanson
                    % calculation of set union of two closed vectors [E = SF stop, A = SF start], first and second (see above), repetition via loops above
                    % https://de.mathworks.com/matlabcentral/answers/366626-overlapping-time-intervals
                    
                    x = [first;second]; % Nx2 matrix of endpoints x1, x2 of intervals

                    x = sort(x,2); % find union of intervals
                    nrow = size(x,1); 
                    [x ind] = sort(x(:));
                    n = [(1:nrow) (1:nrow)]';
                    n = n(ind);
                    c = [ones(1,nrow) -ones(1,nrow)]';
                    c = c(ind);
                    csc = cumsum(c); % =0 at upper end of new interval(s)
                    irit = find(csc==0);
                    ilef = [1; irit+1];
                    ilef(end) = []; % no new interval starting at the very end 

                    y = [x(ilef) x(irit)]; % y matrix is start and end points of the new intervals, y1,y2
                    ny = [n(ilef) n(irit)]; % ny matrix is the corresponding indices of the start and end points in terms of what row of x they occurred in
                
                    % End Author: David Goodmanson
                    
                end
                
                if size(y,1) > 1 % if row size of "y" is greater than 1 than the investigated vectors do not overlap (the start and stop points of both vectors are transfered to "y") , in this case the union of the vectors is set to [0,0] since the start and stop points of both single vectors are already listed furhter above in "y"
                    y = [0, 0];
                else % in all other cases "y" gives the resulting vector union out
                    y = y;
                end

                cbu_table_union([(2*(xx+ss+yy+subframe_number-1)),(2*(xx+ss+yy+subframe_number-1))-1], ii) = y; % positioning of vector union result on right place in cbu_table_union in order to assure the arrangement described under CBU_TABLE_UNION
            end
        end
    end
    
%% SEGMENT IX - COLOR SECTORS (CALCULATION OF RETINAL POSITION OF CBU RESULTING FROM SUBFRAME POSITION)
% CBU_TABLE (a.k.a. CBU_TABLE_X1, in the following section consequently named CBU_TABLE) shows the retinal positions of all subframe combinations (separated from each other) that result from the previous experimental category and CBU phase provoked by the chosen conditions (like frame rate, eye movement velocitiy ... and so on ...) 
% cbu_table is filled with the correct data by reading out data refering to the intersections and unions of sets for the various subframe combinations from cbu_table_intersect and cbu_table_union and combining them in order to generate an specific arrangement of the cbu_table structure (especially rows)
% the structure of the cbu_table is different for every CBU phase, therefore the codes that fill the cbu_table are different for CBU phase 1 and 2 and 3 (see below) 
% the structure of the cbu_table can be different even within a CBU phase when the number of subframes and/or the relative retinal velocity changes (see below) 
% the structure of the cbu_table is explained under application of the set theory and some of its operations (classic signs of set theory do not work in MATLAB so the following abbreviations are used: i = intersection of sets, u = union of sets, \ = difference of sets)

    % CBU_PHASE 1 calculation of CBU_TABLE 
    % general pattern of the cbu_table is always the same for CBU phase 1 (always one SF more until the maximum SF number is reached, than its always one SF less, see explanation tables below)
    % however, different numbers of subframes lead to a different structure of the cbu_table)

    if cbu_phase == 1 
 
        % general build-up of the cbu_table (row count, column count):
        % -------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        % row count = [(subframe_number-1).*2]+1 
        % row 1/2 = ( SF1 ) \ ( SF2 u SF3 u SF4 u SF5 ) >> row1 = start point (mathematically higher value), row2 = stop point (mathematically lower value), example for 5SF
        % row 3/4 = ( SF1 i SF2 ) \ ( SF3 u SF4 u SF5 ) >> row 3 = start point, row 4 = stop point
        % row 5/6 = ...
        % row 7/8 = ...
        % ...
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        % column_count = frame_number
        % column 1 = F1, column 2 = F2, column 3 = F3, column 4 = F4, ...
        % -------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        
        cbu_table_x1 = zeros(2*(((subframe_number-1).*2)+1), frame_number); % % blank cbu_table filled with zeros at first

        % example of cbu_table structure (rows) for a frame with 4SF included:
        %
        % separated SF combinations                 C 1     C 2    further explanation 
        % -------------------------------------------------------------------------------------------------------
        % ( SF1 ) \ ( SF2 u SF3 u SF4 )             [1]  \  [9]    >> color of SF1
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        % ( SF1 i SF2 ) \ ( SF3 u SF4 )             [5]  \  [7]    >> mix from SF1/SF2
        % ( SF1 i SF2 i SF3 ) \ ( SF4 )             [8]  \  [4)    >> mix from SF1/SF2/SF3
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
        % ( SF1 i SF2 i SF3 i SF4 ) \ ( )           [10] \  [x]    >> white content (all SF mixed together)
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF2 i SF3 i SF4 ) \ ( SF1 )             [9]  \  [1]    >> mix from SF2/SF3/SF4
        % ( SF3 i SF4 ) \ ( SF1 u SF2 )             [7]  \  [5]    >> mix from SF3/SF4
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF4 ) \ ( SF1 u SF2 u SF3 )             [4]  \  [8]    >> color of SF4
        % -------------------------------------------------------------------------------------------------------
        %
        % example of cbu_table structure (rows) for a frame with 5SF included:
        %
        % separated SF combinations                 C 1     C 2    further explanation 
        % -------------------------------------------------------------------------------------------------------
        % ( SF1 ) \ ( SF2 u SF3 u SF4 u SF5 )       [1]  \ [14]    >> color of SF1
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 i SF2 ) \ ( SF3 u SF4 u SF5 )       [6]  \ [12]    >> mix from SF1/SF2
        % ( SF1 i SF2 i SF3 ) \ ( SF4 u SF5 )       [10] \  [9]    >> mix from SF1/SF2/SF3
        % ( SF1 i SF2 i SF3 i SF4 ) \ ( SF5 )       [13] \  [5]    >> mix from SF1/SF2/SF3/SF4
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 i SF2 i SF3 i SF4 i SF5 ) \ ( )     [15] \  [x]    >> white content (all SF mixed together)
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF2 i SF3 i SF4 i SF5 ) \ ( SF1 )       [14] \  [1]    >> mix from SF2/SF3/SF4/SF5
        % ( SF3 i SF4 i SF5 ) \ ( SF1 u SF2 )       [12] \  [6]    >> mix from SF3/SF4/SF5
        % ( SF4 i SF5 ) \ ( SF1 u SF2 u SF3 )       [9]  \ [10]    >> mix from SF4/SF5
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF5 ) \ ( SF1 u SF2 u SF3 u SF4 )       [5]  \ [13]    >> color of SF5
        % -------------------------------------------------------------------------------------------------------
        %
        % aim is to readout the correct SF union and intersection sets from cbu_table_union and cbu_table_intersection in order to combine them and calculate the borders for the separated SF combinations like described in the tables above (see left side), therefore at first the correct positions of the SF combinations within the cbu_table_intersect respectivelly cbu_table_union have to be determined (see number elements on the right side), the determination of the number elements is described as follows:
        % the first column number of the tables middle row element - which is in row position [(subframe_number.*2)-1] in every case - is always [((subframe_number.*subframe_number)+subframe_number)/2] for cbu phase 1, for the condition of 5SF the element would be [15] in row 5, the corresponding second column number is always a blank [x] for the middle row element 
        % for all other table elements there is a simple pattern to fill out the remaining gaps of the table: starting with the first number of the first column in descending row order down to the last element before the already known middle row element, it is the numbers of the first rows (downward movement within table) of the sections of cbu_table_intersect (see tables in previous chapter), since it is clear how many rows are part of every section in cbu_table_intersect it is pretty easy to calculate the numbers assigned to the first row in every section, for the taken example of 5SF the chronological order is [1]-[6]-[10]-[13]       
        % the corresponding elements in second column in the upper part of the table (up to the last element before the already known middle row element) are similarly easy to determine, starting again with the first number of the second column in descending row order down to the last element before the already known middle element, it is the numbers of the last rows (upward movement within table) of the sections of cbu_table_union in this case (see tables in previous chapter), since it is clear how many rows are part of every section in cbu_table_union it is pretty easy to calculate the numbers assigned to the last row in every section, for the taken example of 5SF the chronological order is [14]-[12]-[9]-[5]
        % the whole lower part of the cbu_table (all elements below the middle row) for cbu phase 1 can be fully derived from the already determined upper part of the same cbu_table, the elements of the upper part have to be mirrored on a fictitious horizontal axis going trough the tables middle row and a fictitious vertical axis between the first and the second column, means that the chronological order of the second column elements in the upper part of the cbu_table are recurrent in the first column elements of the lower part of the cbu_table ([14]-[12]-[9]-[5] for 5SF example), furthermore the element order of the first column in the upper part of the cbu_table is repeated in the second column in the lower part of the cbu_table ([1]-[6]-[10]-[13] for the 5SF example) 
        % after all table elements are determined there is one more step to do, combine the determined SF union and intersection sets (with the number elements it is possible to determine a general pattern and write a general code for various conditions) and finally determine the borders of the separated SF combinations for the present conditions (like described in the tables above on the left side), at first the determined table elements (see tables above, right side) are used to read out the correct positions within cbu_table_intersect and cbu_table_union, since every position within cbu_table_intersect / cbu_table_union has two values (start and stop point of the intersection or union of SF) the correct value of both has to be chosen 
        % therefore for the upper part of the cbu_table (above the middle row) the mathematically higher border value of the separated SF combination is determined by choosing the maximum value of both given values (see code below) within the position of cbu_table_intersect refering to the determined number element in the first column (see table above), the mathematically lower border value of the separated SF combination is determined by choosing the maximum value of both given values (see code below) within the position of cbu_table_union refering to the determined number element in the second column (see table above) 
        % the middle row element of the cbu_table is determined by simply choosing both values within the position of the cbu_table_intersect representing the intersection of all available SF for this condition, the position of the correct values is always at the end of the cbu_table_intersect
        % for the lower part of the cbu_table (below the middle row) the mathematically higher border value of the separated SF combination is determined by choosing the minimum value of both given values (see code below) within the position of cbu_table_union refering to the determined number element in the second column (see table above), the mathematically lower border value of the separated SF combination is determined by choosing the minimum value of both given values (see code below) within the position of cbu_table_intersect refering to the determined number element in the first column (see table above)  
        
        for ii = 1:1:frame_number
            
            xx = 1;
            yy = (((subframe_number.*subframe_number)+subframe_number)/2)-1;
            
            for aa = 0:1:(subframe_number-2)
                
                cbu_table_x1((aa+1)*2-1,ii) = max(cbu_table_intersect(xx*2-1:xx*2,ii)); % upper part of the filled out cbu_table (threshold value 1, mathematically higher value)
                cbu_table_x1((aa+1)*2,ii) = max(cbu_table_union(yy*2-1:yy*2,ii)); % upper part of the filled out cbu_table (threshold value 2, mathematically lower value)
                
                xx = xx+subframe_number-aa;
                yy = yy-2-aa;
                
            end
            
            cbu_table_x1(subframe_number.*2-1:subframe_number.*2,ii) = cbu_table_intersect(((subframe_number.*subframe_number)+subframe_number)-1:(subframe_number.*subframe_number)+subframe_number,ii); % white content consisting of all SF
            
            xx = 1;
            yy = (((subframe_number.*subframe_number)+subframe_number)/2)-1;
            
            for aa = 0:1:(subframe_number-2)
                
                cbu_table_x1((aa+subframe_number+1)*2,ii) = min(cbu_table_intersect(yy*2-1:yy*2,ii)); % lower part of the filled out cbu_table (threshold value 1, mathematically lower value)
                cbu_table_x1((aa+subframe_number+1)*2-1,ii) = min(cbu_table_union(xx*2-1:xx*2,ii)); % lower part of the filled out cbu_table (threshold value 2, mathematically higher value)
                
                xx = xx+subframe_number-aa;
                yy = yy-2-aa;
                
            end
            
        end
        

    % CBU_PHASE 2 calculation of CBU_TABLE 
    % other than for CBU phase 1 and 3 the general pattern of the cbu_table is not always the same for CBU phase 2
    % cbu_table pattern changes when the basic points E(SF1) and A(SF2 ... last SF) overlap (see subdivisions below)
    % in addition, different numbers of subframes lead to a different structure of the cbu_table
    % because of the higher complexity some additional tables (cbu_assi_one, cbu_assi_two) are integrated helping to determine the borders of the separated SF combinations in the cbu_table
        
    elseif cbu_phase == 2 
               
        cbu_assi_one = zeros((subframe_number.*subframe_number+subframe_number)/2, 3);
        
        % CBU_ASSI_ONE helps to calculate CBU_TABLE within CBU_PHASE 2 (since it is a bit more complicated than the calculation for CBU phase 1 and 3) 
        % more precisely cbu_assi_one helps to readout the correct SF union and intersection sets from cbu_table_union and cbu_table_intersection in order to combine them and calculate the borders for the separated SF combinations like described below
        % cbu_assi_one contains all possible SF combinations like cbu_table_union and cbu_table_intersect does
        % however, the difference between cbu_assi_one and the previous tables is that within cbu_assi_one the borders for the SF combinations are calculated under the assumption of exclusion of other SF shares (e.g. row 2 in the first example below (4SF) means that only SF2 is available in this area and no other SF like SF1 or SF3 to SF4 are available) 
        % every row within cbu_assi_one consists of three parts (left side), the middle part is always the separated SF combination that has to be calculated, the two parts to the left and to the right are all other SFs that have to be excluded when explicit occurance of the middle part has to be assured 
        % the numbers (right side) behind the description of the three parts to combine (left side) are the positions (rows) of the explicit subframe combinations within cbu_table_intersect and cbu_table_union (see tables above)
        % these numbers make it possible to find periodic patterns for all possible cases (e.g. variable number of subframes) in order to write a general code as compact as possible and consistent for all possible cases 
        %
        % example for 4SF (generally valid for all other condition variations, except SF count)
        %
        % separated SF combinations                            C 1     C 2     C 3
        % ------------------------------------------------------------------------------> section 1 (number of subdivisions = subframe_number)
        % ( ) / ( SF1 ) \ ( SF2 u SF3 u SF4 )               -- [x]  /  [1]  \  [9] -2
        % ( SF1 ) / ( SF2 ) \ ( SF3 u SF4 )                 +4 [1]  /  [2]  \  [7] -3
        % ( SF1 u SF2 ) / ( SF3 ) \ ( SF4 )                 +3 [5]  /  [3]  \  [4] --
        % ( SF1 u SF2 u SF3 ) / ( SF4 ) \ ( )               -- [8]  /  [4]  \  [x] --
        % ------------------------------------------------------------------------------> section 2 (number of subdivisions = subframe_number - 1)
        % ( ) / ( SF1 i SF2) \ ( SF3 u SF4 )                -- [x]  /  [5]  \  [7] -3
        % ( SF1 ) / ( SF2 i SF3) \ ( SF4 )                  +4 [1]  /  [6]  \  [4] --
        % ( SF1 u SF2 ) / ( SF3 i SF4) \ ( )                -- [5]  /  [7]  \  [x] --
        % ------------------------------------------------------------------------------> section 3 (number of subdivisions = subframe_number - 2)
        % ( ) / ( SF1 i SF2 i SF3) \ ( SF4 )                -- [x]  /  [8]  \  [4] --
        % ( SF1 ) / ( SF2 i SF3 i SF4) \ ( )                -- [1]  /  [9]  \  [x] --
        % ------------------------------------------------------------------------------> section 4 (number of subdivisions = subframe_number - 3)
        % ( ) / ( SF1 i SF2 i SF3 i SF4 ) \ ( )             -- [x]  / [10]  \  [x] --
        % ------------------------------------------------------------------------------
        %
        %
        % example for 5SF (generally valid for all other condition variations, except SF count)
        %
        % separated SF combinations                            C 1     C 2     C 3
        % ------------------------------------------------------------------------------> section 1 (number of subdivisions = subframe_number)
        % ( ) / ( SF1 ) \ ( SF2 u SF3 u SF4 u SF5 )         -- [x]  /  [1]  \ [14] -2
        % ( SF1 ) / ( SF2 ) \ ( SF3 u SF4 u SF5 )           +5 [1]  /  [2]  \ [12] -3
        % ( SF1 u SF2 ) / ( SF3 ) \ ( SF4 u SF5 )           +4 [6]  /  [3]  \  [9] -4
        % ( SF1 u SF2 u SF3 ) / ( SF4 ) \ ( SF5 )           +3 [10] /  [4]  \  [5] --
        % ( SF1 u SF2 u SF3 u SF4 ) / ( SF5 ) \ ( )         -- [13] /  [5]  \  [x] --
        % ------------------------------------------------------------------------------> section 2 (number of subdivisions = subframe_number - 1)
        % ( ) / ( SF1 i SF2 ) \ ( SF3 u SF4 u SF5 )         -- [x]  /  [6]  \ [12] -3
        % ( SF1 ) / ( SF2 i SF3 ) \ ( SF4 u SF5 )           +5 [1]  /  [7]  \  [9] -4
        % ( SF1 u SF2 ) / ( SF3 i SF4 ) \ ( SF5 )           +4 [6]  /  [8]  \  [5] --
        % ( SF1 u SF2 u SF3 ) / ( SF4 i SF5 ) \ ( )         -- [10] /  [9]  \  [x] --
        % ------------------------------------------------------------------------------> section 3 (number of subdivisions = subframe_number - 2)
        % ( ) / ( SF1 i SF2 i SF3 ) \ ( SF4 u SF5 )         -- [x]  / [10]  \  [9] -4
        % ( SF1 ) / ( SF2 i SF3 i SF4 ) \ ( SF5 )           +5 [1]  / [11]  \  [5] --
        % ( SF1 u SF2 ) / ( SF3 i SF4 i SF5 ) \ ( )         -- [6]  / [12]  \  [x] --
        % ------------------------------------------------------------------------------> section 4 (number of subdivisions = subframe_number - 3)
        % ( ) / ( SF1 i SF2 i SF3 i SF4 ) \ ( SF5 )         -- [x]  / [13]  \  [5] --        
        % ( SF1 ) / ( SF2 i SF3 i SF4 i SF5 ) \ ( )         -- [1]  / [14]  \  [x] --
        % ------------------------------------------------------------------------------> section 5 (number of subdivisions = subframe_number - 4)
        %  ( ) / ( SF1 i SF2 i SF3 i SF4 i SF5 ) \ ( )      -- [x]  / [15]  \  [x] --
        % ------------------------------------------------------------------------------
        %
        % the periodic patterns of the numbers (right side, three columns) will be explained as follow for the example with 5SF (see above):
        % the second column numbers are consecutive numbers from [1] to [((subframe_number.*subframe_number)+subframe_number)/2] representing the count of all possible combinations for the current SF count
        % the first column numbers ([x]-[1]-[6]-[10]-[13]) of section 1 stand for the SF combinations that have to be excluded (SF combinations with lower numbers), this series of numbers always starts with a blank [x] followed by a [1] for every section, this is always the case no matter how many subframes are included into the calculation, the pattern always starts with the [1] in the second row, to calculate the number of the third row the number of subframes (which is +5 here) has to be added resulting in [6] for the third row, the following rows are calculated by adding numbers that always decreased by 1 until the end of the section is reached
        % the third column numbers ([14]-[12]-[9]-[5]-[x]) of section 1 stand for the SF combinations that have to be excluded (SF combinations with higher numbers), this series of numbers always starts with [(((subframe_number.*subframe_number)+subframe_number)/2)-1], this is always the case without exception, in this case the result is a [14] which the pattern also starts with, to calculate the number of the second row the last number added on the left side (see description above) has to be decreased again by 1 (which is 2, since the last number added on the left side was +3) and than subtracted from [14], consequently the number of the second row is [12], the following rows are calculated by subtracting numbers that continuously decrease by 1 until the end of the section is reached, it has to be considered that the last number of every section on the right side is always a blank [x]
        % this pattern is repeated with little deviation in the following section, reason for the deviation is the reduced row number for every following section, these little deviations can be followed without any further explanation by viewing the examples above
        
        cbu_assi_one(1:end,2) = linspace(1,(subframe_number.*subframe_number+subframe_number)/2,(subframe_number.*subframe_number+subframe_number)/2);
        
        cbu_assi_one(1,3) = (((subframe_number.*subframe_number+subframe_number)/2)-1);
        cbu_assi_one(2,1) = 1;
        
        cc = 0;
        
        for aa = 1:1:subframe_number-2
            
            for bb = 3:1:subframe_number-aa+1
                
                cbu_assi_one(cc+bb,1) = cbu_assi_one(cc+bb-1,1)+subframe_number-(bb-3);
                cbu_assi_one(cc+bb-1,3) = cbu_assi_one(cc+bb-2,3)+(2-bb-aa);
                
            end
            
            cc = cc + (subframe_number-aa+1);
            
            cbu_assi_one(cc+1,3) = cbu_assi_one(cc+1-(subframe_number-aa+1),3)-aa-1;
            cbu_assi_one(cc+2,1) = 1;
            
        end
        
        cbu_assi_two = zeros(((subframe_number.*2)-1),3);
        
        % CBU_ASSI_TWO also helps to calculate CBU_TABLE within CBU_PHASE 2 (since it is a bit more complicated than the calculation for CBU phase 1 and 3) 
        % cbu_assi_two only reads out the relevant SF combinations for the current conditions (cbu phase, subframe count ...) from cbu_assi_one (where all possible SF combinations are included)
        % the chronology of the table rows below respresenting the retinal position of the separated SF combinations that occur for this conditions can also be determined with a fix pattern that is stable for all conditions, this makes it possible to generated the chronology for all conditions with a short but universally valid code (see code below)
        % to explain the chronological pattern of the row, only the second column numbers are relevant since the corresponding first and third column numbers depend on the second column number (connecting pattern between the columns has already been described under the explanation of cbu_assi_one)
        % for cbu_assi_two it is important that the chronological order of the second column numbers is not only different SF counts but rather the some additional sub-categories of the cbu phase 2 
        % these sub-categories are defined by the intersections of the basic points E(SF1) and A(SF2...last SF), the number of sub-categories depends on the subframe count and can be calculated with [subframe_number - 2], means that a presentation of 4SF respectivelly 5SF leads to two resspectivelly three sub-categories (these five sub-categories in total are illustrated exemplarily in the tables below)
        % however, as already stated above, even if there are some cases do distinguish, there is a universal pattern that makes it possible to determine the chronology of the number elements of the second column for all conditions (see description below the tables)
        %
        % example for condition of  4SF ...
        %
        % ...and requirement of A(SF3) > E(SF1) > A(SF4):
        %
        % separated SF combination                           C 1     C 2     C 3        further explanation
        % -------------------------------------------------------------------------------------------------------
        % ( ) / ( SF1 ) \ ( SF2 u SF3 u SF4 )                [x]  /  [1]  \  [9]        >> color of SF1
        % ( ) / ( SF1 i SF2 ) \ ( SF3 u SF4 )                [x]  /  [5]  \  [7]        >> mix from SF1/SF2
        % ( ) / ( SF1 i SF2 i SF3 ) \ ( SF4 )                [x]  /  [8]  \  [4]        >> mix from SF1/SF2/SF3
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        % ( SF1 ) / ( SF2 i SF3) \ ( SF4 )                   [1]  /  [6]  \  [4]        >> mix from SF2/SF3
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 ) / ( SF2 i SF3 i SF4 ) \ ( )                [1]  /  [9]  \  [x]        >> mix from SF2/SF3/SF4
        % ( SF1 u SF2 ) / ( SF3 i SF4 ) \ ( )                [5]  /  [7]  \  [x]        >> mix from SF3/SF4
        % ( SF1 u SF2 u SF3 ) / ( SF4 ) \ ( )                [8]  /  [4]  \  [x]        >> color of SF4
        % -------------------------------------------------------------------------------------------------------
        %
        % ... and requirement of A(SF2) > E(SF1) > A(SF3):
        %
        % separated SF combinations                          C 1     C 2     C 3        further explanation
        % -------------------------------------------------------------------------------------------------------
        % ( ) / ( SF1 ) \ ( SF2 u SF3 u SF4 )                [x]  /  [1]  \  [9]        >> color of SF1
        % ( ) / ( SF1 i SF2 ) \ ( SF3 u SF4 )                [x]  /  [5]  \  [7]        >> mix from SF1/SF2
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        % ( SF1 ) / ( SF2 ) \ ( SF3 u SF4 )                  [1]  /  [2]  \  [7]        >> color of SF2
        % ( SF1 ) / ( SF2 i SF3 ) \ ( SF4 )                  [1]  /  [6]  \  [4]        >> mix from SF2/SF3
        % ( SF1 u SF2 ) / ( SF3 ) \ ( SF4 )                  [5]  /  [3]  \  [4]        >> color of SF3
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 u SF2 ) / ( SF3 i SF4 ) \ ( )                [5]  /  [7]  \  [x]        >> mix from SF3/SF4
        % ( SF1 u SF2 u SF3 ) / ( SF4 ) \ ( )                [8]  /  [4]  \  [x]        >> color of SF4
        % -------------------------------------------------------------------------------------------------------
        %        
        % example for condition of 5SF ...
        %
        % ...and requirement of A(SF4) > E(SF1) > A(SF5):
        %
        % separated SF combinations                          C 1     C 2     C 3        further explanation                                                      
        % -------------------------------------------------------------------------------------------------------
        % ( ) / ( SF1 ) \ ( SF2 u SF3 u SF4 u SF5 )          [x]  /  [1]  \ [14]        >> color of SF1
        % ( ) / ( SF1 i SF2 ) \ ( SF3 u SF4 u SF5 )          [x]  /  [6]  \ [12]        >> mix from SF1/SF2
        % ( ) / ( SF1 i SF2 i SF3 ) \ ( SF4 u SF5 )          [x]  / [10]  \  [9]        >> mix from SF1/SF2/SF3
        % ( ) / ( SF1 i SF2 i SF3 i SF4 ) \ ( SF5 )          [x]  / [13]  \  [5]        >> mix from SF1/SF2/SF3/SF4
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 ) / ( SF2 i SF3 i SF4 ) \ ( SF5 )            [1]  / [11]  \  [5]        >> mix from SF2/SF3/SF4
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 ) / ( SF2 i SF3 i SF4 i SF5 ) \ ( )          [1]  / [14]  \  [x]        >> mix from SF2/SF3/SF4/SF5
        % ( SF1 u SF2 ) / ( SF3 i SF4 i SF5 ) \ ( )          [6]  / [12]  \  [x]        >> mix from SF3/SF4/SF5
        % ( SF1 u SF2 u SF3 ) / ( SF4 i SF5 ) \ ( )          [10] /  [9]  \  [x]        >> mix from SF4/SF5
        % ( SF1 u SF2 u SF3 u SF4 ) / ( SF5 ) \ ( )          [13] /  [5]  \  [x]        >> color of SF5
        % -------------------------------------------------------------------------------------------------------
        %
        % ... and requirement of A(SF3) > E(SF1) > A(SF4)
        %
        % separated SF combinations                          C 1     C 2     C 3        further explanation
        % -------------------------------------------------------------------------------------------------------
        % ( ) / ( SF1 ) \ ( SF2 u SF3 u SF4 u SF5 )          [x]  /  [1]  \ [14]        >> color of SF1
        % ( ) / ( SF1 i SF2 ) \ ( SF3 u SF4 u SF5 )          [x]  /  [6]  \ [12]        >> mix from SF1/SF2
        % ( ) / ( SF1 i SF2 i SF3 ) \ ( SF4 u SF5 )          [x]  / [10]  \  [9]        >> mix from SF1/SF2/SF3
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 ) / ( SF2 i SF3 ) \ ( SF4 u SF5 )            [1]  /  [7]  \  [9]        >> mix from SF2/SF3
        % ( SF1 ) / ( SF2 i SF3 i SF4 ) \ ( SF5 )            [1]  / [11]  \  [5]        >> mix from SF2/SF3/SF4
        % ( SF1 u SF2 ) / ( SF3 i SF4 ) \ ( SF5 )            [6]  /  [8]  \  [5]        >> mix from SF3/SF4
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 u SF2 ) / ( SF3 i SF4 i SF5 ) \ ( )          [6]  / [12]  \  [x]        >> mix from SF3/SF4/SF5
        % ( SF1 u SF2 u SF3 ) / ( SF4 i SF5 ) \ ( )          [10] /  [9]  \  [x]        >> mix from SF4/SF5
        % ( SF1 u SF2 u SF3 u SF4 ) / ( SF5 ) \ ( )          [13] /  [5]  \  [x]        >> color of SF5
        % -------------------------------------------------------------------------------------------------------
        %
        % ... and requirement of A(SF2) > E(SF1) > A(SF3)
        %
        % separated SF combinations                          C 1     C 2     C 3        further explanation
        % -------------------------------------------------------------------------------------------------------
        % ( ) / ( SF1 ) \ ( SF2 u SF3 u SF4 u SF5 )          [x]  /  [1]  \ [14]        >> color of SF1
        % ( ) / ( SF1 i SF2 ) \ ( SF3 u SF4 u SF5 )          [x]  /  [6]  \ [12]        >> mix from SF1/SF2
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 ) / ( SF2 ) \ ( SF3 u SF4 u SF5 )            [1]  /  [2]  \ [12]        >> color of SF2
        % ( SF1 ) / ( SF2 i SF3 ) \ ( SF4 u SF5 )            [1]  /  [7]  \  [9]        >> mix from SF2/SF3
        % ( SF1 u SF2 ) / ( SF3 ) \ ( SF4 u SF5 )            [6]  /  [3]  \  [9]        >> color of SF3
        % ( SF1 u SF2 ) / ( SF3 i SF4 ) \ ( SF5 )            [6]  /  [8]  \  [5]        >> mix from SF3/SF4
        % ( SF1 u SF2 u SF3 ) / ( SF4 ) \ ( SF5)             [10]  / [4]  \  [5]        >> color of SF4
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( SF1 u SF2 u SF3 ) / ( SF4 i SF5 ) \ ( )          [10] /  [9]  \  [x]        >> mix from SF4/SF5
        % ( SF1 u SF2 u SF3 u SF4 ) / ( SF5 ) \ ( )          [13] /  [5]  \  [x]        >> color of SF5
        % -------------------------------------------------------------------------------------------------------
        %
        % the general pattern to calculate the second column numbers (as the basis for the transfer of first and third column number into the table) is explained as follows:
        % at first sight there is no apparent general pattern that fits the different conditions (illustrated exemplarily for 4SF and 5SF in the five tables above), however, with some further thinking it becomes clear that the first numbers of the second column are the first rows of the sections of cbu_table_intersect, since it is clear how many rows are part of every section of cbu_table_intersect it is pretty easy to calculate the numbers assigned to the first row in every section, it depends on the sub-category how many number can be determined like this, for the first sub-category it is [subframe - 1 ] numbers that can be determined like this, e.g. the first four numbers of the second column within cbu_table can be determined for the 5SF condition and the requirement of [A(SF4) > E(SF1) > A(SF5)] for the first sub-category (see table above), for this example the numbers are [1]-[6]-[10]-[13], every following sub-category allows one number less to be calculated like explained, means that the second sub-category of the 5SF condition allows to determine the first three numbers like this ... and so on       
        % luckily almost exactly the same procedure can be done when determining the last numbers of the cbu_tables second column, the only difference is that the numbers are transfered from the cbu_table_intersection by taking not the numbers of the first but the last rows of every section of cbu_table_intersect, again these numbers can be calculated without any further effort since it is well known how many rows each section of cbu_table_intersect has (and also established that the numbers of this column are in chronological order in upwards direction starting at [1]), for the previous example of the 5SF condition and the first sub-category the numbers would be [5]-[9]-[12]-[14] starting from the last position upwards within the table, since the cbu_table for this condition has nine rows in total there is only one row left in the middle of the table which has no number assigned (for other conditions the gap of not assigned numbers it larger since the amount of numbers that can be determined with the explaind procedure decreases with further sub-categories (like already explained), to get a better overview the numbers that can be determined with the explained procedure are separated from the remaining previously unknown numbers in the middle visually by a dashed line in the cbu_table (see above)   
        % the remaining numbers in the middle of the table can be calculated with another procedure based on another pattern, the number of the first middle row element directly below the dashed line can be determined by going two rows up, taking the already established number and increasing it by the value of 1, for the previous example of 5SF condition / first sub-category the previously unkown number in row 5 would be [11] since the number assigned two rows above in row 3 is [10], now this cbu_table is completed, all numbers in second column are determined, as already mentioned some cbu_tables refering to other conditions than in the example taken up do have larger gaps, in this case the same procedure is repeated until every row within the gap is filled with a number, means that the second position in the gap is assigned with a number by again going up two rows, taking the assigned number in this row and increasing the number again by the value of 1, in case of more than two gap positions the procedure is repeated until all posiions in the gap are filled with a number, the described two procedures are sufficient to fill all second columns of all cbu_tables for all available conditions (subframe number, cbu phase and their sub-categories ...) 
        % after all numbers in second column are assigned the first and third column of the cbu_table can be determined without further effort since the every second column number has a fix neighbour that can be determined via cbu_table_intersect / cbu_table_union, the corresponding first and third column numbers to the number in second column are simply looked up in cbu_table_intersect / cbu_table_union and transfered to cbu_table.
        
        cbu_assi_two(1,1:end) = cbu_assi_one(1,1:end);
        cbu_assi_two(end,1:end) = cbu_assi_one(subframe_number,1:end);
        
        for dd = subframe_number:-1:3
            
            if subframe_location_table_hor_x1((dd-1)*4+1,1) < subframe_location_table_hor_x1(3,1) && subframe_location_table_hor_x1(3,1) < subframe_location_table_hor_x1((dd-2)*4+1,1)
                
                for ee = 1:1:(dd-2)
                    
                    cbu_assi_two(ee+1,1) = cbu_assi_one(cbu_assi_two(ee,2)+subframe_number-(ee-1),1);
                    cbu_assi_two(ee+1,2) = cbu_assi_two(ee,2)+subframe_number-(ee-1);
                    cbu_assi_two(ee+1,3) = cbu_assi_one(cbu_assi_two(ee,2)+subframe_number-(ee-1),3);
                    
                    cbu_assi_two(end-ee,1) = cbu_assi_one(cbu_assi_two(end-(ee-1),2)+subframe_number-ee,1);
                    cbu_assi_two(end-ee,2) = cbu_assi_two(end-(ee-1),2)+subframe_number-ee;
                    cbu_assi_two(end-ee,3) = cbu_assi_one(cbu_assi_two(end-(ee-1),2)+subframe_number-ee,3);
                    
                end
                
                for gg = 1:1:((subframe_number.*2)-1)-((ee+1)*2)
                    
                    cbu_assi_two(ee+2+(gg-1),1) = cbu_assi_one(cbu_assi_two(ee+gg-1,2)+1,1);
                    cbu_assi_two(ee+2+(gg-1),2) = cbu_assi_two(ee+gg-1,2)+1;
                    cbu_assi_two(ee+2+(gg-1),3) = cbu_assi_one(cbu_assi_two(ee+gg-1,2)+1,3);
                    
                end
                
            end
            
        end
                 
        cbu_table_x1 = zeros(2*((subframe_number.*2)-1), frame_number); % blank cbu_table filled with zeros at first
        
        % CBU_TABLE stores the borders (mathematically higher and lower value) of the separated SF combinations for the previous condition (whereas cbu_assi_one and cbu_assi_two do only contain the number elements that are necessary to read out the correct SF combinations on the way to calculate the borders of the separated SF combinations)
        % by reading out the correct SF combinations (with the help of the determined element numbers within cbu_assi_two) and combining them in order to calculate the borders for the separted SF combinations the observer perception representing the different SF combination stimulations on differnt areas of the observers retina can be simulated
        %
        % general build-up of the cbu_table (row count, column count):
        % ---------------------------------------------------------------------------------------------------------------------------------------------------
        % row count = (subframe_number.*2)-1 
        % row 1/2 = ( SF1 ) \ ( SF2 u SF3 u SF4 u SF5 ) >> row1 = start point (mathematically higher value), row2 = stop point (mathematically lower value)
        % row 3/4 = ( SF1 i SF2 ) \ ( SF3 u SF4 u SF5 ) >> row 3 = start point, row 4 = stop point
        % row 5/6 = ...
        % row 7/8 = ...
        % ...
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % column_count = frame_number
        % column 1 = F1, column 2 = F2, column 3 = F3, column 4 = F4, ...
        % ---------------------------------------------------------------------------------------------------------------------------------------------------        
        %       
        % example tables for determination of separated SF combinations for 4SF and 5SF with its sub-categories have already been illustrated in the previous section during description of cbu_assi_two, see tables above for illustration of the tables and their structure
        
        for xx = 1:1:frame_number
            
            for aa = 1:1:(subframe_number.*2)-1
                
                if cbu_assi_two(aa,1) ~= 0 && cbu_assi_two(aa,3) ~= 0
                    
                    cbu_table_x1(aa+aa-1,xx) = min(cbu_table_union(2*cbu_assi_two(aa,1),xx),cbu_table_intersect(2*cbu_assi_two(aa,2)-1,xx));
                    cbu_table_x1(aa+aa,xx) = max(cbu_table_intersect(2*cbu_assi_two(aa,2),xx),cbu_table_union(2*cbu_assi_two(aa,3)-1,xx));
                    
                elseif cbu_assi_two(aa,1) == 0
                    
                    cbu_table_x1(aa+aa-1,xx) = cbu_table_intersect(2*cbu_assi_two(aa,2)-1,xx);
                    cbu_table_x1(aa+aa,xx) = max(cbu_table_intersect(2*cbu_assi_two(aa,2),xx),cbu_table_union(2*cbu_assi_two(aa,3)-1,xx));
                    
                elseif cbu_assi_two(aa,3) == 0
                    
                    cbu_table_x1(aa+aa-1,xx) = min(cbu_table_union(2*cbu_assi_two(aa,1),xx),cbu_table_intersect(2*cbu_assi_two(aa,2)-1,xx));
                    cbu_table_x1(aa+aa,xx) = cbu_table_intersect(2*cbu_assi_two(aa,2),xx);
                    
                end
                
            end
            
        end
       
    % CBU_PHASE 3 calculation of CBU_TABLE 
    % general pattern of the cbu_table is always the same for CBU phase 3 (SF - GAP - SF - GAP - SF ...) 
    % however, different numbers of subframes change the length (rows) of the cbu_table     
        
    elseif cbu_phase == 3 
  
        % general build-up of the cbu_table (row count, column count):
        % -----------------------------------------------------------------------------------------------------------------------
        % row count = (subframe_number.*2)-1 
        % row 1/2 = ( SF1 ) >> row1 = start point (mathematically higher value), row2 = stop point (mathematically lower value)
        % row 3/4 = ( GP1/2 ) >> row 3 = start point, row 4 = stop point
        % row 5/6 = ...
        % row 7/8 = ...
        % ...
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        % column_count = frame_number
        % column 1 = F1, column 2 = F2, column 3 = F3, column 4 = F4, ...
        % -----------------------------------------------------------------------------------------------------------------------        
        
        cbu_table_x1 = zeros(2*((subframe_number.*2)-1), frame_number); % blank cbu_table filled with zeros at first

        % example of cbu_table structure (rows) for a frame with 4SF included (creation of structure for SF counts other than SF4 is trivial and therefore not documented): 
        %
        % separated SF combinations        C 1          further explanation
        % -----------------------------------------------------------------------------------------------------------------------
        % ( ) / ( SF1 ) \ ( )              [1]          >> color of SF1
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( ) / ( GP1/2 ) \ ( )             -           >> gap between SF1/SF2
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( ) / ( SF2 ) \ ( )              [2]          >> color of SF2
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( ) / ( GP2/3 ) \ ( )             -           >> gap between SF2/SF3
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( ) / ( SF3 ) \ ( )              [3]          >> color of SF3
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( ) / ( GP3/4 ) \ ( )             -           >> gap between SF3/SF4
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % ( ) / ( SF4 ) \ ( )              [4]          >> color of SF4
        % -----------------------------------------------------------------------------------------------------------------------
        %
        % the determination of separated SF combinations for cbu phase 3 is trivial because there are no overlappings between the single SF, means that there are actually no combinations of SF to calculate but rather SF coordinates to assign from location_table 
        % the only step is to choose the correct basic start and stop points of all available SF from the location_table and also assigning the correct start and stop points to the gaps between the SF that certainly have a direct connection to the start respectivelly stop points of the SF, means that these values can also be assigned from the location_table (see code below)
        
        for aa = 1:1:subframe_number
            
            cbu_table_x1(aa*4-3:1:aa*4-2,1:1:frame_number) = location_table_short(aa*2-1:1:aa*2,1:1:frame_number); % ...
            
        end
        
        for bb = 1:1:subframe_number-1
            
            cbu_table_x1(bb*4-1,1:1:frame_number) = cbu_table_x1(bb*4-2,1:1:frame_number); % ...
            cbu_table_x1(bb*4,1:1:frame_number) = cbu_table_x1(bb*4+1,1:1:frame_number); % ...
            
        end
        
    end
   
    % CBU_TABLE_MM >> transition of retinal position specification from angle value in [DEG] (see cbu_table_x1 a.k.a cbu_table in the sections above) to distance value in [MM] (see cbu_table_mm below)
    % the calculations for retinal positions of the displayed content above are made under the assumption that the optical path goes through the first (N) and second nodal point (N´) of the viewers eye which has the benefit that the angle of the incident light ray is similar to the angle of the outgoing light ray - no further calculation is necessary (see cbu_table_x1), this works perfectly well as long as relevant retinal positions can to be expressed in [DEG]
    % however, if the retinal position has to be stated in [MM] a conversion from the angle in [DEG] to the distance in [MM] for the retinal position is not trivial, the problem is that the retinal position in [DEG] (see cbu_table_x1) can not be transfered into retinal position in [MM] via simple calculation of the circular arc by using the the parameter (1) retinal radius (see RR_m) and (2) angle x1 (see cbu_table_x1) starting from the second nodal point (N´), the reason is that the midpoint of the circle (Z) representing the retinal shape (approximation!) is not equal to the second nodal point (N´), which is relevant for the transition of the optical path from the display unit through the optics of the eye to the retina, instead angle x1 has to be converted into angle x2 (see code for cbu_table_x2 below) at first, afterwards the calculation of the circular arc can be operated by using the parameter (1) retinal radius (see RR_m) and (2) angle x2 (see cbu_table_x2) starting from the midpoint of the circle representing the retinal shape, the result is the retinal position as distance to the point of origin of the coordinate system (fovea centralis) in [MM] (see cbu_table_mm) >> there is also a graphical illustration of the geometrical relations >> see confluence >> theoretical cbu model >> display-to-retina correction (helpful for a better understanding)
    
    for ii = 1:1:size(cbu_table_x1,2)
        
        for jj = 1:1:size(cbu_table_x1,1)
            
            if cbu_table_x1(jj,ii) < 0
                
                cbu_table_x2(jj,ii) = (90-(180-(180-(180-abs(cbu_table_x1(jj,ii))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(cbu_table_x1(jj,ii))-90))))).*(tan((deg2rad(abs(cbu_table_x1(jj,ii))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))))).*(-1);
                cbu_table_mm(jj,ii) = ((cbu_table_x2(jj,ii)./360).*2.*pi.*RR_m).*1000;
                
            else
                
                cbu_table_x2(jj,ii) = 90-(180-(180-(180-abs(cbu_table_x1(jj,ii))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(cbu_table_x1(jj,ii))-90))))).*(tan((deg2rad(abs(cbu_table_x1(jj,ii))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))));
                cbu_table_mm(jj,ii) = ((cbu_table_x2(jj,ii)./360).*2.*pi.*RR_m).*1000;
            
            end
            
        end
        
    end
    
%% SEGMENT X - TEMPORAL SUMMATION (OF VISUAL STIMULATION BY SINGLE SUBFRAMES)
% THRESHOLD VALUE FOR TEMPORAL SUMMATION of visual stimuli described by Bloch's Law is not determined at a fix value - e.g. 0.05sec for cones and 0.10sec for rods (Hood and Finkelstein, as cited in Blake and Sekuler, 2006, p.98) - instead, temporal summation is executed for all subframes of one full FRAME CYCLE only in order to calculate the color characteristic and intensity profile of the resulting retinal stimulation >> see also GENERAL ASSUMPTIONS
% FRAME CYCLE that is chosen for execution of temporal summation is the FRAME CYCLE with the largest CBU effect within the defined sequence (which is the frame cycle with the highest eye movement velocity in [DEG/SEC])
% number of INTEGRATED_SUBFRAMES is equal to the number of subframes included within the chosen FRAME CYCLE
% SCAN_START defines the retinal start point of temporal summation and is defined as the mathematically highest value from all basic point's of the INTEGRATED_SUBFRAMES rhomboids
% SCAN_STOP defines the retinal stop point of temporal summation and is defined as the mathematically lowest value from all basic point's of the INTEGRATED_SUBFRAMES rhomboids
% SCAN_RESOLUTION determines the number of retinal spots the temporal summation is executed within the start point of scanning (SCAN_START) and stop point of scanning (STOP_POINT)
% final tables - XYL_TABLE and XYZ_TABLE - show the result of temporal summation (color, intensity = Y_VALUES) for all investigated retinal spots (X_VALUES)
%
% COLUMNS of the final tables are structured as follows ... 
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% column 1 = X_VALUE = retinal position of temporal summation, range of X_VALUES goes from SCAN_START to SCAN_STOP, number of X_VALUES is defined by SCAN_RESOLUTION
% column 2/3/4 = first integrated subframe's chromaticity coordinates (x = column 2 ,y = column 3) and intensity value (L = column 4) for XYL_TABLE respectivelly tristimulus values (X = column 2, Y = column 3, Z = column 4) for XYZ_TABLE at a certain retinal position (X_VALUE) 
% column 5/6/7 = second integrated subframe's chromaticity coordinates (x = column 2 ,y = column 3) and intensity value (L = column 4) for XYL_TABLE respectivelly tristimulus values (X = column 2, Y = column 3, Z = column 4) for XYZ_TABLE at a certain retinal position (X_VALUE)
% ...
% column(end-2)/(end-1)/(end) = chromaticity coordinates (x = column 2 ,y = column 3) and intensity value (L = column 4) for XYL_TABLE respectivelly tristimulus values (X = column 2, Y = column 3, Z = column 4) for XYZ_TABLE after temporal summation of all relevant integrated subframes at a certain retinal position (X_VALUE)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% column count = ( integrated_subframes .* 3 ) + 4
%
% ROWS of the final tables are structured as follows ...
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% every row stands for one retinal position (X_VALUE) at which temporal summation is executed 
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% row count (amount of investigated X_VALUES) = ( scan_start - scan_stop ) / scan_resolution (approximately, basic points of integrated subframes rhomboids as important turning points are included additionally) 


for bb = 1:1:frame_number % SUBFRAME_LOCATION_TABLE_HOR (shows retinal position of basis points of rhomboid that describes subframe --> A´s, C´s, E´s, G´s) is restructured to SUBFRAME_LOCATION_TABLE_HOR_SORTED1; values of the basic points of every subframe are arranged in descending order regarding mathematical value --> A´s, G´s, C´s, E´s 
    
    for aa = 1:1:subframe_number
        
        subframe_location_table_hor_sorted1((aa.*4)-3:(aa.*4),bb) = sortrows(subframe_location_table_hor_x1((aa.*4)-3:(aa.*4),bb),-1);
        
    end
    
end

subframe_time_table_sorted = zeros(subframe_number.*4, frame_number); % SUBFRAME_TIME_TABLE (shows time of retinal presentation of basis points of rhomboid that describes subframe --> A´s, C´s, E´s, G´s) is restructured to SUBFRAME_TIME_TABLE_SORTED; values of the basic points of every subframe are arranged in corresponding order to the new order of SUBFRAME_LOCATION_TABLE_HOR_SORTED1 --> A´s, G´s, C´s, E´s  

for bb = 1:1:frame_number
    
    for aa = 1:1:subframe_number
        
        subframe_time_table_sorted((aa.*4)-2:(aa.*4)-1,bb) = subframe_on_sec;
        subframe_integrated_luminance_table((aa.*4)-3:(aa.*4),bb) = (subframe_time_table_sorted((aa.*4)-3:(aa.*4),bb)).*sf_luminance(aa,1); % SUBFRAME_TIME_TABLE_SORTED is converted to SUBFRAME_INTEGRATED_LUMINANCE_TABLE by further multiplication of subframe presentation time and subframe intensity (luminance); results in the perceived visual intensity for one subframe after summation of the intensity over time (temporal summation); assumption that the threshold value for temporal summation after Bloch's law does not fall below the presentation time for one single subframe; this is highly unlikely since the subframe rate would have to fall below 6.0HZ for an assumed duty cycle of 0.3 and a fix Bloch Time of 0.05sec (such subframe rates are blocked, no calculation with model possible)
        
    end
    
end

subframe_location_table_hor_sorted2 = subframe_location_table_hor_sorted1(:); % SUBFRAME_LOCATION_TABLE_HOR_SORTED1 is restructured again (to SUBFRAME_LOCATION_TABLE_HOR_SORTED2) by putting all columns in one column (better for further calculation)
subframe_integrated_luminance_table_sorted = subframe_integrated_luminance_table(:); % SUBFRAME_INTEGRATED_LUMINANCE_TABLE is restructured (to SUBFRAME_INTEGRATED_LUMINANCE_TABLE_SORTED) by putting all columns in one column (better for further calculation)

integrated_subframes = subframe_number; % calculation of the number of subframes whose single intensities have to be temporally summarized

% alternative calculation INTEGRATED_SUBFRAMES (deactivated) 
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% integrated_subframes = floor((bloch_time_cones_sec+subframe_off_sec)./subframe_duration_sec); % calculation of the number of subframes whichs single intensities have to be temporally summarized under consideration of the defined threshold limit of temporal summation of cones (see GENERAL ASSUMPTIONS section at the start); calculated value of INTEGRATED_SUBFRAMES is rounded down if it is not an integer number in order to only include subframes that lie fully within the borders of temporal summation
%  
% if integrated_subframes > (frame_number.*subframe_number) % if calculated number for integrated subframes is higher than the actual number of subframes relevant for the defined conditions than the number of integrated subframes is equalized with the actual number of subframes relevant for the defined conditions (maximum number of subframes)
%      
%     integrated_subframes = (frame_number.*subframe_number);
%  
% end
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

scan_resolution = 0.0005; % defines the distance between the scanning lines on the retina >> the lower the scan resolution the higher the degree of precision for calculation of intensity and color profile of the CBU stimulus

scan_start = max(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1)+1):(integrated_subframes.*4.*eye_velo_hor_variable_degs_table_index),1)); % mathematically highest value from all basic point's retinal positions of the INTEGRATED_SUBFRAMES rhomboids
scan_stop = min(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1)+1):(integrated_subframes.*4.*eye_velo_hor_variable_degs_table_index),1)); % mathematically lowest value from all basic point's retinal positions of the INTEGRATED_SUBFRAMES rhomboids

x_value_table = [scan_start:-scan_resolution:scan_stop]'; % X_VALUE_TABLE defines the X_VALUES (retinal positions) for which the determination of intensity and color profile is executed (1-column matrix); starting from mathematically highest value of all basic points of the rhomboids of the INTEGRATED SUBFRAMES; stepwise reduction of mathematical value of the following X_VALUES by subtraction of SCAN_RESOLUTION; stopping at mathematically lowest value of all basic points of the rhomboids of the INTEGRATED SUBFRAMES
x_value_table =  sortrows([x_value_table; (subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1)+2):(integrated_subframes.*4.*eye_velo_hor_variable_degs_table_index),1))],-1); % further integration of all basic points of the rhomboids of the INTEGRATED SUBFRAMES (so far only basic point with mathematically highest value definitely included) in order to guarantee that these points as important turning points (especially for calculation of intensity) are part of the X_VALUES

uu = 0;

for xx = 1:1:size(x_value_table,1) % downwards selection of X_VALUES out of X_VALUE_TABLE one by one from row to row   
    
    x_value = x_value_table(xx,1);
    
    uu = uu + 1;
    
    for ss = 1:1:integrated_subframes % determination of the perceived luminance level (named Y_VALUE) for all X_VALUES within X_VALUE_TABLE
        
        if x_value > subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-3,1) % in order to determine the perceived luminance levels (Y_VALUE) at the actual retinal position (X_VALUE) different areas have to be distinguished (areas with zero, increasing, stable and decreasing perceived luminance level behaviour); non-zero perceived luminance level areas are described with a simple linear equation ( y = mx + t ) >> see code below 
            
            y_value_table(uu,ss) = 0.0;
            
        elseif x_value <= subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-3,1) && x_value > subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1)
            
            y_value_table(uu,ss) = ((subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-3,1)-subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1))./(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-3,1)-subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1))).*x_value+(subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1)-((subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-3,1)-subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1))./(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-3,1)-subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1))).*subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1));
            
        elseif  x_value <= subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1) && x_value > subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1)
            
            y_value_table(uu,ss) = ((subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1)-subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1))./(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1)-subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1))).*x_value+(subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1)-((subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1)-subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1))./(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-2,1)-subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1))).*subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1));
            
        elseif  x_value <= subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1) && x_value > subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1)
            
            y_value_table(uu,ss) = ((subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1)-subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1))./(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1)-subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1))).*x_value+(subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1)-((subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1)-subframe_integrated_luminance_table_sorted((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1))./(subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4)-1,1)-subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1))).*subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1));
            
        elseif x_value <= subframe_location_table_hor_sorted2((integrated_subframes.*4.*(eye_velo_hor_variable_degs_table_index-1))+(ss.*4),1)
            
            y_value_table(uu,ss) = 0.0;
            
        end
        
    end
    
end

xyL_table = zeros(uu,3.*integrated_subframes+4); % XYL_TABLE includes X_VALUES (first column >> retinal position) and the corresponding chromaticity coordinates x and y (taken from the DATA_ENTRY_MASK) as well as the perceived luminance L (Y_VALUES, calculated with the code above) for all INLCUDED SUBFRAMES (detailled structure of matrix ... see notes at the beginning of the section)  

xyL_table(1:end,1) = x_value_table;

for aa = 1:1:integrated_subframes
    
    if aa <= subframe_number
        
        bb = aa;
    
    elseif aa > subframe_number
    
        bb = aa-(floor(aa./subframe_number).*subframe_number);
        
        if bb == 0
           
            bb = 3;
        
        end
    
    end
        
    xyL_table(1:end,aa.*3-1) = sf_xcoordinate(bb,1);
    
    xyL_table(1:end,aa.*3) = sf_ycoordinate(bb,1);
    
    xyL_table(1:end,aa.*3+1) = y_value_table(1:end,aa);
    
end

XYZ_table = zeros(uu,3.*integrated_subframes+4); % XYZ_TABLE includes X_VALUES (first column >> retinal position) and the corresponding tristimulus values X, Y and Z (calculated from the existing xyL values, see code below) for all INLCUDED SUBFRAMES (detailled structure of matrix ... see notes at the beginning of the section) 

XYZ_table(1:end,1) = x_value_table; 

% determination of the color and luminance specification of the resulting perception from blending of the relevant single SFs at a certain retinal spot (which color and intensity is perceived by the observer after temporal summation of the defined single subframes of the chosen frame cycle at a certain retinal spot), for determination of the final perception (from temporally summing up the relevant subframes) a transition of chromaticity coordinates (xyY) to tristimulus values (XYZ) back to chromaticity coordinates (xyY) is necessary (see calculation below), color and luminance specification for the final perception are positioned in the last three columns of the XYZ_TABLE (also valid for XYL_TABLE)

for aa = 1:1:integrated_subframes % transfer xyY (or xyL) values for single SFs (defined above) to XYZ values for single SFs
    
    XYZ_table(1:end,aa.*3-1) = xyL_table(1:end,aa.*3-1).*(xyL_table(1:end,aa.*3+1)./xyL_table(1:end,aa.*3));
    
    XYZ_table(1:end,aa.*3) = xyL_table(1:end,aa.*3+1);
    
    XYZ_table(1:end,aa.*3+1) = (xyL_table(1:end,aa.*3+1)./xyL_table(1:end,aa.*3)).*(1-xyL_table(1:end,aa.*3-1)-xyL_table(1:end,aa.*3));

end

for bb = 1:1:integrated_subframes % transfer XYZ values for single SFs to XYZ values for mixed SFs (via simple addition of single XYZ values) 
   
    XYZ_table(1:end,end-2) = XYZ_table(1:end,end-2) + XYZ_table(1:end,bb.*3-1);
    
    XYZ_table(1:end,end-1) = XYZ_table(1:end,end-1) + XYZ_table(1:end,bb.*3);
    
    XYZ_table(1:end,end) = XYZ_table(1:end,end) + XYZ_table(1:end,bb.*3+1);
    
end

for cc = 1:1:integrated_subframes % transfer XYZ values of mixed SFs to xyY values of mixed SFs
   
    xyL_table(1:end,end-2) = XYZ_table(1:end,end-2)./(XYZ_table(1:end,end-2)+XYZ_table(1:end,end-1)+XYZ_table(1:end,end));   
    
    xyL_table(1:end,end-1) =  XYZ_table(1:end,end-1)./(XYZ_table(1:end,end-2)+XYZ_table(1:end,end-1)+XYZ_table(1:end,end));
    
    xyL_table(1:end,end) = XYZ_table(1:end,end-1);
    
end

%% SEGMENT XI - MODEL INDICES 
% compress the previous calculations on cbu perception (for the chosen frame cycle with the highest cbu effect) down to some single indices (cbu score, non-cbu score, reference score)
% necessary to make the output of the theoretical model comparable to the output of the empirical investigations (e.g. COBUS participants evaluated cbu score on a continuous scale ranging from 1 to 5) 
% two options for index calculation are presented (a third option is not worked out so far), see following code including further description
% cbu score, non-cbu score and reference score of OPTION 1 and OPTION 2 are calculated on basis of [MM] and also [DEG] expression for retinal stimulation

% OPTION 1 - area of retinal stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description

    % CBU_SCORE >> gives out the area of retinal stimulation by calculating the surface of cbu provokation in [MM2] respectively [DEG2]
    
    cbu_height_opt1_mm = subframe_location_table_ver_mm(1,1)-subframe_location_table_ver_mm(2,1); % vertical dimension of the retinal cbu stimulation in [MM]
    cbu_height_opt1_deg = subframe_location_table_ver_x1(1,1)-subframe_location_table_ver_x1(2,1); % vertical dimension of the retinal cbu stimulation [DEG]

    if cbu_phase == 1 % horizontal dimension of the retinal cbu stimulation for cbu phase 1 in [MM] respectively [DEG] (consideration of borders between color biased cbu stimulation and unbiased stimulation) 
        
        cbu_width_opt1_mm = (cbu_table_mm(1,eye_velo_hor_variable_degs_table_index)-cbu_table_mm(subframe_number.*2-1,eye_velo_hor_variable_degs_table_index))+(cbu_table_mm(subframe_number.*2,eye_velo_hor_variable_degs_table_index)-cbu_table_mm(end,eye_velo_hor_variable_degs_table_index));        
        cbu_width_opt1_deg = (cbu_table_x1(1,eye_velo_hor_variable_degs_table_index)-cbu_table_x1(subframe_number.*2-1,eye_velo_hor_variable_degs_table_index))+(cbu_table_x1(subframe_number.*2,eye_velo_hor_variable_degs_table_index)-cbu_table_x1(end,eye_velo_hor_variable_degs_table_index));
        
    elseif cbu_phase == 2 % horizontal dimension of the retinal cbu stimulation for cbu phase 2 in [MM] respectively [DEG] (no consideration of borders between color biased cbu stimulation and unbiased stimulation since there is no more unbiased content during cbu phase 2)
        
        cbu_width_opt1_mm = cbu_table_mm(1,eye_velo_hor_variable_degs_table_index)-cbu_table_mm(end,eye_velo_hor_variable_degs_table_index); 
        cbu_width_opt1_deg = cbu_table_x1(1,eye_velo_hor_variable_degs_table_index)-cbu_table_x1(end,eye_velo_hor_variable_degs_table_index); 
        
    elseif cbu_phase == 3 % horizontal dimension of the retinal cbu stimulation for cbu phase 3 in [MM] respectively [DEG] (color biased cbu stimulation surface and also whole gap surface w/o weighting factor is considered >> see confluence for further information)
        
        gap_weight_opt1_phase3 = 1.0; % calculation of the the widths of cbu during cbu phase 3 must include an approach that defines how to handle the gaps between the retinal subframe stimulation, at present the gaps are fully integrated into the determination of the widths of cbu which means that the gaps have the same weight as the subframe stimulation (weighting factor is set to 1.0 by default, see confluence for further description)
        
        aa = 0; % calculation on basis of [MM] expression
        
        for jj = 1:1:subframe_number
            
            aa = aa + cbu_table_mm(jj.*4-3,eye_velo_hor_variable_degs_table_index)-cbu_table_mm(jj.*4-2,eye_velo_hor_variable_degs_table_index); % aa represents the horizontal share of retinal cbu stimulation by the subframes
        
        end
        
        bb = (cbu_table_mm(1,eye_velo_hor_variable_degs_table_index)-cbu_table_mm(end,eye_velo_hor_variable_degs_table_index))-aa; % bb represents the horizontal share of gaps between the retinal subframe stimulation
        
        cbu_width_opt1_mm = aa.*(1+(bb/aa).*gap_weight_opt1_phase3); % calculation of the cbu width of (spatially separated) subframes (see aa, gaps do no contribute to this widths calculation), however, it is assumed that gaps have an effect on the perceived cbu effect, additional multiplication with factor (1+(bb/aa).*gap_weight_opt1_phase3) that allows to include gaps for determination of the cbu widths and following the final cbu score, setting gap_weight_opt1_phase3 allows to define different weightings of the gaps in comparison to the subframe stimulation (e.g. the default value of 1.0 means that gaps are considered for the calculation of the cbu widths and concluding for the final cbu score equally to the subframe stimulation whereas a value of 0.5 would mean that gaps are only expected to contribute with 50% in comparison to the subframes stimulation to the final cbu score) 

        aa = 0; % calculation on basis of [DEG] expression
        
        for jj = 1:1:subframe_number
            
            aa = aa + cbu_table_x1(jj.*4-3,eye_velo_hor_variable_degs_table_index)-cbu_table_x1(jj.*4-2,eye_velo_hor_variable_degs_table_index); % aa represents the horizontal share of retinal cbu stimulation by the subframes
        
        end
        
        bb = (cbu_table_x1(1,eye_velo_hor_variable_degs_table_index)-cbu_table_x1(end,eye_velo_hor_variable_degs_table_index))-aa; % bb represents the horizontal share of gaps between the retinal subframe stimulation
        
        cbu_width_opt1_deg = aa.*(1+(bb/aa).*gap_weight_opt1_phase3); % calculation of the cbu width of (spatially separated) subframes (see aa, gaps do no contribute to this widths calculation), however, it is assumed that gaps have an effect on the perceived cbu effect, additional multiplication with factor (1+(bb/aa).*gap_weight_opt1_phase3) that allows to include gaps for determination of the cbu widths and following the final cbu score, setting gap_weight_opt1_phase3 allows to define different weightings of the gaps in comparison to the subframe stimulation (e.g. the default value of 1.0 means that gaps are considered for the calculation of the cbu widths and concluding for the final cbu score equally to the subframe stimulation whereas a value of 0.5 would mean that gaps are only expected to contribute with 50% in comparison to the subframes stimulation to the final cbu score) 
                      
    end
    
    cbu_score_opt1_mm2 = cbu_width_opt1_mm.*cbu_height_opt1_mm; % final cbu score for option 1 in [MM2] is calculated by multiplying width and height of retinal cbu stimulation 
    cbu_score_opt1_deg2 = cbu_width_opt1_deg.*cbu_height_opt1_deg; % final cbu score for option 1 in [DEG2] is calculated by multiplying width and height of retinal cbu stimulation 
       
    % NON_CBU_SCORE >> gives out the area of retinal stimulation with unbiased content by calculating the surface of unbiased content provocation in [MM2] respectively [DEG2]
    
    non_cbu_height_opt1_mm = subframe_location_table_ver_mm(1,1)-subframe_location_table_ver_mm(2,1); % vertical dimension of the retinal stimulation with unbiased content in [MM]
    non_cbu_height_opt1_deg = subframe_location_table_ver_x1(1,1)-subframe_location_table_ver_x1(2,1); % vertical dimension of the retinal stimulation with unbiased content in [DEG]
    
    if cbu_phase == 1 % horizontal dimension of the retinal stimulation with unbiased content for cbu phase 1 in [MM] respectively [DEG] (consideration of borders between color biased cbu stimulation and unbiased stimulation)
        
        non_cbu_width_opt1_mm = (cbu_table_mm(subframe_number.*2-1,eye_velo_hor_variable_degs_table_index))-(cbu_table_mm(subframe_number.*2,eye_velo_hor_variable_degs_table_index));        
        non_cbu_width_opt1_deg = (cbu_table_x1(subframe_number.*2-1,eye_velo_hor_variable_degs_table_index))-(cbu_table_x1(subframe_number.*2,eye_velo_hor_variable_degs_table_index)); 
        
    elseif cbu_phase == 2 % horizontal dimension of the retinal stimulation with unbiased content for cbu phase 2 in [MM] respectively [DEG] is zero
        
        non_cbu_width_opt1_mm = 0.0;
        non_cbu_width_opt1_deg = 0.0;
        
    elseif cbu_phase == 3 % horizontal dimension of the retinal stimulation with unbiased content for cbu phase 3 in [MM] respectively [DEG] is zero
        
        non_cbu_width_opt1_mm = 0.0;
        non_cbu_width_opt1_deg = 0.0;
        
    end
    
    non_cbu_score_opt1_mm2 = non_cbu_width_opt1_mm.*non_cbu_height_opt1_mm; % final unbiased non-cbu score for option 1 in [MM2] is calculated by multiplying width and height of the unbiased content on the viewer's retina
    non_cbu_score_opt1_deg2 = non_cbu_width_opt1_deg.*non_cbu_height_opt1_deg; % final unbiased non-cbu score for option 1 in [DEG2] is calculated by multiplying width and height of the unbiased content on the viewer's retina
    
    % REFERENCE_SCORE >> gives out the retinal area in [MM2] respectively [DEG2] that would be stimulated with unbiased content if no relative retinal movement was existent 
    
    reference_height_opt1_mm = subframe_location_table_ver_mm(1,1)-subframe_location_table_ver_mm(2,1); % vertical dimension of the retinal stimulation with unbiased content for the reference assumption in [MM]
    reference_height_opt1_deg = subframe_location_table_ver_x1(1,1)-subframe_location_table_ver_x1(2,1); % vertical dimension of the retinal stimulation with unbiased content for the reference assumption in [DEG]
    
    reference_width_opt1_mm = subframe_location_table_hor_mm(1,eye_velo_hor_variable_degs_table_index)-subframe_location_table_hor_mm(2,eye_velo_hor_variable_degs_table_index); % horizontal dimension of the retinal stimulation with unbiased content for the reference assumption in [MM]
    reference_width_opt1_deg = subframe_location_table_hor_x1(1,eye_velo_hor_variable_degs_table_index)-subframe_location_table_hor_x1(2,eye_velo_hor_variable_degs_table_index); % horizontal dimension of the retinal stimulation with unbiased content for the reference assumption in [DEG]
    
    reference_score_opt1_mm2 = reference_width_opt1_mm.*reference_height_opt1_mm; % final reference score in [MM2] is calculated by multiplying width and height of the assumed unbiased content for the defined conditions
    reference_score_opt1_deg2 = reference_width_opt1_deg.*reference_height_opt1_deg; % final reference score in [DEG2] is calculated by multiplying width and height of the assumed unbiased content for the defined conditions
  
% OPTION 2 - area & intensity of retinal stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description

    % CBU_SCORE >> gives out the retinal stimulation by calculating the volume of cbu provokation (inclusion of retinal 2D area in [MM2] respectively [DEG2] and additionally intensity of retinal stimulation) 
    
    cbu_height_opt2_mm = subframe_location_table_ver_mm(1,1)-subframe_location_table_ver_mm(2,1); % vertical dimension of the retinal cbu stimulation in [MM]
    cbu_height_opt2_deg = subframe_location_table_ver_x1(1,1)-subframe_location_table_ver_x1(2,1); % vertical dimension of the retinal cbu stimulation in [DEG]
    
    if cbu_phase == 1 % calculation of the integral of the 3D body front surface that represents cbu stimulation during cbu phase 1 under consideration of [MM] and [DEG] expression for retinal stimulation (consideration of borders between color biased cbu stimulation and unbiased stimulation, see confluence for further description)
      
        for ii = 1:1:size(x_value_table,1) % conversion of X_VALUE_TABLE from [DEG] to [MM]
            
            if x_value_table(ii,1) < 0
                
                x_value_table_x2(ii,1) = (90-(180-(180-(180-abs(x_value_table(ii,1))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(x_value_table(ii,1))-90))))).*(tan((deg2rad(abs(x_value_table(ii,1))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))))).*(-1);
                x_value_table_mm(ii,1) = ((x_value_table_x2(ii,1)./360).*2.*pi.*RR_m).*1000;
                
            else
                
                x_value_table_x2(ii,1) = 90-(180-(180-(180-abs(x_value_table(ii,1))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(x_value_table(ii,1))-90))))).*(tan((deg2rad(abs(x_value_table(ii,1))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))));
                x_value_table_mm(ii,1) = ((x_value_table_x2(ii,1)./360).*2.*pi.*RR_m).*1000;
            
            end
            
        end
        
        border1_mm = cbu_table_mm(subframe_number.*2-1,eye_velo_hor_variable_degs_table_index); % during cbu phase 1 the borders between color biased cbu stimulation and unbiased stimulation have to be considered (refering to [MM] expression)
        border2_mm = cbu_table_mm(subframe_number.*2,eye_velo_hor_variable_degs_table_index);                   
        
        index1_mm = find(flipud(x_value_table_mm)==border1_mm); % indices of border1 and border2 within x_value_table_mm need to be determined (refering to [MM] expression)
        index2_mm = find(flipud(x_value_table_mm)==border2_mm);
        
        cbu_front_surface_mm = trapz(flipud(x_value_table_mm(1:index2_mm,1)),flipud(xyL_table(1:index2_mm,end)))+trapz(flipud(x_value_table_mm(index1_mm:end,1)),flipud(xyL_table(index1_mm:end,end))); % integral of the 3D body front surface is finally calculated under usage of the function "trapz" since the geometrical shapes allow the use without calculation error (refering to [MM] expression)    

        border1_deg = cbu_table_x1(subframe_number.*2-1,eye_velo_hor_variable_degs_table_index); % during cbu phase 1 the borders between color biased cbu stimulation and unbiased stimulation have to be considered (refering to [DEG] expression)
        border2_deg = cbu_table_x1(subframe_number.*2,eye_velo_hor_variable_degs_table_index);                   
        
        index1_deg = find(flipud(x_value_table)==border1_deg); % indices of border1 and border2 within x_value_table need to be determined (refering to [DEG] expression)
        index2_deg = find(flipud(x_value_table)==border2_deg);
        
        cbu_front_surface_deg = trapz(flipud(x_value_table(1:index2_deg,1)),flipud(xyL_table(1:index2_deg,end)))+trapz(flipud(x_value_table(index1_deg:end,1)),flipud(xyL_table(index1_deg:end,end))); % integral of the 3D body front surface is finally calculated under usage of the function "trapz" since the geometrical shapes allow the use without calculation error (refering to [DEG] expression)
               
    elseif cbu_phase == 2 % calculation of the integral of the 3D body front surface that represents cbu stimulation during cbu phase 2 under consideration of [MM] and [DEG] expression for retinal stimulation (no consideration of borders between color biased cbu stimualation and unbiased stimulation since there is no more unbiased content within cbu phase 2, see confluence for further description)
        
        for ii = 1:1:size(x_value_table,1) % conversion of X_VALUE_TABLE from [DEG] to [MM]
            
            if x_value_table(ii,1) < 0
                
                x_value_table_x2(ii,1) = (90-(180-(180-(180-abs(x_value_table(ii,1))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(x_value_table(ii,1))-90))))).*(tan((deg2rad(abs(x_value_table(ii,1))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))))).*(-1);
                x_value_table_mm(ii,1) = ((x_value_table_x2(ii,1)./360).*2.*pi.*RR_m).*1000;
                
            else
                
                x_value_table_x2(ii,1) = 90-(180-(180-(180-abs(x_value_table(ii,1))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(x_value_table(ii,1))-90))))).*(tan((deg2rad(abs(x_value_table(ii,1))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))));
                x_value_table_mm(ii,1) = ((x_value_table_x2(ii,1)./360).*2.*pi.*RR_m).*1000;
            
            end
            
        end
          
        cbu_front_surface_mm = trapz(flipud(x_value_table_mm(1:end,1)),flipud(xyL_table(1:end,end))); % integral of the 3D body front surface is finally calculated under usage of the function "trapz" since the geometrical shapes allow the use without calculation error (refering to [MM] expression)
        cbu_front_surface_deg = trapz(flipud(x_value_table(1:end,1)),flipud(xyL_table(1:end,end))); % integral of the 3D body front surface is finally calculated under usage of the function "trapz" since the geometrical shapes allow the use without calculation error (refering to [DEG] expression)       
        
    elseif cbu_phase == 3 % calculation of the integral of the 3D body front surface that represents retinal cbu stimulation during cbu phase 3 under consideration of [MM] and [DEG] expression for retinal stimulation
        
        for ii = 1:1:size(x_value_table,1) % conversion of X_VALUE_TABLE from [DEG] to [MM]
            
            if x_value_table(ii,1) < 0
                
                x_value_table_x2(ii,1) = (90-(180-(180-(180-abs(x_value_table(ii,1))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(x_value_table(ii,1))-90))))).*(tan((deg2rad(abs(x_value_table(ii,1))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))))).*(-1);
                x_value_table_mm(ii,1) = ((x_value_table_x2(ii,1)./360).*2.*pi.*RR_m).*1000;
                
            else
                
                x_value_table_x2(ii,1) = 90-(180-(180-(180-abs(x_value_table(ii,1))-90))-(rad2deg(asin(sin((deg2rad((180-(180-abs(x_value_table(ii,1))-90))))).*(tan((deg2rad(abs(x_value_table(ii,1))).*(CV_F2_relax_m-CV_N2_relax_m-RR_m))))/RR_m))));
                x_value_table_mm(ii,1) = ((x_value_table_x2(ii,1)./360).*2.*pi.*RR_m).*1000;
            
            end
            
        end
        
        gap_weight_opt2_phase3 = 1.0; % calculation of the 3D body front surface during cbu phase 3 must include an approach that defines how to handle the gaps between the retinal subframe stimulation, at present the gaps are fully integrated into the determination of the 3D body front surface which means that the gaps have the same weight as the subframe stimulation (weighting factor is set to 1.0 by default, see confluence for further description)
        
        aa = 0; % calculation on basis of [MM] expression
        
        for jj = 1:1:subframe_number   
            
            aa = aa + cbu_table_mm(jj.*4-3,eye_velo_hor_variable_degs_table_index)-cbu_table_mm(jj.*4-2,eye_velo_hor_variable_degs_table_index); % aa represents the horizontal share of retinal cbu stimulation by the subframes 
        
        end
        
        bb = (cbu_table_mm(1,eye_velo_hor_variable_degs_table_index)-cbu_table_mm(end,eye_velo_hor_variable_degs_table_index))-aa; % bb represents the horizontal share of gaps between the retinal subframe stimulation
        
        cbu_front_surface_mm = trapz(flipud(x_value_table_mm(1:end,1)),flipud(xyL_table(1:end,end))).*(1+(bb/aa).*gap_weight_opt2_phase3); % calculation of the 3D body front surface of (spatially separated) subframes with function "trapz" (gaps do no contribute to this volume calculation since within the gaps there is no stimulation with light), however, it is assumed that gaps have an effect on the perceived cbu effect, additional multiplication with factor (1+(bb/aa).*gap_weight_opt2_phase3) that allows to include gaps for determination of the 3D body front surface and following the final cbu score, setting gap_weight_opt2_phase3 allows to define different weightings of the gaps in comparison to the subframe stimulation (e.g. the default value of 1.0 means that gaps are considered for the calculation of the 3D body front surface and concluding for the final cbu score equally to the subframe stimulation whereas a value of 0.5 would mean that gaps are only expected to contribute with 50% in comparison to the subframes stimulation to the final cbu score)   

        aa = 0; % calculation on basis of [DEG] expression
        
        for jj = 1:1:subframe_number   
            
            aa = aa + cbu_table_x1(jj.*4-3,eye_velo_hor_variable_degs_table_index)-cbu_table_x1(jj.*4-2,eye_velo_hor_variable_degs_table_index); % aa represents the horizontal share of retinal cbu stimulation by the subframes 
        
        end
        
        bb = (cbu_table_x1(1,eye_velo_hor_variable_degs_table_index)-cbu_table_x1(end,eye_velo_hor_variable_degs_table_index))-aa; % bb represents the horizontal share of gaps between the retinal subframe stimulation
        
        cbu_front_surface_deg = trapz(flipud(x_value_table(1:end,1)),flipud(xyL_table(1:end,end))).*(1+(bb/aa).*gap_weight_opt2_phase3); % calculation of the 3D body front surface of (spatially separated) subframes with function "trapz" (gaps do no contribute to this volume calculation since within the gaps there is no stimulation with light), however, it is assumed that gaps have an effect on the perceived cbu effect, additional multiplication with factor (1+(bb/aa).*gap_weight_opt2_phase3) that allows to include gaps for determination of the 3D body front surface and following the final cbu score, setting gap_weight_opt2_phase3 allows to define different weightings of the gaps in comparison to the subframe stimulation (e.g. the default value of 1.0 means that gaps are considered for the calculation of the 3D body front surface and concluding for the final cbu score equally to the subframe stimulation whereas a value of 0.5 would mean that gaps are only expected to contribute with 50% in comparison to the subframes stimulation to the final cbu score)
                
    end
    
    % calculation of frame number that is not rounded (as factor for scaling the scores)    
    
    if cat == 2 % final cbu_score_opt2 (also non_cbu_score_opt2 and reference_score_opt2, see below) has to be calculated by multiplying with frame_number_not_rounded in order to make conditions with variable frame rates comparable by scaling the cbu score up (see confluence for detailed description and explanation)
    
        frame_number_not_rounded = (eye_distance_hor_px./eye_velo_hor_pxframe)+1; % frame number within CAT2 is calculated on basis of the time eye movement is executed (since no content movement is executed at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; calculated frame number is not rounded off (like it is done in the section DETERMINATION OF EXPERIMENTAL CATEGORY) in order to calculate the exact multiplication factor

    elseif cat == 3
        
        if (content_distance_hor_px/content_velo_hor_pxframe) >= (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is equal or longer than the time of eye movement (eye_distance_px/eye_velo_pxframe) than the frame number to be considered is calculated on basis of the time eye movement is executed (even when content movement is executed for a longer period of time eye movement time is still basis for calculation since CBU only occurs when eye movement is in progress, otherwise it will be CAT1 (eye fix, content variable) that does not provoke CBU at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; calculated frame number is not rounded off (like it is done in the section DETERMINATION OF EXPERIMENTAL CATEGORY) in order to calculate the exact multiplication factor; statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning) 
        
            frame_number_not_rounded = (eye_distance_hor_px./eye_velo_hor_pxframe)+1;
        
        elseif (content_distance_hor_px/content_velo_hor_pxframe) < (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is shorter than the time of eye movement (eye_distance_px/eye_velo_pxframe) than the frame number is calculated on basis of the time content movement is executed (since the content disappears directly after reaching the stop point of its movement path, without presentation of content no CBU can occur); calculation of time in [FRAME] by the term (content_distance_px./content_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; calculated frame number is not rounded off (like it is done in the section DETERMINATION OF EXPERIMENTAL CATEGORY) in order to calculate the exact multiplication factor (in this case the frame number is generally expected to be an integer since content distance and content movement velocity are expected to be positive integers); statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning)
            
            frame_number_not_rounded = (content_distance_hor_px./content_velo_hor_pxframe)+1;
    
        end
        
    end
       
    
    cbu_score_opt2_mm2 = cbu_front_surface_mm.*cbu_height_opt2_mm.*frame_number_not_rounded; % final cbu score for option 2 (refering to [MM2] expression of retinal stimulation) is calculated by multiplying 3D body front surface and height of retinal cbu stimulation (and also frame number, see explanation in code description above)
    cbu_score_opt2_deg2 = cbu_front_surface_deg.*cbu_height_opt2_deg.*frame_number_not_rounded; % final cbu score for option 2 (refering to [DEG2] expression of retinal stimulation) is calculated by multiplying 3D body front surface and height of retinal cbu stimulation (and also frame number, see explanation in code description above)
    
    % NON_CBU_SCORE >> gives out the value for retinal stimulation with unbiased content by calculating the volume of unbiased content provocation (refering to [MM2] respectively [DEG2] expression for retinal stimulation)
    
    non_cbu_height_opt2_mm = subframe_location_table_ver_mm(1,1)-subframe_location_table_ver_mm(2,1); % vertical dimension of the retinal stimulation with unbiased content in [MM]
    non_cbu_height_opt2_deg = subframe_location_table_ver_x1(1,1)-subframe_location_table_ver_x1(2,1); % vertical dimension of the retinal stimulation with unbiased content in [DEG]
    
    if cbu_phase == 1 % integral of the 3D body front surface (calculated with "trapz" function) of unbiased content provocation for cbu phase 1 under consideration of [MM] and [DEG] expression for retinal stimulation (borders between color biased cbu stimulation and unbiased stimulation have to be considered - see index1 and index2 - both indices have been calculated during calculation of cbu score for option 2 above)
          
        non_cbu_front_surface_mm = trapz(flipud(x_value_table_mm(index2_mm:index1_mm,1)),flipud(xyL_table(index2_mm:index1_mm,end)));
        non_cbu_front_surface_deg = trapz(flipud(x_value_table(index2_deg:index1_deg,1)),flipud(xyL_table(index2_deg:index1_deg,end)));
        
    elseif cbu_phase == 2 % no retinal stimulation with unbiased content within cbu phase 2 under consideration of [MM] and [DEG] expression for retinal stimulation
        
        non_cbu_front_surface_mm = 0.0;
        non_cbu_front_surface_deg = 0.0;
        
    elseif cbu_phase == 3 % no retinal stimulation with unbiased content within cbu phase 3 under consideration of [MM] and [DEG] expression for retinal stimulation
        
        non_cbu_front_surface_mm = 0.0;
        non_cbu_front_surface_deg = 0.0;
    
    end
    
    non_cbu_score_opt2_mm2 = non_cbu_front_surface_mm.*non_cbu_height_opt2_mm.*frame_number_not_rounded; % final non-cbu score for retinal stimulation with unbiased content (refering to [MM2] expression of retinal stimulation) is determined by multiplying 3D body front surface and height (and also frame number, see explanation in code description above)
    non_cbu_score_opt2_deg2 = non_cbu_front_surface_deg.*non_cbu_height_opt2_deg.*frame_number_not_rounded; % final non-cbu score for retinal stimulation with unbiased content (refering to [DEG2] expression of retinal stimulation) is determined by multiplying 3D body front surface and height (and also frame number, see explanation in code description above)     
    
    % REFERENCE_SCORE >> gives out the volume of stimulation for unbiased content if no relative retinal movement was existent (refering to [MM2] respectively [DEG2] expression for retinal stimulation)
    
    reference_height_opt2_mm = subframe_location_table_ver_mm(1,1)-subframe_location_table_ver_mm(2,1); % vertical dimension of the retinal stimulation with unbiased content under the made reference assumption in [MM]
    reference_height_opt2_deg = subframe_location_table_ver_x1(1,1)-subframe_location_table_ver_x1(2,1); % vertical dimension of the retinal stimulation with unbiased content under the made reference assumption in [DEG]
    
    reference_front_surface_mm = (subframe_location_table_hor_mm(1,eye_velo_hor_variable_degs_table_index)-subframe_location_table_hor_mm(2,eye_velo_hor_variable_degs_table_index)).*(subframe_on_sec.*sum(sf_luminance)); % 3D body front surface of unbiased content provocation under the made reference assumption (refering to [MM] expression for retinal stimulation)
    reference_front_surface_deg = (subframe_location_table_hor_x1(1,eye_velo_hor_variable_degs_table_index)-subframe_location_table_hor_x1(2,eye_velo_hor_variable_degs_table_index)).*(subframe_on_sec.*sum(sf_luminance)); % 3D body front surface of unbiased content provocation under the made reference assumption (refering to [DEG] expression for retinal stimulation)
    
    reference_score_opt2_mm2 = reference_front_surface_mm.*reference_height_opt2_mm.*frame_number_not_rounded; % final reference score (refering to [MM2] expression for retinal stimulation) is calculated by multiplying 3D body front surface and height of the assumed unbiased content for the defined conditions (and also frame number, see explanation in code description above)
    reference_score_opt2_deg2 = reference_front_surface_deg.*reference_height_opt2_deg.*frame_number_not_rounded; % final reference score (refering to [DEG2] expression for retinal stimulation) is calculated by multiplying 3D body front surface and height of the assumed unbiased content for the defined conditions (and also frame number, see explanation in code description above)
    
% OPTION 3 - area & intensity & color of retinal stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description
    
    % CBU_SCORE >> ...
    
    % NON_CBU_SCORE >> ...
    
    % REFERENCE_SCORE >> ...
    
% inclusion of the final model output into the LOOP_TABLE (cbu score, non-cbu score, reference score, both options 1 and 2, consideration of [DEG] values)   
    
    loop_table(loop_count,7)= cbu_score_opt1_deg2;
    loop_table(loop_count,8)= non_cbu_score_opt1_deg2;
    loop_table(loop_count,9)= reference_score_opt1_deg2;
    loop_table(loop_count,10)= cbu_score_opt2_deg2;
    loop_table(loop_count,11)= non_cbu_score_opt2_deg2;
    loop_table(loop_count,12)= reference_score_opt2_deg2;
    
end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% OVERVIEW TABLE >> lists the loop function's input (e.g. subframe number, duty cycle, content width, frame rate, eye/content movement velocity) and it's output (cbu score, non-cbu score, reference score) for every run of the loop (loop number)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
overview_table = array2table(loop_table);
overview_table.Properties.VariableNames(1:12) = {'SUBFRAME_NUMBER','DUTY_CYCLE','CONTENT_WIDTH_PX','FRAME_RATE_HZ', 'EYE_MOVEMENT_VELOCITY_PXFR', 'CONTENT_MOVEMENT_VELOCITY_PXFR', 'CBU_SCORE_OPT1_DEG2', 'NON_CBU_SCORE_OPT1_DEG2', 'REFERENCE_SCORE_OPT1_DEG2', 'CBU_SCORE_OPT2', 'NON_CBU_SCORE_OPT2', 'REFERENCE_SCORE_OPT2'};

overview_table_add = [1:loop_count]';
overview_table_add = array2table(overview_table_add);
overview_table_add.Properties.VariableNames(1:1) = {'LOOP_NUMBER'};

overview_table = [overview_table_add overview_table]
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% THRESHOLD CALCULATION >> gives out the values for the variable parameter (see VARIABLE above) for which a phase transition from phase 1 to phase 2 and from phase 2 to phase 3 is accomplished
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
frame_rate_hz_threshold_ph1_ph2 = (max(eye_velo_hor_variable_degs_table).*(2-duty_cycle))./(content_width_variable_deg_table(eye_velo_hor_variable_degs_table_index,1).*subframe_number);
frame_rate_hz_threshold_ph2_ph3 = (max(eye_velo_hor_variable_degs_table).*(1-duty_cycle))./(content_width_variable_deg_table(eye_velo_hor_variable_degs_table_index,1).*subframe_number);
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (1) >> illustrates the dependency between the variable parameter (x-axis) and cbu score, non-cbu score and reference score as target parameter (y-axis)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
close all; 

f1 = figure (1); hold on; set(gca,'FontSize',15); % figure 1 refers to model output calculation with option 1 

p(1) = plot(loop_table(1:end,4),loop_table(1:end,7), '--o','Color',[0.75 0.0 0.0],'LineWidth', 1.0,'MarkerEdgeColor',[0.75 0.0 0.0],'MarkerFaceColor','w','MarkerSize',6.0); % plot cbu score
p(2) = plot(loop_table(1:end,4),loop_table(1:end,8), '--o','Color',[0.0 0.5 0.0],'LineWidth', 1.0,'MarkerEdgeColor',[0.0 0.5 0.0],'MarkerFaceColor','w','MarkerSize',6.0); % plot non-cbu score
p(3) = plot(loop_table(1:end,4),loop_table(1:end,9), ':ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot reference score

p(4) = line([frame_rate_hz_threshold_ph1_ph2 frame_rate_hz_threshold_ph1_ph2],[0 max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])+max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])/10],'Color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth', 1.0); % plot phase transition threshold (phase 1 to phase 2)
p(5) = line([frame_rate_hz_threshold_ph2_ph3 frame_rate_hz_threshold_ph2_ph3],[0 max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])+max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])/10],'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth', 1.0); % plot phase transition threshold (phase 2 to phase 3)

xlim([frame_rate_hz_min frame_rate_hz_max]);
set(gca,'xTick',frame_rate_hz_min:(frame_rate_hz_max-frame_rate_hz_min)/13:frame_rate_hz_max);
xlabel('FRAME RATE [HZ]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])+max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])/10]);
set(gca,'yTick',0:max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])/10:max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])+max([loop_table(1:end,7);loop_table(1:end,8);loop_table(1:end,9)])/10);
ylabel('MODEL INDEX [DEG2]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('MODEL INDICES IN DEPENDENY OF FRAME RATE (OPTION 1)', 'Fontsize', 21);
legend('theoretical CBU index (model, option 1)','theoretical non-CBU index (model, option 1)', 'theoretical reference index (model, option 1)','CBU phase transition threshold (PH1/PH2)','CBU phase transition threshold (PH2/PH3)','Location','east','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f2 = figure (2); hold on; set(gca,'FontSize',15); % figure 2 refers to model output calculation with option 2

p(1) = plot(loop_table(1:end,4),loop_table(1:end,10), '--o','Color',[0.75 0.0 0.0],'LineWidth', 1.0,'MarkerEdgeColor',[0.75 0.0 0.0],'MarkerFaceColor','w','MarkerSize',6.0); % plot cbu score
p(2) = plot(loop_table(1:end,4),loop_table(1:end,11), '--o','Color',[0.0 0.5 0.0],'LineWidth', 1.0,'MarkerEdgeColor',[0.0 0.5 0.0],'MarkerFaceColor','w','MarkerSize',6.0); % plot non-cbu score
p(3) = plot(loop_table(1:end,4),loop_table(1:end,12), ':ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot reference score

p(4) = line([frame_rate_hz_threshold_ph1_ph2 frame_rate_hz_threshold_ph1_ph2],[0 max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])+max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])/10],'Color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth', 1.0); % plot phase transition threshold (phase 1 to phase 2)
p(5) = line([frame_rate_hz_threshold_ph2_ph3 frame_rate_hz_threshold_ph2_ph3],[0 max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])+max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])/10],'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth', 1.0); % plot phase transition threshold (phase 2 to phase 3)

xlim([frame_rate_hz_min frame_rate_hz_max]);
set(gca,'xTick',frame_rate_hz_min:(frame_rate_hz_max-frame_rate_hz_min)/13:frame_rate_hz_max);
xlabel('FRAME RATE [HZ]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])+max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])/10]);
set(gca,'yTick',0:max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])/10:max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])+max([loop_table(1:end,10);loop_table(1:end,11);loop_table(1:end,12)])/10);
ylabel('MODEL INDEX [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('MODEL INDICES IN DEPENDENY OF FRAME RATE  (OPTION 2)', 'Fontsize', 21);
legend('theoretical CBU index (model, option 2)','theoretical non-CBU index (model, option 2)', 'theoretical reference index (model, option 2)','CBU phase transition threshold (PH1/PH2)','CBU phase transition threshold (PH2/PH3)','Location','east','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%{ 
% ... NORMALIZATION PROCESS DEACTIVATED ...

% CONDITION (2) >> for normalization process 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if frame_rate_hz_min ~= 30.0 || frame_rate_hz_max ~= 420.0 || ((frame_rate_hz_max-frame_rate_hz_min)/frame_rate_hz_step) ~= fix((frame_rate_hz_max-frame_rate_hz_min)/frame_rate_hz_step)
    
    cprintf('err','ATTENTION - frame rate range not adequate for normalization process (see VARIABLE), full frame rate range from 30.0 to 420.0Hz necessary \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% MIN-MAX-NORMALIZATION OF DATA >> normalization of theoretical model cbu data (empirical data is taken from data set without any alteration, see below)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if round(max(eye_velo_hor_variable_degs_table)) == 18.0 

    emp_data_cat3_table_normalization = [30,3.38710000000000,0.852800000000000;60,2.23580000000000,0.715600000000000;90,1.73430000000000,0.639700000000000;120,1.49510000000000,0.538600000000000;150,1.28360000000000,0.402700000000000;180,1.22680000000000,0.374800000000000;210,1.13240000000000,0.287000000000000;240,1.12550000000000,0.269500000000000;300,1.11540000000000,0.276000000000000;360,1.08640000000000,0.235400000000000;420,1.08060000000000,0.233000000000000];
    
elseif round(max(eye_velo_hor_variable_degs_table)) == 36.0
    
    emp_data_cat3_table_normalization = [30,3.84740000000000,0.804400000000000;60,2.87780000000000,0.784900000000000;90,2.41140000000000,0.742300000000000;120,1.97880000000000,0.678900000000000;150,1.74320000000000,0.618800000000000;180,1.53980000000000,0.546500000000000;210,1.42350000000000,0.503300000000000;240,1.34710000000000,0.417100000000000;300,1.26230000000000,0.443200000000000;360,1.18920000000000,0.329600000000000;420,1.18870000000000,0.332400000000000];
    
elseif round(max(eye_velo_hor_variable_degs_table)) == 54.0
    
    emp_data_cat3_table_normalization = [30,4.12080000000000,0.762200000000000;60,3.31650000000000,0.792000000000000;90,2.74910000000000,0.812000000000000;120,2.33010000000000,0.772000000000000;150,1.98550000000000,0.708200000000000;180,1.83170000000000,0.666500000000000;210,1.59910000000000,0.607900000000000;240,1.52240000000000,0.592900000000000;300,1.31770000000000,0.416500000000000;360,1.31480000000000,0.459100000000000;420,1.22050000000000,0.383600000000000];
                
else
    
    error('ATTENTION - there is no empirical data the theoretical model output can be compared with for the entered value of horizontal eye movement velocity');
    
end
    
lower_norm_limit = min(emp_data_cat3_table_normalization(1:end,2));
upper_norm_limit = max(emp_data_cat3_table_normalization(1:end,2));

% normalization process of theoretical cbu model data refering to option 1 (see confluence)

cbu_score_opt1_deg2_mean_normalized = ((loop_table(1:end,7)-min(loop_table(1:end,7)))./(max(loop_table(1:end,7))-min(loop_table(1:end,7)))).*(upper_norm_limit-lower_norm_limit)+lower_norm_limit;

cbu_score_opt1_deg2_variance_normalized = max(cbu_score_opt1_deg2_mean_normalized).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4);

cbu_score_opt1_deg2_table_normalized = [loop_table(1:end,4),cbu_score_opt1_deg2_mean_normalized,cbu_score_opt1_deg2_variance_normalized];

% normalization process of theoretical cbu model data refering to option 2 (see confluence)

cbu_score_opt2_deg2_mean_normalized = ((loop_table(1:end,10)-min(loop_table(1:end,10)))./(max(loop_table(1:end,10))-min(loop_table(1:end,10)))).*(upper_norm_limit-lower_norm_limit)+lower_norm_limit;

cbu_score_opt2_deg2_variance_normalized = max(cbu_score_opt2_deg2_mean_normalized).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4);

cbu_score_opt2_deg2_table_normalized = [loop_table(1:end,4),cbu_score_opt2_deg2_mean_normalized,cbu_score_opt2_deg2_variance_normalized];
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% EXPONENTIAL FUNCTION >> generates exponential function for curve fitting regarding empirical cbu data values
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
x = emp_data_cat3_table_normalization(1:end,1);
y = emp_data_cat3_table_normalization(1:end,2);
    
genexp = @(b,x) b(1).*exp(b(2).*(x(:,1)))+b(3);
b0 = [0, 0, 0];
framepred = linspace(0, 840, 1000)';
    
mod1 = fitnlm(x, y, genexp, b0);
cbupred = predict(mod1,framepred);
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (2) >> comparison of empirical cbu data with theoretical model cbu data (after normalization of the theoretical model cbu data)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
f3 = figure (3); hold on; set(gca,'FontSize',15); % figure 3 refers to cbu score calculation with option 1

p(1) = plot(cbu_score_opt1_deg2_table_normalized(1:end,1),cbu_score_opt1_deg2_table_normalized(1:end,2),'--','Color',[0.9 0.0 0.0],'LineWidth',1.0); % plot theoretical cbu score data (only mean value)
% p(1) = boundedline(cbu_score_opt1_deg2_table_normalized(1:end,1),cbu_score_opt1_deg2_table_normalized(1:end,2),cbu_score_opt1_deg2_table_normalized(1:end,3),'r','alpha'); % plot theoretical cbu score data (mean and spread values)

p(2) = plot(framepred,cbupred,'b--','LineWidth',1.0); % plot empirical cbu score data (curve fitted exponential function)
p(3) = plot(x,y,'o','MarkerFaceColor','w','MarkerEdgeColor','b','LineWidth',1.0,'MarkerSize',4.0); % plot empirical cbu score data (mean values)
p(4) = boundedline(emp_data_cat3_table_normalization(1:end,1), emp_data_cat3_table_normalization(1:end,2), emp_data_cat3_table_normalization(1:end,3),':','alpha'); % plot empirical cbu score data (standard deviation)

p(5) = line([frame_rate_hz_threshold_ph1_ph2 frame_rate_hz_threshold_ph1_ph2],[1 5],'Color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth', 1.0); % plot phase threshold (phase 1 to phase 2)
p(6) = line([frame_rate_hz_threshold_ph2_ph3 frame_rate_hz_threshold_ph2_ph3],[1 5],'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth', 1.0); % plot phase threshold (phase 2 to phase 3)

xlim([frame_rate_hz_min frame_rate_hz_max]);
set(gca,'xTick',frame_rate_hz_min:(frame_rate_hz_max-frame_rate_hz_min)/13:frame_rate_hz_max);
xlabel('FRAME RATE [HZ]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([1 5]);
set(gca,'yTick',1:0.5:5);
ylabel('NORMALIZED CBU [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('EMPIRICAL VS. THEORETICAL CBU (OPTION 1)', 'Fontsize', 21);
legend('theoretical CBU index (model, option 1)','empirical CBU score (fitted exponential curve)', 'empirical CBU score (MN)', 'empirical CBU score (SD)', 'empirical CBU score (linear interpolation)','CBU phase transition threshold (PH1/PH2)','CBU phase transition threshold (PH2/PH3)','Fontsize',9);

box off;  
set(gcf,'WindowStyle','docked');
hold off;

f4 = figure (4); hold on; set(gca,'FontSize',15); % figure 4 refers to cbu score calculation with option 2

p(1) = plot(cbu_score_opt2_deg2_table_normalized(1:end,1),cbu_score_opt2_deg2_table_normalized(1:end,2),'--','Color',[0.9 0.0 0.0],'LineWidth',1.0); % plot theoretical cbu score data (only mean value)
% p(1) = boundedline(cbu_score_opt2_deg2_table_normalized(1:end,1),cbu_score_opt2_deg2_table_normalized(1:end,2),cbu_score_opt2_deg2_table_normalized(1:end,3),'g','alpha'); % plot theoretical cbu score data (mean and spread values)

p(2) = plot(framepred,cbupred,'b--','LineWidth',1.0); % plot empirical cbu score data (curve fitted exponential function)
p(3) = plot(x,y,'o','MarkerFaceColor','w','MarkerEdgeColor','b','LineWidth',1.0,'MarkerSize',4.0); % plot empirical cbu score data (mean values)
p(4) = boundedline(emp_data_cat3_table_normalization(1:end,1), emp_data_cat3_table_normalization(1:end,2), emp_data_cat3_table_normalization(1:end,3),':','alpha'); % plot empirical cbu score data (standard deviation)

p(5) = line([frame_rate_hz_threshold_ph1_ph2 frame_rate_hz_threshold_ph1_ph2],[1 5],'Color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth', 1.0); % plot phase threshold (phase 1 to phase 2)
p(6) = line([frame_rate_hz_threshold_ph2_ph3 frame_rate_hz_threshold_ph2_ph3],[1 5],'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth', 1.0); % plot phase threshold (phase 2 to phase 3)

xlim([frame_rate_hz_min frame_rate_hz_max]);
set(gca,'xTick',frame_rate_hz_min:(frame_rate_hz_max-frame_rate_hz_min)/13:frame_rate_hz_max);
xlabel('FRAME RATE [HZ]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([1 5]);
set(gca,'yTick',1:0.5:5);
ylabel('NORMALIZED CBU [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('EMPIRICAL VS. THEORETICAL CBU (OPTION 2)', 'Fontsize', 21);
legend('theoretical CBU index (model, option 2)','empirical CBU score (fitted exponential curve)', 'empirical CBU score (MN)', 'empirical CBU score (SD)', 'empirical CBU score (linear interpolation)','CBU phase transition threshold (PH1/PH2)','CBU phase transition threshold (PH2/PH3)','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% ... NORMALIZATION PROCESS DEACTIVATED ...
%} 

% CONDITION (3) >> for standardization process 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if frame_rate_hz_min ~= 30.0 || frame_rate_hz_max ~= 420.0 || frame_rate_hz_step ~= 5.0
    
    cprintf('err','ATTENTION - frame rate range and steps not adequate for standardization process (see VARIABLE), full frame rate range from 30.0 to 420.0Hz and frame rate steps of 5.0Hz necessary \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% STANDARDIZATION OF DATA >> standardization of both data sets - theoretical cbu model data set and empirical cbu data set
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if round(max(eye_velo_hor_variable_degs_table)) == 18.0

    emp_data_cat3_standardization = [30,3.22170000000000,4.06330000000000,2.35170000000000,4.71670000000000,2.95500000000000,3.76170000000000,4.41170000000000,2.53670000000000,2.81830000000000,3.28170000000000,3.43830000000000,2.04830000000000,2.34170000000000,3.71170000000000,2.84170000000000,3.37000000000000,2.70500000000000,2.83170000000000,3.07670000000000,3.65330000000000,2.77830000000000,4.36830000000000,3.95000000000000,3.96830000000000,4.01670000000000,4.48500000000000,2.76000000000000,2.79170000000000,4.10500000000000,4.25170000000000,3.38710000000000;60,2.02670000000000,3.45000000000000,1.41170000000000,2.04500000000000,1.90330000000000,3.02830000000000,1.73830000000000,1.95000000000000,1.94000000000000,2.72670000000000,2.19830000000000,1.65670000000000,1.80500000000000,1.66670000000000,1.94170000000000,1.89330000000000,1.66000000000000,1.51830000000000,2.18500000000000,2.10830000000000,1.98000000000000,2.86330000000000,3.08330000000000,2.05500000000000,3.14670000000000,3.26330000000000,2.05500000000000,1.66170000000000,3.19170000000000,2.92000000000000,2.23580000000000;90,1.66830000000000,2.45330000000000,1.19830000000000,1.37330000000000,1.51830000000000,2.09670000000000,1.55170000000000,1.25830000000000,1.51000000000000,2.18000000000000,1.78830000000000,1.57170000000000,1.52330000000000,1.15330000000000,1.71330000000000,1.28830000000000,1.21670000000000,1.16000000000000,1.90670000000000,1.20170000000000,1.42170000000000,2.23000000000000,2.09830000000000,1.50830000000000,2.69170000000000,2.73500000000000,1.48670000000000,1.32830000000000,2.93000000000000,2.26670000000000,1.73430000000000;120,1.37500000000000,2.17830000000000,1.06170000000000,1,1.41330000000000,2.03170000000000,1.29170000000000,1.16330000000000,1.25500000000000,1.94330000000000,1.74830000000000,1.25670000000000,1.58170000000000,1.03170000000000,1.52670000000000,1.28500000000000,1.11670000000000,1.15830000000000,1.75330000000000,1.07170000000000,1.26500000000000,1.89500000000000,1.40330000000000,1.32670000000000,1.97170000000000,2.37500000000000,1.12670000000000,1.14830000000000,2.16670000000000,1.93000000000000,1.49510000000000;150,1.05670000000000,1.37830000000000,1.02830000000000,1,1.34330000000000,1.68000000000000,1,1.08830000000000,1.48670000000000,1.56170000000000,1.33170000000000,1.16500000000000,1.52670000000000,1.11000000000000,1.14330000000000,1.10830000000000,1.00830000000000,1.07330000000000,1.54330000000000,1.04170000000000,1.11170000000000,1.78670000000000,1.11830000000000,1.13670000000000,1.75670000000000,2.08000000000000,1,1.08000000000000,1.30330000000000,1.45830000000000,1.28360000000000;180,1.06170000000000,1.32000000000000,1.00670000000000,1,1.13670000000000,1.50000000000000,1,1.01830000000000,1.37670000000000,1.46000000000000,1.22000000000000,1.06670000000000,1.57830000000000,1.00830000000000,1.20670000000000,1.17830000000000,1.00830000000000,1.05500000000000,1.16330000000000,1.03670000000000,1.05830000000000,1.63670000000000,1.05330000000000,1.07000000000000,1.44830000000000,2.14000000000000,1,1.03330000000000,1.56670000000000,1.39500000000000,1.22680000000000;210,1,1.08830000000000,1.00330000000000,1,1.02170000000000,1.22500000000000,1,1.01330000000000,1.19500000000000,1.45170000000000,1.15170000000000,1.10000000000000,1.35670000000000,1.00330000000000,1.02000000000000,1.04000000000000,1,1.00170000000000,1.21330000000000,1.02830000000000,1.01170000000000,1.23500000000000,1.06670000000000,1.12500000000000,1.16830000000000,2.01000000000000,1,1,1.10830000000000,1.33500000000000,1.13240000000000;240,1.00830000000000,1.03670000000000,1.05830000000000,1.01500000000000,1.06670000000000,1.22330000000000,1.01500000000000,1.01830000000000,1.25830000000000,1.07830000000000,1.12000000000000,1.05330000000000,1.38500000000000,1.00670000000000,1,1.03000000000000,1.00670000000000,1.05330000000000,1.01830000000000,1.01500000000000,1,1.14670000000000,1.04170000000000,1.22830000000000,1.24830000000000,2.09330000000000,1,1,1.20330000000000,1.33670000000000,1.12550000000000;300,1.00830000000000,1.29670000000000,1.00170000000000,1.01000000000000,1.01330000000000,1.08830000000000,1,1,1.39670000000000,1.06500000000000,1.15330000000000,1.06830000000000,1.16500000000000,1.02000000000000,1,1.13170000000000,1.01330000000000,1.01170000000000,1.13670000000000,1.00830000000000,1,1.20170000000000,1.03170000000000,1.13330000000000,1.27000000000000,2.06830000000000,1,1,1,1.16830000000000,1.11540000000000;360,1.00500000000000,1.07330000000000,1.01330000000000,1,1.08330000000000,1.06500000000000,1,1.04670000000000,1.37500000000000,1.13170000000000,1,1.00830000000000,1.04670000000000,1,1,1,1.02670000000000,1,1.08000000000000,1.00500000000000,1.01670000000000,1.17000000000000,1.04670000000000,1.01500000000000,1.03670000000000,2.10830000000000,1,1,1.16830000000000,1.07000000000000,1.08640000000000;420,1.02170000000000,1,1,1.09500000000000,1,1.05500000000000,1,1,1.23830000000000,1.07670000000000,1.06330000000000,1.00500000000000,1.29500000000000,1,1,1.02670000000000,1.00830000000000,1.00330000000000,1.06670000000000,1.00330000000000,1,1.24500000000000,1.02830000000000,1.04170000000000,1.01330000000000,2.09500000000000,1,1.00330000000000,1.02170000000000,1.01170000000000,1.08060000000000];
       
elseif round(max(eye_velo_hor_variable_degs_table)) == 36.0
    
    emp_data_cat3_standardization = [30,3.90500000000000,4.20830000000000,3.05330000000000,4.84170000000000,3.92170000000000,4.32170000000000,4.13500000000000,3.41170000000000,3.16830000000000,3.53670000000000,3.86330000000000,2.44330000000000,2.37170000000000,4.11500000000000,3.61000000000000,4.23170000000000,3.30500000000000,3.87830000000000,3.35830000000000,4.37500000000000,3.22330000000000,4.69500000000000,4.17000000000000,3.99670000000000,4.38330000000000,4.60330000000000,3.41330000000000,3.40830000000000,4.56000000000000,4.91500000000000,3.84740000000000;60,2.91500000000000,3.93830000000000,1.88830000000000,3.84500000000000,2.68500000000000,3.59670000000000,2.67170000000000,2.23170000000000,2.24670000000000,3.05830000000000,3.13670000000000,1.94500000000000,1.99170000000000,2.78170000000000,2.69000000000000,3.13830000000000,2.20330000000000,2.36830000000000,2.64670000000000,2.90830000000000,2.49830000000000,3.27830000000000,3.85670000000000,2.58500000000000,3.17830000000000,3.64500000000000,2.49170000000000,2.50500000000000,3.86830000000000,3.54170000000000,2.87780000000000;90,2.31500000000000,3.71000000000000,1.50000000000000,2.72500000000000,2.13830000000000,3.19330000000000,1.80670000000000,1.75500000000000,2.08170000000000,2.85170000000000,2.59000000000000,1.66330000000000,2.10000000000000,1.72830000000000,2.18830000000000,2.46500000000000,1.83000000000000,1.84500000000000,2.30670000000000,2.23500000000000,1.85500000000000,3.06670000000000,3.46670000000000,2.47830000000000,2.91670000000000,3.26000000000000,1.81500000000000,1.95670000000000,3.34500000000000,3.15500000000000,2.41140000000000;120,1.75500000000000,3.01330000000000,1.22170000000000,1.68330000000000,1.71670000000000,2.59170000000000,1.69330000000000,1.35000000000000,1.67330000000000,2.60330000000000,2.20500000000000,1.43670000000000,1.86670000000000,1.53830000000000,1.88670000000000,1.73000000000000,1.34170000000000,1.46670000000000,2.18330000000000,1.49000000000000,1.76330000000000,2.46330000000000,2.69170000000000,1.81830000000000,2.42330000000000,2.86830000000000,1.54330000000000,1.76500000000000,2.96500000000000,2.61500000000000,1.97880000000000;150,1.68830000000000,2.33170000000000,1.16000000000000,1.62170000000000,1.60330000000000,2.34500000000000,1.59500000000000,1.23830000000000,1.62830000000000,2.08330000000000,1.91500000000000,1.23000000000000,1.62330000000000,1.11330000000000,1.89170000000000,1.41000000000000,1.20500000000000,1.37500000000000,2.04670000000000,1.24500000000000,1.50000000000000,2.06670000000000,2.23830000000000,1.66170000000000,2.08170000000000,2.45330000000000,1.27330000000000,1.64500000000000,2.71830000000000,2.30670000000000,1.74320000000000;180,1.48500000000000,1.59330000000000,1.12670000000000,1.29670000000000,1.31830000000000,2.02170000000000,1.50170000000000,1.07830000000000,1.49670000000000,1.98500000000000,1.68670000000000,1.31670000000000,1.69170000000000,1.14330000000000,1.55670000000000,1.46830000000000,1.16000000000000,1.24330000000000,1.53500000000000,1.10170000000000,1.22830000000000,1.97170000000000,1.76830000000000,1.39670000000000,1.84170000000000,2.24500000000000,1.11500000000000,1.34830000000000,2.43170000000000,2.04170000000000,1.53980000000000;210,1.24830000000000,1.86330000000000,1.05830000000000,1.08000000000000,1.37170000000000,1.68500000000000,1.26000000000000,1.02500000000000,1.46170000000000,1.79500000000000,1.62330000000000,1.09670000000000,1.62830000000000,1.12170000000000,1.46830000000000,1.74170000000000,1.07330000000000,1.07830000000000,1.34170000000000,1.06330000000000,1.15330000000000,2.17670000000000,1.43500000000000,1.39500000000000,1.59500000000000,2.14830000000000,1.04830000000000,1.20170000000000,1.86500000000000,1.60170000000000,1.42350000000000;240,1.16330000000000,1.15330000000000,1.07000000000000,1.27500000000000,1.28000000000000,1.36500000000000,1.35670000000000,1,1.66830000000000,1.56500000000000,1.56170000000000,1.03830000000000,1.64330000000000,1.11330000000000,1.38330000000000,1.30670000000000,1.08000000000000,1.19670000000000,1.47670000000000,1.04170000000000,1.14500000000000,1.56670000000000,1.49500000000000,1.36000000000000,1.70500000000000,2.18500000000000,1.02670000000000,1.04500000000000,1.62170000000000,1.52500000000000,1.34710000000000;300,1,1.02000000000000,1.00830000000000,1.05170000000000,1.19000000000000,1.63830000000000,1.01000000000000,1,1.51000000000000,1.35170000000000,1.40170000000000,1.04170000000000,1.42000000000000,1.01170000000000,1,1.41000000000000,1.14670000000000,1.03500000000000,1.33670000000000,1.03830000000000,1.05500000000000,1.40000000000000,1.61000000000000,1.43670000000000,1.40670000000000,2.23670000000000,1,1.09500000000000,1.37170000000000,1.63670000000000,1.26230000000000;360,1.04330000000000,1.04500000000000,1.03670000000000,1.08170000000000,1.05500000000000,1.18670000000000,1.05000000000000,1,1.57500000000000,1.53000000000000,1.18170000000000,1.05170000000000,1.41170000000000,1.01330000000000,1.12670000000000,1.00330000000000,1.08330000000000,1.00330000000000,1.17000000000000,1.00670000000000,1,1.37000000000000,1.23670000000000,1.23330000000000,1.32830000000000,2.23330000000000,1,1.02500000000000,1.35500000000000,1.23830000000000,1.18920000000000;420,1.11170000000000,1.22170000000000,1.00500000000000,1.08500000000000,1.02000000000000,1.16500000000000,1,1,1.46670000000000,1.32170000000000,1.36330000000000,1.12670000000000,1.28670000000000,1.00330000000000,1.11000000000000,1.00500000000000,1.13500000000000,1.01500000000000,1.10170000000000,1.00830000000000,1,1.28330000000000,1.24170000000000,1.31500000000000,1.13670000000000,2.20000000000000,1,1,1.61170000000000,1.32170000000000,1.18870000000000];
      
elseif round(max(eye_velo_hor_variable_degs_table)) == 54.0
    
    emp_data_cat3_standardization = [30,3.81670000000000,4.79000000000000,3.24500000000000,4.86000000000000,4.56170000000000,4.46830000000000,4.73500000000000,4.33170000000000,3.60330000000000,3.76170000000000,4.19830000000000,2.60330000000000,2.61670000000000,4.43170000000000,4.20670000000000,4.76830000000000,3.68670000000000,4.12330000000000,3.24330000000000,4.18500000000000,3.32330000000000,4.71670000000000,4.24170000000000,4.65830000000000,4.34330000000000,4.54830000000000,4.20670000000000,3.96170000000000,4.80500000000000,4.58170000000000,4.12080000000000;60,3.68500000000000,4.45330000000000,2.18330000000000,4.57830000000000,3.16000000000000,3.79670000000000,4.06670000000000,2.91170000000000,2.70170000000000,3.37330000000000,3.53330000000000,2.01000000000000,2.08500000000000,3.95830000000000,3.27670000000000,3.73000000000000,2.87670000000000,3.11500000000000,2.91500000000000,2.90170000000000,2.65670000000000,3.46330000000000,3.84330000000000,3.18000000000000,3.40830000000000,4.32830000000000,2.74330000000000,2.86000000000000,4.03670000000000,3.66330000000000,3.31650000000000;90,3.04500000000000,4.08670000000000,1.71330000000000,3.76170000000000,2.58330000000000,3.30170000000000,2.96670000000000,2.33170000000000,1.90000000000000,3.06830000000000,2.93670000000000,1.83000000000000,2.25170000000000,2.57000000000000,3.06830000000000,3.16500000000000,2.63500000000000,1.98670000000000,2.35830000000000,2.31830000000000,2.18330000000000,2.99000000000000,3.52670000000000,2.54000000000000,2.86170000000000,3.39000000000000,2.14330000000000,2.22330000000000,3.53330000000000,3.20330000000000,2.74910000000000;120,2.56170000000000,3.36830000000000,1.57500000000000,2.39000000000000,2.02670000000000,3.22170000000000,2.20330000000000,1.71170000000000,1.61830000000000,2.93000000000000,2.47170000000000,1.52500000000000,2.04830000000000,1.80170000000000,2.34670000000000,3.03000000000000,2.07670000000000,1.77670000000000,2.13670000000000,1.78500000000000,1.66830000000000,2.69170000000000,3.29170000000000,2.10670000000000,2.62500000000000,3.17330000000000,1.63830000000000,1.99170000000000,3.25670000000000,2.85330000000000,2.33010000000000;150,1.72330000000000,2.86670000000000,1.27330000000000,2.07170000000000,1.73670000000000,2.14670000000000,1.79670000000000,1.52170000000000,1.75000000000000,2.90000000000000,2.33830000000000,1.53330000000000,1.68500000000000,1.36330000000000,2.33830000000000,2.33170000000000,1.66500000000000,1.53330000000000,1.72330000000000,1.56500000000000,1.38670000000000,2.20830000000000,2.44330000000000,1.76500000000000,2.29830000000000,2.60830000000000,1.61500000000000,1.77670000000000,2.68500000000000,2.91500000000000,1.98550000000000;180,1.79000000000000,2.41170000000000,1.36000000000000,2.13500000000000,1.59330000000000,2.29000000000000,1.84000000000000,1.36830000000000,1.69830000000000,2.43000000000000,2.41170000000000,1.30170000000000,1.66830000000000,1.25000000000000,1.86000000000000,2.31670000000000,2.41000000000000,1.33500000000000,1.74330000000000,1.26330000000000,1.16670000000000,2.08670000000000,1.84000000000000,1.30670000000000,2.25830000000000,2.47000000000000,1.28670000000000,1.53670000000000,2.31330000000000,2.21000000000000,1.83170000000000;210,1.67830000000000,2.11670000000000,1.18670000000000,1.67000000000000,1.38170000000000,2.14500000000000,1.40000000000000,1.11500000000000,1.63000000000000,2.41170000000000,1.77500000000000,1.18670000000000,1.52000000000000,1.10500000000000,1.68330000000000,1.54500000000000,1.76670000000000,1.13670000000000,1.31670000000000,1.17670000000000,1.01000000000000,2.11670000000000,2.08330000000000,1.16830000000000,1.92170000000000,2.32670000000000,1.24000000000000,1.44670000000000,2.16670000000000,1.54670000000000,1.59910000000000;240,1.66500000000000,1.45830000000000,1.19330000000000,1.26330000000000,1.56830000000000,1.76170000000000,1.24670000000000,1.10500000000000,1.56170000000000,1.65170000000000,1.54500000000000,1.05670000000000,1.57330000000000,1.34500000000000,1.71330000000000,1.57000000000000,1.46330000000000,1.27000000000000,1.69670000000000,1.13830000000000,1.09670000000000,1.86170000000000,2.18500000000000,1.20830000000000,1.75830000000000,2.48500000000000,1.06830000000000,1.03500000000000,2.18170000000000,1.94500000000000,1.52240000000000;300,1.16170000000000,1.35170000000000,1.11500000000000,1.24670000000000,1.22670000000000,1.58000000000000,1.24830000000000,1.03500000000000,1.71670000000000,1.60170000000000,1.52170000000000,1.02500000000000,1.39330000000000,1.04670000000000,1.29330000000000,1.37830000000000,1.24830000000000,1.17830000000000,1.29670000000000,1.05830000000000,1,1.47830000000000,1.54170000000000,1.29500000000000,1.46000000000000,2.30000000000000,1.02330000000000,1,1.37670000000000,1.33330000000000,1.31770000000000;360,1.40170000000000,1.33670000000000,1.12170000000000,1.38170000000000,1.06500000000000,1.46000000000000,1.16000000000000,1,1.67000000000000,1.50830000000000,1.31000000000000,1.05830000000000,1.68170000000000,1.01330000000000,1.08830000000000,1.34830000000000,1.19000000000000,1.17500000000000,1.27330000000000,1.02000000000000,1,1.60830000000000,1.25830000000000,1.39830000000000,1.49330000000000,2.35830000000000,1.01830000000000,1.06670000000000,1.49330000000000,1.48670000000000,1.31480000000000;420,1.07170000000000,1.17000000000000,1.04830000000000,1.05500000000000,1,1.28830000000000,1.00500000000000,1,1.43000000000000,1.53330000000000,1.31830000000000,1.01830000000000,1.56500000000000,1.16330000000000,1.19170000000000,1.21170000000000,1.26670000000000,1.03830000000000,1.25500000000000,1.00330000000000,1.02330000000000,1.42830000000000,1.21170000000000,1.21830000000000,1.21330000000000,2.24330000000000,1.00500000000000,1,1.44170000000000,1.19670000000000,1.22050000000000];
    
else
    
    error('ATTENTION - there is no empirical cbu data the theoretical model cbu data can be compared with for the entered value of horizontal eye movement velocity');
    
end

% standardization process of empirical cbu data (see confluence)

emp_data_cat3_mean = sum(emp_data_cat3_standardization(1:end,end))/size(emp_data_cat3_standardization(1:end,end),1);
emp_data_cat3_sd = std(emp_data_cat3_standardization(1:end,end)); 
emp_data_cat3_mean_standardized = (emp_data_cat3_standardization(1:end,end)-emp_data_cat3_mean)/emp_data_cat3_sd; % standardized empirical mean values

emp_data_cat3_sub_standardized = (emp_data_cat3_standardization(1:end,2:end-1)-emp_data_cat3_mean)/emp_data_cat3_sd;
emp_data_cat3_sd_standardized = std(emp_data_cat3_sub_standardized')'; % standardized standard deviation of the empirical data

emp_data_cat3_table_standardized = [emp_data_cat3_standardization(1:end,1),emp_data_cat3_mean_standardized, emp_data_cat3_sd_standardized]; % final standardized matrix for empirical cbu data (first column = frame rates; second column = mean values; third column = standard deviation)

% standardization process of theoretical cbu model data refering to option 1 (see confluence)

cbu_score_opt1_deg2_mean = sum(loop_table(1:end,7))/size(loop_table(1:end,7),1);
cbu_score_opt1_deg2_sd = std(loop_table(1:end,7));
cbu_score_opt1_deg2_mean_standardized = (loop_table(1:end,7)-cbu_score_opt1_deg2_mean)/cbu_score_opt1_deg2_sd; % standardized theoretical mean values

cbu_score_opt1_deg2_variance_standardized = max(cbu_score_opt1_deg2_mean_standardized).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4); % definition of the standardized variance for the theoretical cbu score

cbu_score_opt1_deg2_table_standardized = [loop_table(1:end,4),cbu_score_opt1_deg2_mean_standardized, cbu_score_opt1_deg2_variance_standardized]; % final standardized matrix for theoretical cbu data (first column = frame rates; second column = mean values; third column = variance)

% standardization process of theoretical cbu model data refering to option 2 (see confluence)

cbu_score_opt2_deg2_mean = sum(loop_table(1:end,10))/size(loop_table(1:end,10),1);
cbu_score_opt2_deg2_sd = std(loop_table(1:end,10));
cbu_score_opt2_deg2_mean_standardized = (loop_table(1:end,10)-cbu_score_opt2_deg2_mean)/cbu_score_opt2_deg2_sd; % standardized theoretical mean values

cbu_score_opt2_deg2_variance_standardized = max(cbu_score_opt2_deg2_mean_standardized).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4); % definition of the standardized variance for the theoretical cbu score

cbu_score_opt2_deg2_table_standardized = [loop_table(1:end,4),cbu_score_opt2_deg2_mean_standardized, cbu_score_opt2_deg2_variance_standardized]; % final standardized matrix for theoretical cbu data (first column = frame rates; second column = mean values; third column = variance)

% standardization process of theoretical cbu model data refering to option 2 after applying stevens' psychophysical power function filter

stevens1 = (1.0*(loop_table(1:end,10)).^0.3); % applying stevens' psychophysical power function filter (variable beta values: 0.3, 0.5, 0.7 and 1.0) to model output 
stevens2 = (1.0*(loop_table(1:end,10)).^0.5);
stevens3 = (1.0*(loop_table(1:end,10)).^0.7);
stevens4 = (1.0*(loop_table(1:end,10)).^1.0);

cbu_score_opt2_deg2_mean1 = sum(stevens1)/size(stevens1,1);
cbu_score_opt2_deg2_sd1 = std(stevens1);
cbu_score_opt2_deg2_mean_standardized1 = (stevens1-cbu_score_opt2_deg2_mean1)/cbu_score_opt2_deg2_sd1; % standardized theoretical mean values
cbu_score_opt2_deg2_variance_standardized1 = max(cbu_score_opt2_deg2_mean_standardized1).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4); % definition of the standardized variance for the theoretical cbu score
cbu_score_opt2_deg2_table_standardized1 = [loop_table(1:end,4),cbu_score_opt2_deg2_mean_standardized1, cbu_score_opt2_deg2_variance_standardized1]; % final standardized matrix for theoretical cbu data (first column = frame rates; second column = mean values; third column = variance)

cbu_score_opt2_deg2_mean2 = sum(stevens2)/size(stevens2,1);
cbu_score_opt2_deg2_sd2 = std(stevens2);
cbu_score_opt2_deg2_mean_standardized2 = (stevens2-cbu_score_opt2_deg2_mean2)/cbu_score_opt2_deg2_sd2; % standardized theoretical mean values
cbu_score_opt2_deg2_variance_standardized2 = max(cbu_score_opt2_deg2_mean_standardized2).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4); % definition of the standardized variance for the theoretical cbu score
cbu_score_opt2_deg2_table_standardized2 = [loop_table(1:end,4),cbu_score_opt2_deg2_mean_standardized2, cbu_score_opt2_deg2_variance_standardized2]; % final standardized matrix for theoretical cbu data (first column = frame rates; second column = mean values; third column = variance)

cbu_score_opt2_deg2_mean3 = sum(stevens3)/size(stevens3,1);
cbu_score_opt2_deg2_sd3 = std(stevens3);
cbu_score_opt2_deg2_mean_standardized3 = (stevens3-cbu_score_opt2_deg2_mean3)/cbu_score_opt2_deg2_sd3; % standardized theoretical mean values
cbu_score_opt2_deg2_variance_standardized3 = max(cbu_score_opt2_deg2_mean_standardized3).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4); % definition of the standardized variance for the theoretical cbu score
cbu_score_opt2_deg2_table_standardized3 = [loop_table(1:end,4),cbu_score_opt2_deg2_mean_standardized3, cbu_score_opt2_deg2_variance_standardized3]; % final standardized matrix for theoretical cbu data (first column = frame rates; second column = mean values; third column = variance)

cbu_score_opt2_deg2_mean4 = sum(stevens4)/size(stevens4,1);
cbu_score_opt2_deg2_sd4 = std(stevens4);
cbu_score_opt2_deg2_mean_standardized4 = (stevens4-cbu_score_opt2_deg2_mean4)/cbu_score_opt2_deg2_sd4; % standardized theoretical mean values
cbu_score_opt2_deg2_variance_standardized4 = max(cbu_score_opt2_deg2_mean_standardized4).*0.45.*min(loop_table(1:end,4))./loop_table(1:end,4); % definition of the standardized variance for the theoretical cbu score
cbu_score_opt2_deg2_table_standardized4 = [loop_table(1:end,4),cbu_score_opt2_deg2_mean_standardized4, cbu_score_opt2_deg2_variance_standardized4]; % final standardized matrix for theoretical cbu data (first column = frame rates; second column = mean values; third column = variance)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% EXPONENTIAL FUNCTION >> generates exponential function for curve fitting regarding empirical cbu data values (after standardization of the empirical data, see above)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
x = emp_data_cat3_table_standardized(1:end,1);
y = emp_data_cat3_table_standardized(1:end,2)+abs(min(emp_data_cat3_table_standardized(1:end,2)));
    
genexp = @(b,x) b(1).*exp(b(2).*(x(:,1)))+b(3);
b0 = [0, 0, 0];
framepred = linspace(0, 840, 1000)';
    
mod1 = fitnlm(x, y, genexp, b0);
cbupred = predict(mod1,framepred)-abs(min(emp_data_cat3_table_standardized(1:end,2)));
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (3) >> comparison of empirical cbu data with theoretical model cbu data (after standardization of the theoretical cbu model data and empirical cbu data)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
f5 = figure (5); hold on; set(gca,'FontSize',15); % figure 5 refers to cbu score calculation with option 1

p(1) = plot(cbu_score_opt1_deg2_table_standardized(1:end,1),cbu_score_opt1_deg2_table_standardized(1:end,2),'--','Color',[0.9 0.0 0.0],'LineWidth',1.0); % plot theoretical cbu score data (only mean value)
% p(1) = boundedline(cbu_score_opt1_deg2_table_standardized(1:end,1),cbu_score_opt1_deg2_table_standardized(1:end,2),cbu_score_opt1_deg2_table_standardized(1:end,3),'r','alpha'); % plot theoretical cbu score data (mean and spread values) 

p(2) = plot(framepred,cbupred,'b--','LineWidth',1.0); % plot empirical cbu score data (curve fitted exponential function)
p(3) = plot(emp_data_cat3_table_standardized(1:end,1),emp_data_cat3_table_standardized(1:end,2),'o','MarkerFaceColor','w','MarkerEdgeColor','b','LineWidth',1.0,'MarkerSize',4.0); % plot empirical cbu score data (mean values)
p(4) = boundedline(emp_data_cat3_table_standardized(1:end,1), emp_data_cat3_table_standardized(1:end,2), emp_data_cat3_table_standardized(1:end,3),':','alpha'); % plot empirical cbu score data (standard deviation)

p(5) = line([frame_rate_hz_threshold_ph1_ph2 frame_rate_hz_threshold_ph1_ph2],[[min([(cbu_score_opt1_deg2_table_standardized(1:end,2)-cbu_score_opt1_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]) max([(cbu_score_opt1_deg2_table_standardized(1:end,2)+cbu_score_opt1_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])]],'Color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth', 1.0); % plot phase threshold (phase 1 to phase 2)
p(6) = line([frame_rate_hz_threshold_ph2_ph3 frame_rate_hz_threshold_ph2_ph3],[[min([(cbu_score_opt1_deg2_table_standardized(1:end,2)-cbu_score_opt1_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]) max([(cbu_score_opt1_deg2_table_standardized(1:end,2)+cbu_score_opt1_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])]],'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth', 1.0); % plot phase threshold (phase 2 to phase 3)

xlim([frame_rate_hz_min frame_rate_hz_max]);
set(gca,'xTick',frame_rate_hz_min:(frame_rate_hz_max-frame_rate_hz_min)/13:frame_rate_hz_max);
xlabel('FRAME RATE [HZ]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))])-(max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10 max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])+(max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10]);
set(gca,'yTick',min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))])-(max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10:((max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])+(max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10)-(min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))])-(max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10))/10:max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])+(max([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt1_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10);
ylabel('STANDARDIZED CBU [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('EMPIRICAL VS. THEORETICAL CBU (OPTION 1, W/O STEVENS)', 'Fontsize', 21);
legend('theoretical CBU index (model, option 1, w/o stevens)','empirical CBU score (fitted exponential curve)', 'empirical CBU score (MN)', 'empirical CBU score (SD)', 'empirical CBU score (linear interpolation)','CBU phase transition threshold (PH1/PH2)','CBU phase transition threshold (PH2/PH3)','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f6 = figure (6); hold on; set(gca,'FontSize',15); % figure 6 refers to cbu score calculation with option 2

p(1) = plot(cbu_score_opt2_deg2_table_standardized(1:end,1),cbu_score_opt2_deg2_table_standardized(1:end,2),'--','Color',[0.9 0.0 0.0],'LineWidth',1.0); % plot theoretical cbu score data (only mean value) 
% p(1) = boundedline(cbu_score_opt2_deg2_table_standardized(1:end,1),cbu_score_opt2_deg2_table_standardized(1:end,2),cbu_score_opt2_deg2_table_standardized(1:end,3),'g','alpha'); % plot theoretical cbu score data (mean and spread values) 

p(2) = plot(framepred,cbupred,'b--','LineWidth',1.0); % plot empirical cbu score data (curve fitted exponential function)
p(3) = plot(emp_data_cat3_table_standardized(1:end,1),emp_data_cat3_table_standardized(1:end,2),'o','MarkerFaceColor','w','MarkerEdgeColor','b','LineWidth',1.0,'MarkerSize',4.0); % plot empirical cbu score data (mean values)
p(4) = boundedline(emp_data_cat3_table_standardized(1:end,1), emp_data_cat3_table_standardized(1:end,2), emp_data_cat3_table_standardized(1:end,3),':','alpha'); % plot empirical cbu score data (standard deviation)

p(5) = line([frame_rate_hz_threshold_ph1_ph2 frame_rate_hz_threshold_ph1_ph2],[min([(cbu_score_opt2_deg2_table_standardized(1:end,2)-cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]) max([(cbu_score_opt2_deg2_table_standardized(1:end,2)+cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])],'Color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth', 1.0); % plot phase threshold (phase 1 to phase 2)
p(6) = line([frame_rate_hz_threshold_ph2_ph3 frame_rate_hz_threshold_ph2_ph3],[min([(cbu_score_opt2_deg2_table_standardized(1:end,2)-cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]) max([(cbu_score_opt2_deg2_table_standardized(1:end,2)+cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])],'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth', 1.0); % plot phase threshold (phase 2 to phase 3)

xlim([frame_rate_hz_min frame_rate_hz_max]);
set(gca,'xTick',frame_rate_hz_min:(frame_rate_hz_max-frame_rate_hz_min)/13:frame_rate_hz_max);
xlabel('FRAME RATE [HZ]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))])-(max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10 max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])+(max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10]);
set(gca,'yTick',min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))])-(max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10:((max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])+(max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10)-(min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))])-(max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10))/10:max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])+(max([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])-min([cbu_score_opt2_deg2_table_standardized(1:end,2);(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]))/10);
ylabel('STANDARDIZED CBU [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('EMPIRICAL VS. THEORETICAL CBU (OPTION 2, W/O STEVENS)', 'Fontsize', 21);
legend('theoretical CBU index (model, option 2, w/o stevens)','empirical CBU score (fitted exponential curve)', 'empirical CBU score (MN)', 'empirical CBU score (SD)', 'empirical CBU score (linear interpolation)','CBU phase transition threshold (PH1/PH2)','CBU phase transition threshold (PH2/PH3)','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f7 = figure (7); hold on; set(gca,'FontSize',15); % figure 7 refers to cbu score calculation with option 2 (including stevens filter; applied beta: 0.3 to 0.7)

p(1) = plot(framepred,cbupred,'b--','LineWidth',1.0); % plot empirical cbu score data (curve fitted exponential function)
p(2) = plot(emp_data_cat3_table_standardized(1:end,1),emp_data_cat3_table_standardized(1:end,2),'o','MarkerFaceColor','w','MarkerEdgeColor','b','LineWidth',1.0,'MarkerSize',4.0); % plot empirical cbu score data (mean values)
p(3) = boundedline(emp_data_cat3_table_standardized(1:end,1), emp_data_cat3_table_standardized(1:end,2), emp_data_cat3_table_standardized(1:end,3),':','alpha'); % plot empirical cbu score data (standard deviation)

p(4) = plot(cbu_score_opt2_deg2_table_standardized4(1:end,1),cbu_score_opt2_deg2_table_standardized4(1:end,2),'--','Color',[0.5 0.5 0.5],'LineWidth',0.75); % plot theoretical cbu score data (only mean value)

p(5) = plot(cbu_score_opt2_deg2_table_standardized2(1:end,1),cbu_score_opt2_deg2_table_standardized2(1:end,2),'--','Color',[0.9 0.0 0.0],'LineWidth',1.0);
p(6) = boundedline(cbu_score_opt2_deg2_table_standardized2(1:end,1),cbu_score_opt2_deg2_table_standardized2(1:end,2),cbu_score_opt2_deg2_table_standardized3(1:end,2)-cbu_score_opt2_deg2_table_standardized2(1:end,2), 'r:','alpha'); % plot theoretical cbu score data (mean and spread values)

%p(7) = line([frame_rate_hz_threshold_ph1_ph2 frame_rate_hz_threshold_ph1_ph2],[min([(cbu_score_opt2_deg2_table_standardized(1:end,2)-cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]) max([(cbu_score_opt2_deg2_table_standardized(1:end,2)+cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])],'Color',[0.5,0.5,0.5],'LineStyle','-.','LineWidth', 1.0); % plot phase threshold (phase 1 to phase 2)
%p(8) = line([frame_rate_hz_threshold_ph2_ph3 frame_rate_hz_threshold_ph2_ph3],[min([(cbu_score_opt2_deg2_table_standardized(1:end,2)-cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)-emp_data_cat3_table_standardized(1:end,3))]) max([(cbu_score_opt2_deg2_table_standardized(1:end,2)+cbu_score_opt2_deg2_table_standardized(1:end,3));(emp_data_cat3_table_standardized(1:end,2)+emp_data_cat3_table_standardized(1:end,3))])],'Color',[0.5,0.5,0.5],'LineStyle','--','LineWidth', 1.0); % plot phase threshold (phase 2 to phase 3)

xlim([frame_rate_hz_min frame_rate_hz_max]);
set(gca,'xTick',frame_rate_hz_min:(frame_rate_hz_max-frame_rate_hz_min)/13:frame_rate_hz_max);
xlabel('FRAME RATE [HZ]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([-1.684 4.875]);
set(gca,'yTick',-1.684:(4.875-(-1.684))/10:4.875);
ylabel('STANDARDIZED CBU [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('EMPIRICAL VS. THEORETICAL CBU (OPTION 2, STEVENS)', 'Fontsize', 21);
legend('empirical CBU score (fitted exponential curve)', 'empirical CBU score (MN)', 'empirical CBU score (SD)', 'empirical CBU score (linear interpolation)','theoretical CBU index (model, option 2, w/o stevens)','theoretical CBU index (model, option 2, stevens with beta=0.5)','theoretical CBU index (model, option 2, stevens with beta=0.3-0.7)','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CONDITION (4) >> for comparison between empirical and theoretical cbu scores  
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if frame_rate_hz_min ~= 30.0 || frame_rate_hz_max ~= 420.0 || (30/frame_rate_hz_step) ~= fix(30/frame_rate_hz_step)
    
    cprintf('err','ATTENTION - discrete frame rate values of choice not adequate for comparison between theoretical and empirical cbu perception (see VARIABLE), frame rate values of 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 300.0, 360.0, 420.0hz necessary \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (4) >> further comparison of empirical cbu data with theoretical model cbu data (after standardization of the theoretical cbu model data and empirical cbu data)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cbu_score_opt1_deg2_table_standardized_short = cbu_score_opt1_deg2_table_standardized(ismember(cbu_score_opt1_deg2_table_standardized(:,1),[emp_data_cat3_table_standardized(:,1)]),:);
cbu_score_opt2_deg2_table_standardized_short = cbu_score_opt2_deg2_table_standardized(ismember(cbu_score_opt2_deg2_table_standardized(:,1),[emp_data_cat3_table_standardized(:,1)]),:);
cbu_score_opt2_deg2_table_standardized1_short = cbu_score_opt2_deg2_table_standardized1(ismember(cbu_score_opt2_deg2_table_standardized1(:,1),[emp_data_cat3_table_standardized(:,1)]),:);

%{
f8 = figure (8); hold on; set(gca,'FontSize',15); % figure 8 (bland-altman-plot with sample values) refers to cbu score calculation with option 1

BlandAltmanPlot(emp_data_cat3_table_standardized(:,2),cbu_score_opt1_deg2_table_standardized_short(:,2),'plotCI',[false])

xlabel('CBU - MEAN  [ ( T + E ) / 2 ] ','Fontsize',21);
xtickformat('%,.1f');

ylabel('CBU - DIFFERENCE  [ T - E ]','Fontsize',21);
ytickformat('%,.1f');

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('BLAND ALTMAN PLOT (SAMPLE, OPTION 1, STANDARDIZED)', 'Fontsize', 21);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f9 = figure (9); hold on; set(gca,'FontSize',15); % figure 9 (bland-altman-plot with sample values) refers to cbu score calculation with option 2

BlandAltmanPlot(emp_data_cat3_table_standardized(:,2),cbu_score_opt2_deg2_table_standardized_short(:,2),'plotCI',[false])

xlabel('CBU - MEAN  [ ( T + E ) / 2 ] ','Fontsize',21);
xtickformat('%,.1f');

ylabel('CBU - DIFFERENCE  [ T - E ]','Fontsize',21);
ytickformat('%,.1f');

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('BLAND ALTMAN PLOT (SAMPLE, OPTION 2, STANDARDIZED)', 'Fontsize', 21);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f10 = figure (10); hold on; set(gca,'FontSize',15); % figure 10 (bland-altman-plot with single subject values) refers to cbu score calculation with option 1

BlandAltmanPlot2(emp_data_cat3_sub_standardized(:),repmat(cbu_score_opt1_deg2_table_standardized_short(:,2),30,1),'plotCI',[false])

xlabel('CBU - MEAN  [ ( T + E ) / 2 ] ','Fontsize',22);
xtickformat('%,.1f');

ylabel('CBU - DIFFERENCE  [ T - E ]','Fontsize',22);
ytickformat('%,.1f');

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('BLAND ALTMAN PLOT (SUBJECT, OPTION 1, STANDARDIZED)', 'Fontsize', 21);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f11 = figure (11); hold on; set(gca,'FontSize',15); % figure 11 (bland-altman-plot with single subject values) refers to cbu score calculation with option 2

BlandAltmanPlot2(emp_data_cat3_sub_standardized(:),repmat(cbu_score_opt2_deg2_table_standardized_short(:,2),30,1),'plotCI',[false])

xlabel('CBU - MEAN  [ ( T + E ) / 2 ] ','Fontsize',21);
xtickformat('%,.1f');

ylabel('CBU - DIFFERENCE  [ T - E ]','Fontsize',21);
ytickformat('%,.1f');

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('BLAND ALTMAN PLOT (SUBJECT, OPTION 2, STANDARDIZED)', 'Fontsize', 21);

box off;
set(gcf,'WindowStyle','docked');
hold off;
%}

f12 = figure (12); hold on; set(gca,'FontSize',15); % figure 12 (difference scatter plot with linear regression) refers to cbu score calculation with option 1

p(1) = plot(cbu_score_opt1_deg2_table_standardized_short(:,2),emp_data_cat3_table_standardized(:,2),'ro','LineWidth',0.5, 'MarkerFaceColor',[1.0 0.85 0.85], 'MarkerSize',5.0);

% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x = cbu_score_opt1_deg2_table_standardized_short(:,2);
y = emp_data_cat3_table_standardized(:,2); 
    
genlin = @(b,x) b(1).*x(:,1)+b(2); % linear function for correlation between luminance level and cbu score
b0 = [0, 0];
xpred = linspace((min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10), (max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10), 1000)';
mod1 = fitnlm(x, y, genlin, b0);
ypred = predict(mod1,xpred);

p(2) = plot(xpred,ypred,'r--','LineWidth',1.0);
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

p(3) = plot([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10],[min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10],'k--','LineWidth',1.0);
p(4) = plot([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10],[0 0],'k-.','LineWidth',0.5);
p(5) = plot([0 0],[min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10],'k-.','LineWidth',0.5);

% axis equal;

xlim([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10]);
set(gca,'xTick',min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10:((max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10)-(min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10))/10:max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10);
xlabel('MODEL CBU INDEX [-]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));

ylim([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10]);
set(gca,'yTick',min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10:((max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10)-(min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10))/10:max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt1_deg2_table_standardized_short(:,2)]))/10);
ylabel('STUDY CBU SCORE [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('CBU COMPARISON (STANDARDIZED, SAMPLE, OPTION 1, W/O STEVENS)', 'Fontsize', 21);
legend('CBU data point','linear CBU regression','bisecting line', 'zero line','Location','northwest','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f13 = figure (13); hold on; set(gca,'FontSize',15); % figure 13 (difference scatter plot with linear regression) refers to cbu score calculation with option 2

p(1) = plot(cbu_score_opt2_deg2_table_standardized_short(:,2),emp_data_cat3_table_standardized(:,2),'ro','LineWidth',0.5, 'MarkerFaceColor',[1.0 0.85 0.85], 'MarkerSize',5.0);

% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x = cbu_score_opt2_deg2_table_standardized_short(:,2); 
y = emp_data_cat3_table_standardized(:,2); 
   
genlin = @(b,x) b(1).*x(:,1)+b(2); % linear function for correlation between luminance level and cbu score
b0 = [0, 0];
xpred = linspace((min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10), (max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10), 1000)';
mod1 = fitnlm(x, y, genlin, b0);
ypred = predict(mod1,xpred);

p(2) = plot(xpred,ypred,'r--','LineWidth',1.0);
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

p(3) = plot([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10],[min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10],'k--','LineWidth',1.0);
p(4) = plot([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10],[0 0],'k-.','LineWidth',0.5);
p(5) = plot([0 0],[min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10],'k-.','LineWidth',0.5);

% axis equal;

xlim([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10]);
set(gca,'xTick',min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10:((max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10)-(min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10))/10:max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10);
xlabel('MODEL CBU INDEX [-]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));

ylim([min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10 max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10]);
set(gca,'yTick',min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10:((max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10)-(min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10))/10:max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])+(max([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)])-min([emp_data_cat3_table_standardized(:,2);cbu_score_opt2_deg2_table_standardized_short(:,2)]))/10);
ylabel('STUDY CBU SCORE [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('CBU COMPARISON (STANDARDIZED, SAMPLE, OPTION 2, W/O STEVENS)', 'Fontsize', 21);
legend('CBU data point','linear CBU regression','bisecting line', 'zero line','Location','northwest','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f14 = figure (14); hold on; set(gca,'FontSize',15); % figure 14 (difference scatter plot with linear regression) refers to cbu score calculation with option 2 after applying stevens filter (beta=0.3)

p(1) = plot(cbu_score_opt2_deg2_table_standardized1_short(:,2),emp_data_cat3_table_standardized(:,2),'ro','LineWidth',0.5, 'MarkerFaceColor',[1.0 0.85 0.85], 'MarkerSize',5.0);

% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
x = cbu_score_opt2_deg2_table_standardized1_short(:,2); 
y = emp_data_cat3_table_standardized(:,2); 
   
genlin = @(b,x) b(1).*x(:,1)+b(2); % linear function for correlation between luminance level and cbu score
b0 = [0, 0];
xpred = linspace(-1.35, 4.85, 1000)';
mod1 = fitnlm(x, y, genlin, b0);
ypred = predict(mod1,xpred);

p(2) = plot(xpred,ypred,'r--','LineWidth',1.0);
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

p(3) = plot([-1.35 4.85],[-1.35 4.85],'k--','LineWidth',1.0);
p(4) = plot([-1.35 4.85],[0 0],'k-.','LineWidth',0.5);
p(5) = plot([0 0],[-1.35 4.85],'k-.','LineWidth',0.5);

% axis equal;

xlim([-1.35 4.85]);
set(gca,'xTick',-1.35:(4.85-(-1.35))/10:4.85);
xlabel('MODEL CBU INDEX [-]','Fontsize',21);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));

ylim([-1.35 4.85]);
set(gca,'yTick',-1.35:(4.85-(-1.35))/10:4.85);
ylabel('STUDY CBU SCORE [-]','Fontsize',21);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('CBU COMPARISON (STANDARDIZED, SAMPLE, OPTION 2, STEVENS)', 'Fontsize', 21);
legend('CBU data point','linear CBU regression','bisecting line', 'zero line','Location','northwest','Fontsize',9);

box off;
set(gcf,'WindowStyle','docked');
hold off;

f15 = figure (15); hold on; set(gca,'FontSize',15); set(gcf,'WindowStyle','docked'); % figure 15 (difference boxplot) refers to cbu score calculation with option 1

p(1) = plot([0 450],[0 0],'r--','LineWidth',1.0);

boxplot(repmat(cbu_score_opt1_deg2_table_standardized_short(:,2)',size(emp_data_cat3_sub_standardized',1),1)-emp_data_cat3_sub_standardized','Labels',{'30.0','60.0','90.0','120.0','150.0','180.0','210.0','240.0','300.0','360.0','420'},'Colors','k','Symbol','k+');
lines1 = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines1, 'Color', 'b');

xlabel('FRAME RATE [HZ]','Fontsize',21);

aa = cbu_score_opt1_deg2_table_standardized_short(:,2)-emp_data_cat3_sub_standardized;
bb = min(aa(:));
cc = max(aa(:));

ylim([bb-(cc-bb)/10 cc+(cc-bb)/10]);
set(gca,'yTick',bb-(cc-bb)/10:((cc+(cc-bb)/10)-(bb-(cc-bb)/10))/10:cc+(cc-bb)/10);
ylabel('CBU - DIFFERENCE [M-S]','Fontsize',21);
ytickformat('%,.2f');

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('CBU COMPARISON (STANDARDIZED, SAMPLE, OPTION 1, W/O STEVENS)', 'Fontsize', 21);
legend('zero line','Location','northeast','Fontsize',9);

box off;
hold off;

f16 = figure (16); hold on; set(gca,'FontSize',15); set(gcf,'WindowStyle','docked');% figure 16 (difference boxplot) refers to cbu score calculation with option 2

p(1) = plot([0 450],[0 0],'r--','LineWidth',1.0);

boxplot(repmat(cbu_score_opt2_deg2_table_standardized_short(:,2)',size(emp_data_cat3_sub_standardized',1),1)-emp_data_cat3_sub_standardized','Labels',{'30.0','60.0','90.0','120.0','150.0','180.0','210.0','240.0','300.0','360.0','420'},'Colors','k','Symbol','k+');
lines2 = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines2, 'Color', 'b');

xlabel('FRAME RATE [HZ]','Fontsize',21);

dd = cbu_score_opt2_deg2_table_standardized_short(:,2)-emp_data_cat3_sub_standardized;
ee = min(dd(:));
ff = max(dd(:));

ylim([ee-(ff-ee)/10 ff+(ff-ee)/10]);
set(gca,'yTick',ee-(ff-ee)/10:((ff+(ff-ee)/10)-(ee-(ff-ee)/10))/10:ff+(ff-ee)/10);
ylabel('CBU - DIFFERENCE [M-S]','Fontsize',21);
ytickformat('%,.2f');

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('CBU COMPARISON (STANDARDIZED, SAMPLE, OPTION 2, W/O STEVENS)', 'Fontsize', 21);
legend('zero line','Location','northeast','Fontsize',9);

box off;
hold off;

f17 = figure (17); hold on; set(gca,'FontSize',15); set(gcf,'WindowStyle','docked');% figure 17 (difference boxplot) refers to cbu score calculation with option 2 after applying stevens filter (beta=0.3)

p(1) = plot([0 450],[0 0],'r--','LineWidth',1.0);

boxplot(repmat(cbu_score_opt2_deg2_table_standardized1_short(:,2)',size(emp_data_cat3_sub_standardized',1),1)-emp_data_cat3_sub_standardized','Labels',{'30.0','60.0','90.0','120.0','150.0','180.0','210.0','240.0','300.0','360.0','420'},'Colors','k','Symbol','k+');
lines2 = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines2, 'Color', 'b');

xlabel('FRAME RATE [HZ]','Fontsize',21);

dd = cbu_score_opt2_deg2_table_standardized1_short(:,2)-emp_data_cat3_sub_standardized;
ee = min(dd(:));
ff = max(dd(:));

ylim([-2.02 4.286]);
set(gca,'yTick',-2.02:(4.286-(-2.02))/10:4.286);
ylabel('CBU - DIFFERENCE [M-S]','Fontsize',21);
ytickformat('%,.2f');

grid on;
ax = gca;
ax.GridColor = [0.4, 0.4, 0.4];

title ('CBU COMPARISON (STANDARDIZED, SAMPLE, OPTION 2, STEVENS)', 'Fontsize', 21);
legend('zero line','Location','northeast','Fontsize',9);

box off;
hold off;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% SEGMENT B - VARIATION OF HORIZONTAL EYE MOVEMENT VELOCITY [DEG/SEC] - CODE IS UNFINISHED AND THEREFORE DEACTIVATED, BEFORE RUNNING THE CODE, THE CODE STRUCTURE OF SEGMENT A MUST BE TRANSFERED 
%{
% first step - determination of the constant parameter (frame rate, duty cycle, content width, subframe number) 
% second step - definition of the range (min./max.) for and the distance (step) between the discrete values of eye movement velocity as the variable parameter can be defined
% third step - run section exclusively, not the whole function (!)

% CONSTANTS >> values that are not changed during the loop
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
frame_rate_hz = 30.0;
duty_cycle = 0.30;
content_width_px = 40.0;
subframe_number = 3.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% VARIABLE >> value that is changed over the determined range below, eye_velo_hor_degs needs to be a positive decimal number (including zero)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
eye_velo_hor_degs_min = 10.0;
eye_velo_hor_degs_max = 300.0;
eye_velo_hor_degs_step = 10.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CONDITION (1) >> for calculation of theoretial cbu scores  
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ((eye_velo_hor_degs_max-eye_velo_hor_degs_min)/eye_velo_hor_degs_step) ~= fix((eye_velo_hor_degs_max-eye_velo_hor_degs_min)/eye_velo_hor_degs_step)
    
    cprintf('err','ATTENTION - lower and upper borders of eye velocity range do not match with eye velocity steps (see VARIABLE) \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CBU SCORE >> gives out the different cbu scores (option 1, option 2 and option 3) calculated by the theoretical cbu model 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for eye_velo_hor_degs = (eye_velo_hor_degs_min):(eye_velo_hor_degs_step):(eye_velo_hor_degs_max) % eye velocity range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % cbu score - option 1: area of retinal CBU stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
       
    cbu_score_opt1_table(aa,1) = eye_velo_hor_degs; % first column includes eye velocities      
    cbu_score_opt1_table(aa,2) = cbu_score_opt1; % second column includes theoretically calculated cbu score (mean value)

    % cbu score - option 2: area & intensity of retinal CBU stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
   
    cbu_score_opt2_table(aa,1) = eye_velo_hor_degs;    
    cbu_score_opt2_table(aa,2)= cbu_score_opt2; 

    % cbu score - option 3: area & intensity & color of retinal CBU stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description

    % cbu_score_opt3_table(aa,1) = ... ;   
    % cbu_score_opt3_table(aa,2)= ... ; 
    
end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% WHITE CONTENT >> gives out the score for retinal stimulation with unbiased white content (option 1, option 2 and option 3)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for eye_velo_hor_degs = (eye_velo_hor_degs_min):(eye_velo_hor_degs_step):(eye_velo_hor_degs_max) % eye velocity range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % white content - option 1: area of retinal white content stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
    
    white_content_opt1_table(aa,1) = eye_velo_hor_degs;        
    white_content_opt1_table(aa,2) = white_content_opt1;

    % white content - option 2: area & intensity of retinal white content stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
    
    white_content_opt2_table(aa,1) = eye_velo_hor_degs;        
    white_content_opt2_table(aa,2)= white_content_opt2;    

    % white content - option 3: area & intensity & color of retinal white content stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description
    
    % white_content_opt3_table(aa,1) = ... ;        
    % white_content_opt3_table(aa,2)= ... ;
    
end

cbu_score_opt1_table, white_content_opt1_table, reference_opt1

cbu_score_opt2_table, white_content_opt2_table, reference_opt2
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% THRESHOLD CALCULATION >> gives out the values for the variable parameter (see above) for which a phase transition from phase 1 to phase 2 and from phase 2 to phase 3 is accomplished
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
eye_velo_hor_degs_threshold12 = ((content_width_px.*subframe_number.*frame_rate_hz)/(2-duty_cycle)).*pixel_pitch_hor_deg
eye_velo_hor_degs_threshold23 = ((content_width_px.*subframe_number.*frame_rate_hz)/(1-duty_cycle)).*pixel_pitch_hor_deg
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (1) >> illustrates the dependency between the variable parameter (x-axis) and cbu score, white content and reference as target parameter (y-axis)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
close all; 

f1 = figure (1); hold on; % figure 1 refers to target paremeter calculation with option 1
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt1_table(1:end,1),cbu_score_opt1_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt1_table(1:end,1),white_content_opt1_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([eye_velo_hor_degs_min eye_velo_hor_degs_max],[reference_opt1 reference_opt1],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([eye_velo_hor_degs_threshold12 eye_velo_hor_degs_threshold12],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([eye_velo_hor_degs_threshold23 eye_velo_hor_degs_threshold23],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([eye_velo_hor_degs_min eye_velo_hor_degs_max]);
set(gca,'xTick',eye_velo_hor_degs_min:(eye_velo_hor_degs_max-eye_velo_hor_degs_min)/10:eye_velo_hor_degs_max);
xlabel('HORIZONTAL EYE MOVEMENT VELOCITY [DEG/SEC]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])]);
set(gca,'yTick',0:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])/10:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1]));
ylabel('COLOR BREAK-UP SCORE [MM2]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF EYE MOVEMENT VELOCITY', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;

f2 = figure (2); hold on; % figure 2 refers to target parameter calculation with option 2
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt2_table(1:end,1),cbu_score_opt2_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt2_table(1:end,1),white_content_opt2_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([eye_velo_hor_degs_min eye_velo_hor_degs_max],[reference_opt2 reference_opt2],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([eye_velo_hor_degs_threshold12 eye_velo_hor_degs_threshold12],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([eye_velo_hor_degs_threshold23 eye_velo_hor_degs_threshold23],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([eye_velo_hor_degs_min eye_velo_hor_degs_max]);
set(gca,'xTick',eye_velo_hor_degs_min:(eye_velo_hor_degs_max-eye_velo_hor_degs_min)/10:eye_velo_hor_degs_max);
xlabel('HORIZONTAL EYE MOVEMENT VELOCITY [DEG/SEC]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])]);
set(gca,'yTick',0:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])/10:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2]));
ylabel('COLOR BREAK-UP SCORE [---]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF EYE MOVEMENT VELOCITY', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;
%}
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% SEGMENT C - VARIATION OF DUTY CYCLE [-] - CODE IS UNFINISHED AND THEREFORE DEACTIVATED, BEFORE RUNNING THE CODE, THE CODE STRUCTURE OF SEGMENT A MUST BE TRANSFERED
%{
% first step - determination of the constant parameter (frame rate, eye movement velocity, content width, subframe number) 
% second step - definition of the range (min./max.) for and the distance (step) between the discrete values of duty cycle as the variable parameter can be defined
% third step - run section exclusively, not the whole function (!)

% CONSTANTS >> values that are not changed during the loop
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
frame_rate_hz = 30.0;
eye_velo_hor_degs = 100.0;
content_width_px = 120.0;
subframe_number = 3.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% VARIABLE >> value that is changed over the determined range below, duty_cycle needs to be a decimal number ranging from 0 to 1
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
duty_cycle_min = 0.1;
duty_cycle_max = 1.0;
duty_cycle_step = 0.1;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CONDITION (1) >> for calculation of theoretical cbu scores  
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ((duty_cycle_max-duty_cycle_min)/duty_cycle_step) ~= fix((duty_cycle_max-duty_cycle_min)/duty_cycle_step)
    
    cprintf('err','ATTENTION - lower and upper borders of duty cycle range do not match with duty cycle steps (see VARIABLE) \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CBU SCORE >> gives out the different cbu scores (option 1, option 2 and option 3) calculated by the theoretical cbu model 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for duty_cycle = (duty_cycle_min):(duty_cycle_step):(duty_cycle_max) % duty cycle range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % cbu score - option 1: area of retinal CBU stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
    
    cbu_score_opt1_table(aa,1) = duty_cycle; % first column includes frame rates        
    cbu_score_opt1_table(aa,2)= cbu_score_opt1; % second column includes theoretically calculated cbu score (mean value)
    
    % cbu score - option 2: area & intensity of retinal CBU stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
   
    cbu_score_opt2_table(aa,1) = duty_cycle;    
    cbu_score_opt2_table(aa,2)= cbu_score_opt2; 

    % cbu score - option 3: area & intensity & color of retinal CBU stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description

    % cbu_score_opt3_table(aa,1) = ... ;   
    % cbu_score_opt3_table(aa,2)= ... ; 

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% WHITE CONTENT >> gives out the score for retinal stimulation with unbiased white content (option 1, option 2 and option 3)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for duty_cycle = (duty_cycle_min):(duty_cycle_step):(duty_cycle_max) % duty cycle range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % white content - option 1: area of retinal white content stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
    
    white_content_opt1_table(aa,1) = duty_cycle;    
    white_content_opt1_table(aa,2)= white_content_opt1;
    
    % white content - option 2: area & intensity of retinal white content stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
    
    white_content_opt2_table(aa,1) = duty_cycle;        
    white_content_opt2_table(aa,2)= white_content_opt2;    

    % white content - option 3: area & intensity & color of retinal white content stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description
    
    % white_content_opt3_table(aa,1) = ... ;        
    % white_content_opt3_table(aa,2)= ... ;    
        
end

cbu_score_opt1_table, white_content_opt1_table, reference_opt1

cbu_score_opt2_table, white_content_opt2_table, reference_opt2
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% THRESHOLD CALCULATION >> gives out the values for the variable parameter (see above) for which a phase transition from phase 1 to phase 2 and from phase 2 to phase 3 is accomplished
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
duty_cycle_threshold12 = 2-(((content_width_px.*pixel_pitch_hor_deg).*subframe_number.*frame_rate_hz)./eye_velo_hor_degs)
duty_cycle_threshold23 = 1-(((content_width_px.*pixel_pitch_hor_deg).*subframe_number.*frame_rate_hz)./eye_velo_hor_degs)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (1) >> illustrates the dependency between the variable parameter (x-axis) and cbu score, white content and reference as target parameter (y-axis)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
close all; 

f1 = figure (1); hold on; % figure 1 refers to target parameter calculation with option 1 
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt1_table(1:end,1),cbu_score_opt1_table(1:end,2), '--ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt1_table(1:end,1),white_content_opt1_table(1:end,2), '--ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([duty_cycle_min duty_cycle_max],[reference_opt1 reference_opt1],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([duty_cycle_threshold12 duty_cycle_threshold12],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([duty_cycle_threshold23 duty_cycle_threshold23],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([duty_cycle_min duty_cycle_max]);
set(gca,'xTick',duty_cycle_min:(duty_cycle_max-duty_cycle_min)/10:duty_cycle_max);
xlabel('DUTY CYCLE [-]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));

ylim([0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])]);
set(gca,'yTick',0:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])/10:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1]));
ylabel('COLOR BREAK-UP SCORE [MM2]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF DUTY CYCLE', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;

f2 = figure (2); hold on; % figure 2 refers to target parameter calculation with option 2 
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt2_table(1:end,1),cbu_score_opt2_table(1:end,2), '--ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt2_table(1:end,1),white_content_opt2_table(1:end,2), '--ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([duty_cycle_min duty_cycle_max],[reference_opt2 reference_opt2],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([duty_cycle_threshold12 duty_cycle_threshold12],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([duty_cycle_threshold23 duty_cycle_threshold23],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([duty_cycle_min duty_cycle_max]);
set(gca,'xTick',duty_cycle_min:(duty_cycle_max-duty_cycle_min)/10:duty_cycle_max);
xlabel('DUTY CYCLE [-]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));

ylim([0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])]);
set(gca,'yTick',0:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])/10:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2]));
ylabel('COLOR BREAK-UP SCORE [---]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF DUTY CYCLE', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;
%}
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% SEGMENT D - VARIATION OF CONTENT WIDTH [PX] - CODE IS UNFINISHED AND THEREFORE DEACTIVATED, BEFORE RUNNING THE CODE, THE CODE STRUCTURE OF SEGMENT A MUST BE TRANSFERED
%{
% first step - determination of the constant parameter (frame rate, duty cycle, eye movement velocity, subframe number) 
% second step - definition of the range (min./max.) for and the distance (step) between the discrete values of content width as the variable parameter can be defined
% third step - run section exclusively, not the whole function (!)

% CONSTANTS >> values that are not changed during the loop
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
frame_rate_hz = 90.0;
duty_cycle = 0.30;
eye_velo_hor_degs = 100.0;
subframe_number = 3.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% VARIABLE >> value that is changed over the determined range below, content_width_px needs to be a positive integer and can not exceed display dimensions
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
content_width_px_min = 5.0;
content_width_px_max = 100.0;
content_width_px_step = 1.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CONDITION (1) >> for calculation of theoretical cbu scores  
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ((content_width_px_max-content_width_px_min)/content_width_px_step) ~= fix((content_width_px_max-content_width_px_min)/content_width_px_step)
    
    cprintf('err','ATTENTION - lower and upper borders of content width range do not match with content width steps (see VARIABLE) \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CBU SCORE >> gives out the different cbu scores (option 1, option 2 and option 3) calculated by the theoretical cbu model 
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for content_width_px = (content_width_px_min):(content_width_px_step):(content_width_px_max) % content width range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % cbu score - option 1: area of retinal CBU stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
    
    cbu_score_opt1_table(aa,1) = content_width_px; % first column includes frame rates        
    cbu_score_opt1_table(aa,2)= cbu_score_opt1; % second column includes theoretically calculated cbu score (mean value)
    
    % cbu score - option 2: area & intensity of retinal CBU stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
   
    cbu_score_opt2_table(aa,1) = content_width_px;    
    cbu_score_opt2_table(aa,2)= cbu_score_opt2; 

    % cbu score - option 3: area & intensity & color of retinal CBU stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description

    % cbu_score_opt3_table(aa,1) = ... ;   
    % cbu_score_opt3_table(aa,2)= ... ; 
     
end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% WHITE CONTENT >> gives out the score for retinal stimulation with unbiased white content (option 1, option 2 and option 3)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for content_width_px = (content_width_px_min):(content_width_px_step):(content_width_px_max) % content width range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % white content - option 1: area of retinal white content stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
    
    white_content_opt1_table(aa,1) = content_width_px;       
    white_content_opt1_table(aa,2)= white_content_opt1;
    
    % white content - option 2: area & intensity of retinal white content stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
    
    white_content_opt2_table(aa,1) = content_width_px;        
    white_content_opt2_table(aa,2)= white_content_opt2;    

    % white content - option 3: area & intensity & color of retinal white content stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description
    
    % white_content_opt3_table(aa,1) = ... ;        
    % white_content_opt3_table(aa,2)= ... ;
    
end

cbu_score_opt1_table, white_content_opt1_table, reference_opt1

cbu_score_opt2_table, white_content_opt2_table, reference_opt2
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% THRESHOLD CALCULATION >> gives out the values for the variable parameter (see above) for which a phase transition from phase 1 to phase 2 and from phase 2 to phase 3 is accomplished
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
content_width_px_threshold12 = ((eye_velo_hor_degs.*(2-duty_cycle))/(subframe_number.*frame_rate_hz))/pixel_pitch_hor_deg
content_width_px_threshold23 = ((eye_velo_hor_degs.*(1-duty_cycle))/(subframe_number.*frame_rate_hz))/pixel_pitch_hor_deg
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (1) >> illustrates the dependency between the variable parameter (x-axis) and cbu score, white content and reference as target parameter (y-axis)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
close all; 

f1 = figure (1); hold on; % figure 1 refers to target parameter calculation with option 1
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt1_table(1:end,1),cbu_score_opt1_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt1_table(1:end,1),white_content_opt1_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([content_width_px_min content_width_px_max],[reference_opt1 reference_opt1],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([content_width_px_threshold12 content_width_px_threshold12],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([content_width_px_threshold23 content_width_px_threshold23],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([content_width_px_min content_width_px_max]);
set(gca,'xTick',content_width_px_min:(content_width_px_max-content_width_px_min)/10:content_width_px_max);
xlabel('CONTENT WIDTH [PX]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])]);
set(gca,'yTick',0:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])/10:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1]));
ylabel('COLOR BREAK-UP SCORE [MM2]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF CONTENT WIDTH', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;

f2 = figure (2); hold on; % figure 2 refers to target parameter calculation with option 2
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt2_table(1:end,1),cbu_score_opt2_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt2_table(1:end,1),white_content_opt2_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([content_width_px_min content_width_px_max],[reference_opt2 reference_opt2],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([content_width_px_threshold12 content_width_px_threshold12],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([content_width_px_threshold23 content_width_px_threshold23],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([content_width_px_min content_width_px_max]);
set(gca,'xTick',content_width_px_min:(content_width_px_max-content_width_px_min)/10:content_width_px_max);
xlabel('CONTENT WIDTH [PX]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])]);
set(gca,'yTick',0:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])/10:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2]));
ylabel('COLOR BREAK-UP SCORE [---]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF CONTENT WIDTH', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;
%}
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%% SEGMENT E - VARIATION OF SUBFRAME NUMBER [-] - CODE IS UNFINISHED AND THEREFORE DEACTIVATED, BEFORE RUNNING THE CODE, THE CODE STRUCTURE OF SEGMENT A MUST BE TRANSFERED 
%{
% first step - determination of the constant parameter (frame rate, duty cycle, content width, eye movement velocity) 
% second step - definition of the range (min./max.) for and the distance (step) between the discrete values of subframe number as the variable parameter can be defined
% third step - run section exclusively, not the whole function (!)
%
% ATTENTION - before running this section some code lines within CBU_MODEL have to be deactivated, without deactivation the calculations within this section can not be executed ...
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% 62 >> frame_rate_hz
% 65 >> subframe_number
% 69 >> duty_cycle
% 174 >> content_width_px
% 212 >> eye_velo_hor_degs
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% 96 - 164 >> data entry (lighting parameter)
% 320 - 390 >> robustness (lighting parameter)
% 559 - 562 >> sf_intensity_table
% 571 - 574 >> sf_color_table
% 1428 - 1543 >> temporal summation (1)
% 1547 - 1575 >> temporal summation (2)
% 1793 - end >> all plotting
% -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% ... it is possible that the code lines change since work on CBU_MODEL is in progress (05/05/2020)


% CONSTANTS >> values that are not changed during the loop
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
frame_rate_hz = 20.0;
duty_cycle = 0.30;
eye_velo_hor_degs = 150.0;
content_width_px = 120.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% VARIABLE >> value that is changed over the determined range below, subframe_number needs to be a positive integer from 2 to 6
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subframe_number_min = 2.0;
subframe_number_max = 6.0;
subframe_number_step = 1.0;
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CONDITION (1) >> for calculation of theoretical cbu scores  
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ((subframe_number_max-subframe_number_min)/subframe_number_step) ~= fix((subframe_number_max-subframe_number_min)/subframe_number_step)
    
    cprintf('err','ATTENTION - lower and upper borders of subframe number range do not match with subframe number steps (see VARIABLE) \n') 
    return

end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% CBU SCORE >> gives out the different cbu scores (option 1, option 2 and option 3) calculated by the theoretical cbu model
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for subframe_number = (subframe_number_min):(subframe_number_step):(subframe_number_max) % subframe number range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % cbu score - option 1: area of retinal CBU stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
    
    cbu_score_opt1_table(aa,1) = subframe_number; % first column includes frame rates        
    cbu_score_opt1_table(aa,2)= cbu_score_opt1; % second column includes theoretically calculated cbu score (mean value)
    
    % cbu score - option 2: area & intensity of retinal CBU stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
   
    cbu_score_opt2_table(aa,1) = subframe_number;    
    cbu_score_opt2_table(aa,2)= cbu_score_opt2; 

    % cbu score - option 3: area & intensity & color of retinal CBU stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description

    % cbu_score_opt3_table(aa,1) = ... ;   
    % cbu_score_opt3_table(aa,2)= ... ;
     
end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% WHITE CONTENT >> gives out the score for retinal stimulation with unbiased white content (option 1, option 2 and option 3)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
aa = 0;

for subframe_number = (subframe_number_min):(subframe_number_step):(subframe_number_max) % subframe number range of choice tested (see VARIABLE above)
    
    aa = aa + 1;   
    
    [pixel_pitch_hor_deg,cbu_score_opt1,white_content_opt1,reference_opt1,cbu_score_opt2,white_content_opt2,reference_opt2] = cbu_model(eye_velo_hor_degs,frame_rate_hz,duty_cycle,content_width_px,subframe_number);

    % white content - option 1: area of retinal white content stimulation (no light intensity considered, no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> surface of stimulation) for a more detailed description
    
    white_content_opt1_table(aa,1) = subframe_number;        
    white_content_opt1_table(aa,2)= white_content_opt1;
    
    % white content - option 2: area & intensity of retinal white content stimulation (no light color considered) >> see confluence (theoretical cbu model >> cbu scoring >> volume of stimulation) for a more detailed description
    
    white_content_opt2_table(aa,1) = subframe_number;        
    white_content_opt2_table(aa,2)= white_content_opt2;    

    % white content - option 3: area & intensity & color of retinal white content stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description
    
    % white_content_opt3_table(aa,1) = ... ;        
    % white_content_opt3_table(aa,2)= ... ;    
   
end

cbu_score_opt1_table, white_content_opt1_table, reference_opt1

cbu_score_opt2_table, white_content_opt2_table, reference_opt2
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% THRESHOLD CALCULATION >> gives out the values for the variable parameter (see above) for which a phase transition from phase 1 to phase 2 and from phase 2 to phase 3 is accomplished
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
subframe_number_threshold12 = ((eye_velo_hor_degs.*(2-duty_cycle))/(content_width_px.*frame_rate_hz))/pixel_pitch_hor_deg
subframe_number_threshold23 = ((eye_velo_hor_degs.*(1-duty_cycle))/(content_width_px.*frame_rate_hz))/pixel_pitch_hor_deg
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% GRAPHICS (1) >> illustrates the dependency between the variable parameter (x-axis) and cbu score, white content and reference as target parameter (y-axis)
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
close all; 

f1 = figure (1); hold on; % figure 1 refers to target parameter calculation with option 1
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt1_table(1:end,1),cbu_score_opt1_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt1_table(1:end,1),white_content_opt1_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([subframe_number_min subframe_number_max],[reference_opt1 reference_opt1],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([subframe_number_threshold12 subframe_number_threshold12],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([subframe_number_threshold23 subframe_number_threshold23],[0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([subframe_number_min subframe_number_max]);
set(gca,'xTick',subframe_number_min:(subframe_number_max-subframe_number_min)/10:subframe_number_max);
xlabel('SUBFRAME NUMBER [-]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])]);
set(gca,'yTick',0:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1])/10:max([cbu_score_opt1_table(1:end,2);white_content_opt1_table(1:end,2);reference_opt1]));
ylabel('COLOR BREAK-UP SCORE [MM2]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF SUBFRAME NUMBER', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;

f2 = figure (2); hold on; % figure 2 refers to target parameter calculation with option 2
set(gca,'FontSize',11);

p(1)=plot(cbu_score_opt2_table(1:end,1),cbu_score_opt2_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','b','MarkerSize',3.0); % plot cbu score
p(2)=plot(white_content_opt2_table(1:end,1),white_content_opt2_table(1:end,2), '-ko','LineWidth', 1.0,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',3.0); % plot white content
p(3)=line([subframe_number_min subframe_number_max],[reference_opt2 reference_opt2],'Color',[0.0,1.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot reference

p(4)=line([subframe_number_threshold12 subframe_number_threshold12],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 1 to phase 2)
p(5)=line([subframe_number_threshold23 subframe_number_threshold23],[0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])],'Color',[1.0,0.0,0.0],'LineStyle','--','LineWidth', 1.25); % plot phase threshold (phase 2 to phase 3)

xlim([subframe_number_min subframe_number_max]);
set(gca,'xTick',subframe_number_min:(subframe_number_max-subframe_number_min)/10:subframe_number_max);
xlabel('SUBFRAME NUMBER [-]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.1f'));

ylim([0 max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])]);
set(gca,'yTick',0:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2])/10:max([cbu_score_opt2_table(1:end,2);white_content_opt2_table(1:end,2);reference_opt2]));
ylabel('COLOR BREAK-UP SCORE [---]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];
title ('CBU SCORE IN DEPENDENY OF SUBFRAME NUMBER', 'Fontsize', 20)
set(gcf,'WindowStyle','docked');
hold off;
%}
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

end