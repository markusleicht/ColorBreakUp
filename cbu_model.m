function cbu_model

% AUTHOR: Markus Leicht
% EMAIL: markusleicht84@gmx.de
% DATA: 2021-02-16  
% VERSION: 1.0.0.0
% DOWNLOAD: https://github.com/markusleicht/ColorBreakUp 
%
% BACKGROUND:
% The function CBU_MODEL was written during the author's doctorate. The associated doctoral thesis refers to the model and provides further information on the subject of color break-up (CBU): Leicht, M. (2022). Perception of Color Break-Up [Doctoral thesis]. Ilmenau University of Technology).
%
% PURPOSE:
% The CBU_MODEL aspires to forecast a viewer's perception of CBU on an objectified basis for predefined scenarios. The model scenario always includes (1) a field-sequential color (FSC) display system presenting digital content and (2) a viewer observing this content. Under certain conditions, such a scenario potentially leads to the occurrence of CBU effects (see Chapter 2 in the associated doctoral thesis). The CBU effects are computed by the CBU_MODEL in a multi-stage process. 
%
% BENEFIT:
% Modeling of potentially CBU-provoking scenarios is advantageous in many ways. First, it is beneficial to assess CBU effects quantitatively and objectively, in contrast to expensive empirical research with its technical and participant-related limitations. Second, both the processes involved in CBU genesis and the impact of certain determinants on CBU characteristics can be determined simply and straighforwardly at the push of a button. And third, the model can be used to establish a classification scheme for CBU variants with different color patterns.
%
% STRUCTURE:
% Segment I - General Model Assumptions
% Segment II - Data Entry Mask
% Segment III - Experimental Category
% Segment IV - Data Entry Rules
% Segment V - Subframe Characteristics
% Segment VI - Color Break-Up Phase
% Segment VII - Phase Transition Thresholds
% Segment VIII - Intersections and Unions
% Segment IX - Color Sectors
% Segment X - Temporal Summation
% Segment XI - Model Indices
% Segment XII - Overview Table
% Segment XIII - Graphical Illustration
%
% INPUT (SEGMENT II - DATA ENTRY MASK):
% Model input parameters must be determined to specify the frame conditions of the model scenario (Segment II: Data Entry Mask; current model input refers to the example of choice that is discussed in Chapter 3 of the associated doctoral thesis). The display unit's technical specifications, the characteristics of the presented content, the viewer and his/her behavior (e.g. eye movements), and other surrounding conditions need to be defined in detail. Based on the model input, light intensity and color characteristics (hue and saturation) of the displayed content (physical stimulus) that potentially provokes CBU are quantified. Furthermore, the content's size, its position at display level, and the exposure time of light emission are defined. Taking the underlying FSC principle of content creation into account, the content characterization is executed at subframe level, i.e., every single subframe is separately specified at first. This allows deconstruction of the CBU genesis to its starting point. Subsequently, the interactions of the individual subframes within a full frame cycle are computed to create a realistic copy of the displayed content. 
%
% DATA OUTPUT (SEGMENT XII - OVERVIEW TABLE):
% The corresponding stimulation of the viewer is calculated by the CBU_MODEL under consideration of the determined content characteristics, leading to the exact description of retinal stimulation within the viewer's eye. The model predicts the retinal zone of stimulation and the intensity and color characteristic profile within this area. This eventually allows to distinguish between CBU-free and CBU-loaded retinal zones by putting out a model-based replication of the perceived scene. In other words: The model output predicts a real viewer's visual impression that is provoked by the defined scenario.
% After running the function, the model's data output is summarized within an overview table (Segment XII) in the command window. The table lists the classification of experimental category and cbu phase and specifies the thresholds for cbu phase transitions (phase 1 to phase 2 to phase 3). Furthermore, the model indices that compress the complexity of the model calculations to single numbers are listed (cbu score, non-cbu score, reference score). The model input parameters are listed as well - giving an overview of the model's input and ouput.
%
% GRAPHIC OUTPUT (SEGMENT XIII - GRAPHICAL ILLUSTRATION):
% Figure 1: Space-Time Plot
% Figure 2: Space-Time-Intensity Plot
% Figure 3: Intensity Profile
% Figure 4: Intensity Map
% Figure 5: Color Profile
% Figure 6: Color Map (Overview)
% Figure 7: Color Map (Detail)
%
% DESCRIPTION:
% A basic code description is embedded in the MATLAB function itself. However, the fundamental structure of the model and its build-up are described in the corresponding doctoral thesis (see Chapter 3).
%
% COMPATIBILITY:
% The function was written with MATLAB Version R2018b (OS: Windows 10, 64bit).      
%
% MATLAB SEARCH PATH
% For the functionality of the CBU_MODEL it is not necessary to add further folders to the search path. External code from Xavier Beudaert for the calculation of intersection sets (line 809-873) and from David Goodmanson for the calculation of union sets (line 994-1015) is directly implemented into the code of the CBU_MODEL.
%
% Copyright 2021 Markus Leicht

%% SEGMENT I - GENERAL MODEL ASSUMPTIONS
%
% COORDINATE SYSTEMS
% - the coordinate system of the display unit has its point of origin in the horizontal/vertical center of the display unit, mathematically positive position values are assigned in rightward respectively upward direction and mathematically negative position values are assigned in leftward respectively downward direction
% - the point of origin of the observers's retinal coordinate system coincides with the fovea centralis, the image inversion on the retina leads to mathematically positive position values in leftward/downward direction and mathematically negative position values in rightward/upward direction
%
% OBSERVER & DISPLAY POSITION 
% - the observer's visual system is represented by a cyclops eye whose fovea centralis coincides with the horizontal/vertical center of the display unit in primary eye position 
% - the display unit is adjusted without any tilt refering the observer's head position
%
% EYE & CONTENT POSITION
% - eye respectively content positions can be determined in horizontal and vertical direction within the limits of the display unit
% - for the horizontal dimension eye and content position have to be defined by determination of two positions - the start and stop point of the eye respectivelly content movement
% - for the vertical dimension eye and content position has to be defined by determination of only one stable position of eye and content since there is no movement in vertical direction
%
% EYE & CONTENT MOVEMENT
% - only left to right movement in horizontal direction for eye and content movement (no right to left movement, no vertical movement, no transversal movement)
% - constant movement velocities for content and eye movement (no positive or negative accelaration)
% - eye and content movement do always start at the same time - no delay (important to determine the number of frame cycles that have to be calculated)
% - content disappears the moment the stop point of the movement path is reached (important to determine the number of frame cycles that have to be calculated)
% - eye movement velocity in [DEG/SEC] is defined as rotation around the eyes center of rotation (Z')
%
% EYE ANATOMY 
% - parameter of the standard observer eye refer to gullstrand-emsley schematic human eye model (Emsley, 1936)
% - retina of the standard observer is assumed to have a simplified spherical shape of r = 12.0mm (Smith & Atchison, 1997, S.678)
% - optical axis, visual axis and fixation line are assumed to be identical, meaning that angle epsilon (btw. optical axis and visual axis) and also angle eta (btw. optical axis and fixation line) are assumed to be zero
%
% OPTICAL IMAGING
% - transition from content position on display / screen in [DEG] to retinal position in [DEG] via nodal points N and N'
% - calculation of retinal image size in [MM] via display-to-retina transition correction
%
% COLOR AND LIGHT INTENSITY
% - light color of all single SF is generally described via chromaticity coordinates (x,y) 
% - light intensity of all single SF is generally described in a physical manner in radiance (Le) in [W/(SR.*M2)] and as well in a photometric manner in luminance (Lv) in [CD/M2], the choice of radiance / luminance is based on the assumption of a projection setup (projector and screen)  
% - light color and intensity of the emerging mixed color (e.g. white) resulting from blending the single SFs (e.g. red, green and blue) is calculated by transitioning the chromaticity coordinates (x,y) and luminance (Lv or Y) of the single SF via tristimulus values (XYZ) to the chromaticity coordinates (x,y) and lumiance (Lv or Y) of the mixed color
%
% COGNITIVE PROCESSING
% - the time threshold for temporal summation of visual stimuli described by Bloch's Law is not determined at a fix value - e.g. 0.05sec for cones and 0.10sec for rods (Hood and Finkelstein, as cited in Blake and Sekuler, 2006, p.98) - instead, temporal summation is executed for all subframes of one full frame cycle only (frame cycle with the largest CBU effect is chosen) in order to calculate the color characteristic and light intensity profile of the resulting retinal stimulation  
% - the areal threshold for spatial summation of visual stimuli described by Ricco's Law is not considered within the model since the characteristic of the applied stimuli (position and size) within the model exceeds the limits for which Ricco's Law is valid 
%
% DEPTH OF MODEL 
% - model calculations of the resulting light intensity and color pattern of retinal stimulation mainly consider geometrical optics for position calculation in combination with the process of temporal summation of light stimuli (one full frame cycle is integrated, for the final model output the number of frame cycles presented for the defined conditions is also taken into account)
% - additional processing steps - e.g. variability of retinal sensitivity from fovea to periphery or cognitive processing regarding saccadic suppression - which lead to the actual subjective perception of cbu are not considered so far

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

frame_rate_hz = 90.0; % number of frames displayed per second
frame_duration_sec = 1./frame_rate_hz; % duration of one single frame cycle

subframe_number = 3; % number of subframes displayed within on single frame cycle, specification of subframes (light intensity & color, see below) is designed for up to 6 subframes (could be extended if necessary)
subframe_rate_hz = frame_rate_hz.*subframe_number; % duration of one single subframe cycle, x-times (depending on number of subframes) the frame rate of the display
subframe_duration_sec = 1./subframe_rate_hz; % duration of one single subframe (on and off time included)

duty_cycle = 0.30; % on-off-ratio of one single subframe cycle (on-time divided by complete subframe-time)

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

content_width_px = 40.0; % horizontal dimension of content, content is assumed to have a simple quadratic or rectangular shape (e.g. vertical bar)

content_velo_hor_pxframe = 8.0; % content movement velocity [PX/FR] in horizontal direction 

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

eye_velo_hor_pxframe = 8.0; % eye movement velocity [PX/FR] in horizontal direction    

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
    frame_number = eye_distance_hor_px./eye_velo_hor_pxframe+1; % frame number within CAT2 is calculated on basis of the time eye movement is executed (since no content movement is executed at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included  

elseif content_velo_hor_pxframe > 0 && eye_velo_hor_pxframe > 0
    cat = 3;
    
    if (content_distance_hor_px/content_velo_hor_pxframe) >= (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is equal or longer than the time of eye movement (eye_distance_px/eye_velo_pxframe), than the frame number to be considered is calculated on basis of the time eye movement is executed, even when content movement is executed for a longer period of time eye movement time is still basis for calculation since CBU only occurs when eye movement is in progress, otherwise it will be CAT1 (eye fix, content variable) that does not provoke CBU at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included, statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning) 
        frame_number = eye_distance_hor_px./eye_velo_hor_pxframe+1;
    elseif (content_distance_hor_px/content_velo_hor_pxframe) < (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is shorter than the time of eye movement (eye_distance_px/eye_velo_pxframe) than the frame number is calculated on basis of the time content movement is executed since the content disappears directly after reaching the stop point of its movement path (without presentation of content no CBU can occur); calculation of time in [FRAME] by the term (content_distance_px./content_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning)
        frame_number = content_distance_hor_px./content_velo_hor_pxframe+1;
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
    
elseif ~isscalar(content_velo_hor_pxframe) || content_velo_hor_pxframe < 0 || content_velo_hor_pxframe ~= fix(content_velo_hor_pxframe)
    error('ATTENTION - content_velo_hor_pxframe needs to be a positive integer (including zero)');

elseif ~isscalar(eye_position_ver_px) || eye_position_ver_px ~= fix(eye_position_ver_px) || abs(eye_position_ver_px) > resolution_ver_px./2
    error('ATTENTION - eye_position_ver_px needs to be a integer and can not exceed display dimensions');    
    
elseif ~isscalar(eye_start_hor_px) || eye_start_hor_px ~= fix(eye_start_hor_px) || abs(eye_start_hor_px) > resolution_hor_px./2
    error('ATTENTION - eye_start_hor_px needs to be a integer and can not exceed display dimensions');

elseif ~isscalar(eye_stop_hor_px) || eye_stop_hor_px ~= fix(eye_stop_hor_px) || abs(eye_stop_hor_px) > resolution_hor_px./2
    error('ATTENTION - eye_stop_hor_px needs to be a integer and can not exceed display dimensions');    
    
elseif ~isscalar(eye_velo_hor_pxframe) || eye_velo_hor_pxframe < 0 || eye_velo_hor_pxframe ~= fix(eye_velo_hor_pxframe)
    error('ATTENTION - eye_velo_hor_pxframe needs to be a positive integer (including zero)');
       
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
    
    % transfer stable [PX]/[PX/FR] values regarding width, movement distance and velocity of content respectively movement distance and velocity of eye to the corresponding variable display position dependent [DEG]/[DEG/SEC] values (before starting to calculate basic points A, C, E and G)
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

if subframe_location_table_hor_x1(1,eye_velo_hor_variable_degs_table_index) == subframe_location_table_hor_x1(5,eye_velo_hor_variable_degs_table_index) % proof assumption that A(SF1) is equal to A(SF2) for verification of PH0; since a movement pattern with constant retinal velocity in [PX/FR] (velocity values are constant for complete sequence) respectively [DEG/SEC] (at least constant velocity values within a frame cycle) is presumed (see SEGMENT I), it is sufficient to check the assumption for the first two subframes of the chosen frame cycle only; it should generally be noted that (at present) it is not possible to confirm this if loop statement that would result in a classification within PH0 since the classification of the experimental category in the code further above (see SEGMENT III) stops the code from running when CAT0 or CAT1 are determined which are the only categories that result in PH0, for reasons of completeness the if loop refering to PH0 is not excluded from the code;
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
                    
                    % Copyright (c) 2011, Xavier Xavier
                    % All rights reserved.
                    % Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
                    % * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
                    % * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution
                    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
    
    % calculation of frame number that is not rounded is deactivated since only positive integers are allowed as frame numbers for code execution (at present) >> see SEGMENT IV    
    %{
    if cat == 2 % final cbu_score_opt2 (also non_cbu_score_opt2 and reference_score_opt2, see below) has to be calculated by multiplying with frame_number_not_rounded in order to make conditions with variable frame rates comparable by scaling the cbu score up (see confluence for detailed description and explanation)
    
        frame_number_not_rounded = (eye_distance_hor_px./eye_velo_hor_pxframe)+1; % frame number within CAT2 is calculated on basis of the time eye movement is executed (since no content movement is executed at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; calculated frame number is not rounded off (like it is done in the section DETERMINATION OF EXPERIMENTAL CATEGORY) in order to calculate the exact multiplication factor

    elseif cat == 3
        
        if (content_distance_hor_px/content_velo_hor_pxframe) >= (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is equal or longer than the time of eye movement (eye_distance_px/eye_velo_pxframe) than the frame number to be considered is calculated on basis of the time eye movement is executed (even when content movement is executed for a longer period of time eye movement time is still basis for calculation since CBU only occurs when eye movement is in progress, otherwise it will be CAT1 (eye fix, content variable) that does not provoke CBU at all); calculation of time in [FRAME] by the term (eye_distance_px./eye_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; calculated frame number is not rounded off (like it is done in the section DETERMINATION OF EXPERIMENTAL CATEGORY) in order to calculate the exact multiplication factor; statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning) 
        
            frame_number_not_rounded = (eye_distance_hor_px./eye_velo_hor_pxframe)+1;
        
        elseif (content_distance_hor_px/content_velo_hor_pxframe) < (eye_distance_hor_px/eye_velo_hor_pxframe) % if the time of content movement (content_distance_px/content_velo_pxframe) within CAT3 is shorter than the time of eye movement (eye_distance_px/eye_velo_pxframe) than the frame number is calculated on basis of the time content movement is executed (since the content disappears directly after reaching the stop point of its movement path, without presentation of content no CBU can occur); calculation of time in [FRAME] by the term (content_distance_px./content_velo_pxframe) is always added up by the value of 1 since the first frame in start position has to be included; calculated frame number is not rounded off (like it is done in the section DETERMINATION OF EXPERIMENTAL CATEGORY) in order to calculate the exact multiplication factor (in this case the frame number is generally expected to be an integer since content distance and content movement velocity are expected to be positive integers); statements regarding the determination of the number of frames that have to be calculated are made under the assumptions that eye and content start their movement at the same time without delay & that the content disappears directly after the stop point is reached (see general assumptions at the beginning)
            
            frame_number_not_rounded = (content_distance_hor_px./content_velo_hor_pxframe)+1;
    
        end
        
    end
    %}   
    
    cbu_score_opt2_mm2 = cbu_front_surface_mm.*cbu_height_opt2_mm.*frame_number; % final cbu score for option 2 (refering to [MM2] expression of retinal stimulation) is calculated by multiplying 3D body front surface and height of retinal cbu stimulation (and also frame number, see explanation in code description above)
    cbu_score_opt2_deg2 = cbu_front_surface_deg.*cbu_height_opt2_deg.*frame_number; % final cbu score for option 2 (refering to [DEG2] expression of retinal stimulation) is calculated by multiplying 3D body front surface and height of retinal cbu stimulation (and also frame number, see explanation in code description above)
    
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
    
    non_cbu_score_opt2_mm2 = non_cbu_front_surface_mm.*non_cbu_height_opt2_mm.*frame_number; % final non-cbu score for retinal stimulation with unbiased content (refering to [MM2] expression of retinal stimulation) is determined by multiplying 3D body front surface and height (and also frame number, see explanation in code description above)
    non_cbu_score_opt2_deg2 = non_cbu_front_surface_deg.*non_cbu_height_opt2_deg.*frame_number; % final non-cbu score for retinal stimulation with unbiased content (refering to [DEG2] expression of retinal stimulation) is determined by multiplying 3D body front surface and height (and also frame number, see explanation in code description above)     
    
    % REFERENCE_SCORE >> gives out the volume of stimulation for unbiased content if no relative retinal movement was existent (refering to [MM2] respectively [DEG2] expression for retinal stimulation)
    
    reference_height_opt2_mm = subframe_location_table_ver_mm(1,1)-subframe_location_table_ver_mm(2,1); % vertical dimension of the retinal stimulation with unbiased content under the made reference assumption in [MM]
    reference_height_opt2_deg = subframe_location_table_ver_x1(1,1)-subframe_location_table_ver_x1(2,1); % vertical dimension of the retinal stimulation with unbiased content under the made reference assumption in [DEG]
    
    reference_front_surface_mm = (subframe_location_table_hor_mm(1,eye_velo_hor_variable_degs_table_index)-subframe_location_table_hor_mm(2,eye_velo_hor_variable_degs_table_index)).*(subframe_on_sec.*sum(sf_luminance)); % 3D body front surface of unbiased content provocation under the made reference assumption (refering to [MM] expression for retinal stimulation)
    reference_front_surface_deg = (subframe_location_table_hor_x1(1,eye_velo_hor_variable_degs_table_index)-subframe_location_table_hor_x1(2,eye_velo_hor_variable_degs_table_index)).*(subframe_on_sec.*sum(sf_luminance)); % 3D body front surface of unbiased content provocation under the made reference assumption (refering to [DEG] expression for retinal stimulation)
    
    reference_score_opt2_mm2 = reference_front_surface_mm.*reference_height_opt2_mm.*frame_number; % final reference score (refering to [MM2] expression for retinal stimulation) is calculated by multiplying 3D body front surface and height of the assumed unbiased content for the defined conditions (and also frame number, see explanation in code description above) 
    reference_score_opt2_deg2 = reference_front_surface_deg.*reference_height_opt2_deg.*frame_number; % final reference score (refering to [DEG2] expression for retinal stimulation) is calculated by multiplying 3D body front surface and height of the assumed unbiased content for the defined conditions (and also frame number, see explanation in code description above)

% OPTION 3 - area & intensity & color of retinal stimulation >> see confluence (theoretical cbu model >> cbu scoring >> color factor) for a more detailed description

    % CBU_SCORE >> ...

    % NON_CBU_SCORE >> ...
    
    % REFERENCE_SCORE >> ...    

%% SEGMENT XII - OVERVIEW TABLE (INCLUDING PARAMETERS OF PRESENT SEQUENCE CONDITIONS AND RESULTING MODEL OUTPUT)

CATEGORY = {'classification';'---';'scanning';'---';'---';'---';'---';'timing';'---';'---';'---';'---';'---';'---';'lighting';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'color';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'hardware';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'content';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'human eye';'---';'---';'---';'---';'---';'---';'---';'---';'---';'viewing';'---';'---';'---';'anatomy';'---';'---';'---';'---';'---';'---';'---';'cognition';'---';'---';'phase threshold';'---';'---';'---';'model index';'---';'---';'---';'---';'---'}; 

PARAMETER = {'experimental category';'cbu phase';'integrated frame cycle';'integrated subframes';'scan resolution';'scan start position';'scan stop position';'frame number';'frame rate';'frame duration';'subframe number';'subframe rate';'subframe duration';'duty cycle';'sf1 radiance';'sf2 radiance';'sf3 radiance';'sf4 radiance';'sf5 radiance';'sf6 radiance';'sf1 luminance';'sf2 luminance';'sf3 luminance';'sf4 luminance';'sf5 luminance';'sf6 luminance';'sf1 x-coordinate';'sf2 x-coordinate';'sf3 x-coordinate';'sf4 x-coordinate';'sf5 x-coordinate';'sf6 x-coordinate';'sf1 y-coordinate';'sf2 y-coordinate';'sf3 y-coordinate';'sf4 y-coordinate';'sf5 y-coordinate';'sf6 y-coordinate';'resolution (horizontal)';'resolution (vertical)';'aspect ratio (horizontal)';'aspect ratio (vertical)';'diagonal dimension';'horizontal dimension';'vertical dimension';'pixel pitch (horizontal)';'---';'pixel pitch (vertical)';'---';'height';'---';'vertical position';'---';'width';'---';'horizontal start position';'---';'horizontal stop position';'---';'horizontal movement distance';'---';'horizontal movement velocity';'---';'vertical position';'---';'horizontal start position';'---';'horizontal stop position';'---';'horizontal movement distance';'---';'horizontal movement velocity';'---';'viewing distance';'diagonal viewing angle';'horizontal viewing angle';'vertical viewing angle';'CV – N (relax)';'CV – N2 (relax)';'CV - F2 (relax)';'CV – N (acco)';'CV – N2 (acco)';'CV - F2 (acco)';'CV - Z';'RR';'bloch time (cones)';'bloch time (rods)';'ricco space';'transition PH1-PH2';'---';'transition PH2-PH3';'---';'cbu score (option 1}';'non-cbu score (option 1)';'reference score (option 1)';'cbu score (option 2}';'non-cbu score (option 2)';'reference score (option 2)'};

VALUE = [cat;cbu_phase;eye_velo_hor_variable_degs_table_index;integrated_subframes;scan_resolution;scan_start;scan_stop;frame_number;frame_rate_hz;frame_duration_sec;subframe_number;subframe_rate_hz;subframe_duration_sec;duty_cycle;sf1_radiance;sf2_radiance;sf3_radiance;sf4_radiance;sf5_radiance;sf6_radiance;sf1_luminance;sf2_luminance;sf3_luminance;sf4_luminance;sf5_luminance;sf6_luminance;sf1_xcoordinate;sf2_xcoordinate;sf3_xcoordinate;sf4_xcoordinate;sf5_xcoordinate;sf6_xcoordinate;sf1_ycoordinate;sf2_ycoordinate;sf3_ycoordinate;sf4_ycoordinate;sf5_ycoordinate;sf6_ycoordinate;resolution_hor_px;resolution_ver_px;display_aspect_ratio_hor;display_aspect_ratio_ver;display_dia_m.*1000;display_hor_m.*1000;display_ver_m.*1000;pixel_pitch_hor_deg;pixel_pitch_hor_m.*1000;pixel_pitch_ver_deg;pixel_pitch_ver_m.*1000;content_height_deg;content_height_px;content_position_ver_deg;content_position_ver_px;max(content_width_variable_deg_table);content_width_px;content_start_hor_deg;content_start_hor_px;content_stop_hor_deg;content_stop_hor_px;content_distance_hor_deg;content_distance_hor_px;max(content_velo_hor_variable_degs_table);content_velo_hor_pxframe;eye_position_ver_deg;eye_position_ver_px;eye_start_hor_deg;eye_start_hor_px;eye_stop_hor_deg;eye_stop_hor_px;eye_distance_hor_deg;eye_distance_hor_px;max(eye_velo_hor_variable_degs_table);eye_velo_hor_pxframe;viewing_distance_m.*1000;viewing_angle_dia_deg;viewing_angle_hor_deg;viewing_angle_ver_deg;CV_N_relax_m.*1000;CV_N2_relax_m.*1000;CV_F2_relax_m.*1000;CV_N_acc_m.*1000;CV_N2_acc_m.*1000;CV_F2_acc_m.*1000;CV_Z_m.*1000;RR_m.*1000;bloch_time_cones_sec;bloch_time_rods_sec;ricco_space_mm;transition_ph1_ph2_degs;transition_ph1_ph2_pxframe;transition_ph2_ph3_degs;transition_ph2_ph3_pxframe;cbu_score_opt1_deg2;non_cbu_score_opt1_deg2;reference_score_opt1_deg2;cbu_score_opt2_deg2;non_cbu_score_opt2_deg2;reference_score_opt2_deg2];

UNIT = {'---';'---';'---';'---';'deg';'deg';'deg';'---';'hz';'sec';'---';'hz';'sec';'---';'w/(sr.*m2)';'w/(sr.*m2)';'w/(sr.*m2)';'w/(sr.*m2)';'w/(sr.*m2)';'w/(sr.*m2)';'cd/m2';'cd/m2';'cd/m2';'cd/m2';'cd/m2';'cd/m2';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'---';'px';'px';'---';'---';'mm';'mm';'mm';'deg';'mm';'deg';'mm';'deg';'px';'deg';'px';'deg';'px';'deg';'px';'deg';'px';'deg';'px';'deg/sec';'px/frame';'deg';'px';'deg';'px';'deg';'px';'deg';'px';'deg/sec';'px/frame';'mm';'deg';'deg';'deg';'mm';'mm';'mm';'mm';'mm';'mm';'mm';'mm';'sec';'sec';'mm';'deg/sec';'px/frame';'deg/sec';'px/frame';'deg2';'deg2';'deg2';'---';'---';'---'};

ANNOTATION = {'experimental category defines movement pattern of eye and content, CAT0 = no content movement and no eye movement existent, CAT1 = content movement but no eye movement existent, CAT2 = eye movement but no content movement existent, CAT3 = content movement as well as eye movement existent';'CBU phase defines cbu characteristics, PH0 = no cbu (no content and/or eye movement), PH1 = cbu occurs with the following color pattern: B/C/W/Y/R, PH2 = cbu occurs (B/C/G/Y/R), PH3 = cbu occurs (B/GAP/G/GAP/R), determination of CBU phase explicitly refers to the frame cycle with the highest eye movement velocity in [DEG/SEC] which is also the frame cycle with the largest cbu effect (variable angular velocity values because of the flat screen characteristic of the display unit)';'frame cycle that is taken into account for temporal summation (=frame cycle for which the highest eye movement velocity in [DEG/SEC] and therefore the largest CBU effect occurs, variable angular velocity values because of the flat screen characteristic of the display unit), frame cycles that are displayed under the defined conditions are numbered consecutively';'number of subframes that are included into calculation for temporal summation';'distance between retinal scanning spots for temporal summation';'horizontal start point for scanning refering to temporal summation';'horizontal stop point for scanning refering to temporal summation';'number of frames presented during current movement pattern';'repetition rate for frame presentation cycle';'time period for one single frame presentation cycle (incl. all on and off times of all included subframes)';'number of subframes included within one frame';'repetition rate for subframe presentation cycle';'time period for one single subframe presentation cycle (incl. on and off time)';'ratio between on time and on + off time within a time period of one single subframe presentation cycle';'radiance of subframe 1';'radiance of subframe 2';'radiance of subframe 3';'radiance of subframe 4';'radiance of subframe 5';'radiance of subframe 6';'luminance of subframe 1';'luminance of subframe 2';'luminance of subframe 3';'luminance of subframe 4';'luminance of subframe 5';'luminance of subframe 6';'chromaticity coordinate x of subframe 1';'chromaticity coordinate x of subframe 2';'chromaticity coordinate x of subframe 3';'chromaticity coordinate x of subframe 4';'chromaticity coordinate x of subframe 5';'chromaticity coordinate x of subframe 6';'chromaticity coordinate y of subframe 1';'chromaticity coordinate y of subframe 2';'chromaticity coordinate y of subframe 3';'chromaticity coordinate y of subframe 4';'chromaticity coordinate y of subframe 5';'chromaticity coordinate y of subframe 6';'native resolution of displaying unit (e.g. projector) in horizontal direction';'native resolution of displaying unit (e.g. projector) in vertical direction';'horizontal aspect ratio of displaying unit (e.g. projector)';'vertical aspect ratio of displaying unit (e.g. projector)';'diagonal dimension of the displaying unit (e.g. projector)';'horizontal dimension of the displaying unit (e.g. projector)';'vertical dimension of the displaying unit (e.g. projector)';'angular distance of two horizontal neighbour pixels measured from both pixel centers (valid for two horizontally centered pixels of displaying unit, reduction of angular pixel pitch in periphery)';'distance of two horizontal neighbour pixels measured from both pixel centers';'angular distance of two vertical neighbour pixels measured from both pixel centers (valid for two vertically centered pixels of displaying unit, reduction of angular pixel pitch in periphery)';'distance of two vertical neighbour pixels measured from both pixel centers';'angular dimension of cbu stimulating content in vertical direction, content is assumed to have a simple quadratic or rectangular shape (e.g. vertical bar)';'dimension of cbu stimulating content in vertical direction, content is assumed to have a simple quadratic or rectangular shape (e.g. vertical bar)';'stable vertical (angular) position of upper edge of the content during execution of the movement pattern, point of origin is the vertical center of the display (down - negative, up - positive)';'stable vertical position of upper edge of the content during execution of the movement pattern, point of origin is the vertical center of the display (down - negative, up - positive)';'horizontal angular dimension of cbu stimulating content (display position dependent, maximum value during sequence), content is assumed to have a simple quadratic or rectangular shape (e.g. vertical bar)';'horizontal dimension of cbu stimulating content, content is assumed to have a simple quadratic or rectangular shape (e.g. vertical bar)';'horizontal angular start point of content movement, position of left edge of the content at time 0, point of origin is the horizontal center of the display (left - negative, right - positive)';'horizontal start point of content movement, position of left edge of the content at time 0, point of origin is the horizontal center of the display (left - negative, right - positive)';'horizontal angular stop point of content movement, position of left edge of the content after movement, point of origin is the horizontal center of the display (left - negative, right - positive)';'horizontal stop point of content movement, position of left edge of the content after movement, point of origin is the horizontal center of the display (left - negative, right - positive)';'angular distance of content movement in horizontal direction';'distance of content movement in horizontal direction';'angular velocity of horizontal content movement, display position dependent parameter, maximum value during sequence (positive values - left to right movement, negative values - right to left movement)';'velocity of horizontal content movement (positive values - left to right movement, negative values - right to left movement)';'stable vertical (angular) position of viewers eye during movement pattern (in horizontal direction), point of origin is the vertical center of the display (down - negative, up - positive)';'stable vertical position of viewers eye during movement pattern (in horizontal direction), point of origin is the vertical center of the display (down - negative, up - positive)';'angular horizontal start point of eye movement, point of origin is the horizontal center of the display (left - negative, right - positive)';'horizontal start point of eye movement, point of origin is the horizontal center of the display (left - negative, right - positive)';'angular horizontal stop point of eye movement, point of origin is the horizontal center of the display (left - negative, right - positive)';'horizontal stop point of eye movement, point of origin is the horizontal center of the display (left - negative, right - positive)';'angular distance of eye movement in horizontal direction';'distance of eye movement in horizontal direction';'angular velocity of horizontal eye movement, display position dependent parameter, maximum value during sequence (positive values - left to right movement, negative values right to left movement)';'velocity of horizontal eye movement (positive values - left to right movement, negative values right to left movement)';'distance between observer (corneal vertex) and level of displaying unit (e.g. projector)';'diagonal viewing angle of the displaying unit (e.g. projector) from observer position';'horizontal viewing angle of the displaying unit (e.g. projector) from observer position';'vertical viewing angle of the displaying unit (e.g. projector) from observer position';'corneal vertex (CV) to first nodal point (N) - relaxed accommodation';'corneal vertex (CV) to second nodal point (N´) - relaxed accommodation';'corneal vertex (CV) to second principal focus (F´) – relaxed accommodation';'corneal vertex (CV) to first nodal point (N) - accommodated (10.9D)';'corneal vertex (CV) to second nodal point (N´) - accommodated (10.9D)';'corneal vertex (CV) to second principal focus (F´) – accommodated (10.9D)';'cornela vertex (CV) to eyes center of rotation (Z)';'radius of retina (RR), simplified assumption';'threshold time for temporal summation referring to cones (inactive, set to zero)';'threshold time for temporal summation referring to rods (inactive, set to zero)';'threshold space for spatial summation (inactive, set to zero)';'angular eye movement velocity threshold in [DEG/SEC] for transition from PH1 to PH2 (values for threshold calculation refer to the frame cycle with the highest eye movement velocity in [DEG/SEC], variable values (specifically content width) because of the flat screen characteristic of the display unit)';'eye movement velocity threshold in [PX/FR] for transition from PH1 to PH2';'angular eye movement velocity threshold in [DEG/SEC] for transition from PH2 to PH3 (values for threshold calculation refer to the frame cycle with the highest eye movement velocity in [DEG/SEC], variable values (specifically content width) because of the flat screen characteristic of the display unit)';'eye movement velocity threshold in [PX/FR] for transition from PH2 to PH3';'final cbu score in [DEG2] for option 1 (area of color biased retinal stimulation for defined conditions, no light intensity considered, no specification of light color considered)';'final non-cbu score in [DEG2] for option 1 (area of unbiased retinal stimulation for defined conditions, no light intensity considered, no specification of light color considered)';'final reference score in [DEG2] for option 1 (retinal area that would be stimulated with unbiased content if no relative retinal movement was existent, no light intensity considered, no specification of light color considered)';'final cbu score for option 2 (volume of color biased retinal stimulation for defined conditions, area of retinal stimulation expressed in [DEG2], no specification of light color considered)';'final non-cbu score for option 2 (volume of unbiased retinal stimulation for defined conditions, area of retinal stimulation expressed in [DEG2], no specification of light color considered)';'final reference score for option 2 (volume of retinal stimulation that corresponds to unbiased content stimulation if no relative retinal movement was existent, area of retinal stimulation expressed in [DEG2], no specification of light color considered)'};

overview_table = table(CATEGORY, PARAMETER, VALUE, UNIT, ANNOTATION)

%% SEGMENT XIII - GRAPHICAL ILLUSTRATION 

close all;

% set XYZ tristimulus values of single primary colors (SFs) and resulting mixed color (harmonized SF mixture) to RGB triplet values (sRGB, white point D65) in order to illustrate colors of primaries and mixed color correctly in following figures
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
sf_XYZ = [sf_Xtristimulus sf_Ytristimulus sf_Ztristimulus]; % merging X,Y and Z tristimulus values of subframes (three arrays) to one matrix

sf_RGB = xyz2rgb(sf_XYZ); % convert subframe XYZ tristimulus values to subframe RGB triplet (sRGB, white point D65)

for bb = 1:1:size(sf_RGB,1) % normalize subframe RGB triplet values to values ranging from 0 to 1
   
    if sf_RGB(bb,1) < 0 || sf_RGB(bb,2) < 0 || sf_RGB(bb,3) < 0
        sf_RGB(bb,1:3) = sf_RGB(bb,1:3) - min(sf_RGB(bb,1:3));
        m = min(sf_RGB(bb,1:3)); range = max(sf_RGB(bb,1:3)) - m; sf_RGB(bb,1:3) = (sf_RGB(bb,1:3) - m) / range;
    end
    
    if sf_RGB(bb,1) > 1 || sf_RGB(bb,2) > 1 || sf_RGB(bb,3) > 1
        sf_RGB(bb,1:3) = sf_RGB(bb,1:3)/max(sf_RGB(bb,1:3));
    end
    
end

mix_XYZ = [mix_Xtristimulus mix_Ytristimulus mix_Ztristimulus]; % merging X,Y and Z tristimulus values of mixed color to one array

mix_RGB = (xyz2rgb(mix_XYZ)); % convert mixed color XYZ tristimulus values to mixed color RGB triplet (sRGB, white point D65)

if mix_RGB(1,1) < 0 || mix_RGB(1,2) < 0 || mix_RGB(1,3) < 0 % normalize mixed color RGB triplet values to values ranging from 0 to 1
    mix_RGB(1,1:3) = mix_RGB(1,1:3) - min(mix_RGB(1,1:3));
    m = min(mix_RGB(1,1:3)); range = max(mix_RGB(1,1:3)) - m; mix_RGB(1,1:3) = (mix_RGB(1,1:3) - m) / range;
end

if mix_RGB(1,1) > 1 || mix_RGB(1,2) > 1 || mix_RGB(1,3) > 1
    mix_RGB(1,1:3) = mix_RGB(1,1:3)/max(mix_RGB(1,1:3));
end
% -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


f1 = figure (1); hold on; set(gca,'FontSize',11); % 2D SPACE-TIME PLOT (x-axis: retinal position in [DEG] // y-axis: presentation time in [SEC])

for aa = 1:1:frame_number % loop for plotting of basic rhomboids 
    
    if aa ~= eye_velo_hor_variable_degs_table_index % all rhomboids corresponding to all frame cycles within sequence
        
        for bb = 0:1:subframe_number-1
            
            x = subframe_location_table_hor_x1(1+(bb*4):1:4+(bb*4),aa);
            x(5,1) = x(1,1);
            
            y = subframe_time_table(1+(bb*4):1:4+(bb*4),aa);
            y(5,1) = y(1,1);
            
            p((aa-1).*subframe_number+bb+1) = plot(x, y, ':o','LineWidth', 0.5, 'Color', [0.75, 0.75, 0.75],'MarkerEdgeColor','k','MarkerFaceColor',[0.75, 0.75, 0.75],'MarkerSize',3.0);
            
        end
        
    elseif aa == eye_velo_hor_variable_degs_table_index % specific rhomboid corresponding to frame cycle with highest cbu score within sequence
        
        for bb = 0:1:subframe_number-1
            
            x = subframe_location_table_hor_x1(1+(bb*4):1:4+(bb*4),aa);
            x(5,1) = x(1,1);
            
            y = subframe_time_table(1+(bb*4):1:4+(bb*4),aa);
            y(5,1) = y(1,1);
            
            p((aa-1).*subframe_number+bb+1) = plot(x, y, ':o','LineWidth', 0.5,'Color',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerSize',4.0);
            
        end
        
        for cc = 1:1:size(cbu_table_x1,1) % loop for plotting of color sector threshold lines (only for rhomboid that corresponds to frame cycle with highest cbu score)
            
            p(frame_number.*subframe_number+cc) = line([cbu_table_x1(cc,aa) cbu_table_x1(cc,aa)],[0 max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10], 'Color', [0.6, 0.6, 0.6],'LineStyle','-.','LineWidth', 0.5);
            
        end
             
    end
     
end

xlim([min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10 max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10]);
set(gca,'xTick',min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10:((max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10)-(min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10))/5:max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10);
xlabel('HORIZONTAL POSITION [DEG]','Fontsize',14)
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));
set(gca,'TickDir','out');
xtickangle(-90);

ylim([0 max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10]);
set(gca,'yTick',0:((max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10)-(0))/5:max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10);
ylabel('PRESENTATION TIME [SEC]','Fontsize',14)
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));
ytickangle(0);

if subframe_number == 2
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+2) (frame_number.*subframe_number+1)]),['2D subframe rhomboid (all SF, all FC)'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 3
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+2) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+3) (frame_number.*subframe_number+1)]),['2D subframe rhomboid (all SF, all FC)'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 4
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+2) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+3) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+4) (frame_number.*subframe_number+1)]),['2D subframe rhomboid (all SF, all FC)'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 5
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+2) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+3) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+4) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+5) (frame_number.*subframe_number+1)]),['2D subframe rhomboid (all SF, all FC)'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 6
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+2) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+3) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+4) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+5) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number+6) (frame_number.*subframe_number+1)]),['2D subframe rhomboid (all SF, all FC)'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D subframe rhomboid (SF' num2str(subframe_number-(subframe_number-6)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5,'Location','northeastoutside');
    
end
    
grid on
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];

title ('SPACE-TIME PLOT', 'Fontsize', 16)
set(gcf,'WindowStyle','docked');

hold off;


%{
f8 = figure (8); f8 = plot3(0,0,0); hold on; set(gca,'FontSize',11); % 3D SPACE-TIME-INTENSITY PLOT (x-axis: retinal position in [DEG] // y-axis: presentation time in [SEC] // z-axis: radiance in [W/(SR.*M2)])

for aa = 1:1:frame_number % loop for plotting of 3D subframe objects
    
    if aa ~= eye_velo_hor_variable_degs_table_index % all subframe objects corresponding to all frame cycles within sequence
              
        for bb = 0:1:subframe_number-1
            
            x = subframe_location_table_hor_x1(1+(bb*4):1:4+(bb*4),aa);
            x(5,1) = x(1,1);
            
            y = subframe_time_table(1+(bb*4):1:4+(bb*4),aa);
            y(5,1) = y(1,1);
            
            z0 = zeros(5,1);
            
            z = subframe_radiance_intensity_table(1+(bb*4):1:4+(bb*4),aa);
            z(5,1) = z(1,1);
            
            plot3(x, y, z0, ':o','LineWidth', 0.5,'Color', [0.75, 0.75, 0.75],'MarkerEdgeColor','k','MarkerFaceColor',[0.75, 0.75, 0.75],'MarkerSize',3.0);
            plot3(x, y, z, ':o','LineWidth', 0.5,'Color', [0.75, 0.75, 0.75],'MarkerEdgeColor','k','MarkerFaceColor',[0.75, 0.75, 0.75],'MarkerSize',3.0);
            
            plot3([x(1,1) x(1,1)],[y(1,1) y(1,1)], [z0(1,1) z(1,1)],':','LineWidth',0.5,'Color', [0.75, 0.75, 0.75]);
            plot3([x(2,1) x(2,1)],[y(2,1) y(2,1)], [z0(2,1) z(2,1)],':','LineWidth', 0.5,'Color', [0.75, 0.75, 0.75]);
            plot3([x(3,1) x(3,1)],[y(3,1) y(3,1)], [z0(3,1) z(3,1)],':','LineWidth', 0.5,'Color', [0.75, 0.75, 0.75]);
            plot3([x(4,1) x(4,1)],[y(4,1) y(4,1)], [z0(4,1) z(4,1)],':','LineWidth', 0.5,'Color', [0.75, 0.75, 0.75]);
            
        end
        
    elseif aa == eye_velo_hor_variable_degs_table_index % specific subframe object corresponding to frame cycle with highest cbu score within sequence
        
        for bb = 0:1:subframe_number-1
            
            x = subframe_location_table_hor_x1(1+(bb*4):1:4+(bb*4),aa);
            x(5,1) = x(1,1);
            
            y = subframe_time_table(1+(bb*4):1:4+(bb*4),aa);
            y(5,1) = y(1,1);
            
            z0 = zeros(5,1);
            
            z = subframe_radiance_intensity_table(1+(bb*4):1:4+(bb*4),aa);
            z(5,1) = z(1,1);
            
            plot3(x, y, z0, ':ko','LineWidth', 0.5,'MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerSize',4.0);
            plot3(x, y, z, ':ko','LineWidth', 0.5,'MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerSize',4.0);
            
            plot3([x(1,1) x(1,1)],[y(1,1) y(1,1)], [z0(1,1) z(1,1)], ':k','LineWidth', 0.5);
            plot3([x(2,1) x(2,1)],[y(2,1) y(2,1)], [z0(2,1) z(2,1)], ':k','LineWidth', 0.5);
            plot3([x(3,1) x(3,1)],[y(3,1) y(3,1)], [z0(3,1) z(3,1)], ':k','LineWidth', 0.5);
            plot3([x(4,1) x(4,1)],[y(4,1) y(4,1)], [z0(4,1) z(4,1)], ':k','LineWidth', 0.5);
            
        end
        
        for cc = 1:1:size(cbu_table_x1,1) % loop for plotting of color sector threshold lines (only for subframe object that corresponds to frame cycle with highest cbu score)
            
            plot3([cbu_table_x1(cc,aa) cbu_table_x1(cc,aa)],[0 max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10],[0,0], 'Color', [0.6, 0.6, 0.6],'LineStyle','-.','LineWidth', 0.5);
            
        end
        
    end
    
end

xlim([min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10 max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10]);
set(gca,'xTick',min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10:((max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10)-(min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10))/5:max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));

ylim([0 max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10]);
set(gca,'yTick',0:((max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10)-(0))/5:max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10);  
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));

zlim([0 max(subframe_radiance_intensity_table(:))+(max(subframe_radiance_intensity_table(:))-min(subframe_radiance_intensity_table(:)))/10]);
set(gca,'zTick', 0:((max(subframe_radiance_intensity_table(:))+(max(subframe_radiance_intensity_table(:))-min(subframe_radiance_intensity_table(:)))/10)-(0))/5:max(subframe_radiance_intensity_table(:))+(max(subframe_radiance_intensity_table(:))-min(subframe_radiance_intensity_table(:)))/10);
set(gca, 'ZTicklabel', num2str(get(gca, 'ZTick')', '%.2f'));

title ('SPACE-TIME-INTENSITY PLOT (1)', 'Fontsize', 16);
text('Units', 'normalized', 'Position', [0.0 -0.11], 'String', '>> X-AXIS: HORIZONTAL POSITION [DEG] // Y-AXIS: PRESENTATION TIME [SEC] // Z-AXIS: RADIANCE [W/(SR.*M2)] <<','Fontsize',14);
set(gcf,'WindowStyle','docked');

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];

view (5,40); 

hold off;
%}


f2 = figure (2); f2 = plot3(0,0,0); hold on; set(gca,'FontSize',11); % 3D SPACE-TIME-INTENSITY PLOT (x-axis: retinal position in [DEG] // y-axis: presentation time in [SEC] // z-axis: luminance in [CD/M2])

for aa = 1:1:frame_number % loop for plotting of 3D subframe objects
    
    if aa ~= eye_velo_hor_variable_degs_table_index % all subframe objects corresponding to all frame cycles within sequence
        
        for bb = 0:1:subframe_number-1
            
            x = subframe_location_table_hor_x1(1+(bb*4):1:4+(bb*4),aa);
            x(5,1) = x(1,1);
            
            y = subframe_time_table(1+(bb*4):1:4+(bb*4),aa);
            y(5,1) = y(1,1);
            
            z0 = zeros(5,1);
            
            z = subframe_luminance_intensity_table(1+(bb*4):1:4+(bb*4),aa);
            z(5,1) = z(1,1);
            
            p((aa-1).*subframe_number.*6+(bb.*6)+1) = plot3(x, y, z0, ':o','LineWidth', 0.5,'Color', [0.75, 0.75, 0.75],'MarkerEdgeColor','k','MarkerFaceColor',[0.75, 0.75, 0.75],'MarkerSize',3.0);
            p((aa-1).*subframe_number.*6+(bb.*6)+2) = plot3(x, y, z, ':o','LineWidth', 0.5,'Color', [0.75, 0.75, 0.75],'MarkerEdgeColor','k','MarkerFaceColor',[0.75, 0.75, 0.75],'MarkerSize',3.0);
            
            p((aa-1).*subframe_number.*6+(bb.*6)+3) = plot3([x(1,1) x(1,1)],[y(1,1) y(1,1)], [z0(1,1) z(1,1)], ':','LineWidth',0.5,'Color', [0.75, 0.75, 0.75]);
            p((aa-1).*subframe_number.*6+(bb.*6)+4) = plot3([x(2,1) x(2,1)],[y(2,1) y(2,1)], [z0(2,1) z(2,1)], ':','LineWidth',0.5,'Color', [0.75, 0.75, 0.75]);
            p((aa-1).*subframe_number.*6+(bb.*6)+5) = plot3([x(3,1) x(3,1)],[y(3,1) y(3,1)], [z0(3,1) z(3,1)], ':','LineWidth',0.5,'Color', [0.75, 0.75, 0.75]);
            p((aa-1).*subframe_number.*6+(bb.*6)+6) = plot3([x(4,1) x(4,1)],[y(4,1) y(4,1)], [z0(4,1) z(4,1)], ':','LineWidth',0.5,'Color', [0.75, 0.75, 0.75]);
            
        end
        
    elseif aa == eye_velo_hor_variable_degs_table_index % specific subframe object corresponding to frame cycle with highest cbu score within sequence
        
        for bb = 0:1:subframe_number-1
            
            x = subframe_location_table_hor_x1(1+(bb*4):1:4+(bb*4),aa);
            x(5,1) = x(1,1);
            
            y = subframe_time_table(1+(bb*4):1:4+(bb*4),aa);
            y(5,1) = y(1,1);
            
            z0 = zeros(5,1);
            
            z = subframe_luminance_intensity_table(1+(bb*4):1:4+(bb*4),aa);
            z(5,1) = z(1,1);
            
            p((aa-1).*subframe_number.*6+(bb.*6)+1) = plot3(x, y, z0, ':o','LineWidth', 0.5,'Color',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerSize',4.0);
            p((aa-1).*subframe_number.*6+(bb.*6)+2) = plot3(x, y, z, ':o','LineWidth', 0.5,'Color',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerSize',4.0);
            
            p((aa-1).*subframe_number.*6+(bb.*6)+3) = plot3([x(1,1) x(1,1)],[y(1,1) y(1,1)], [z0(1,1) z(1,1)], ':','LineWidth', 0.5,'Color',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
            p((aa-1).*subframe_number.*6+(bb.*6)+4) = plot3([x(2,1) x(2,1)],[y(2,1) y(2,1)], [z0(2,1) z(2,1)], ':','LineWidth', 0.5,'Color',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
            p((aa-1).*subframe_number.*6+(bb.*6)+5) = plot3([x(3,1) x(3,1)],[y(3,1) y(3,1)], [z0(3,1) z(3,1)], ':','LineWidth', 0.5,'Color',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
            p((aa-1).*subframe_number.*6+(bb.*6)+6) = plot3([x(4,1) x(4,1)],[y(4,1) y(4,1)], [z0(4,1) z(4,1)], ':','LineWidth', 0.5,'Color',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
            
        end
        
        for cc = 1:1:size(cbu_table_x1,1) % loop for plotting of color sector threshold lines (only for subframe object that corresponds to frame cycle with highest cbu score)
            
            p(frame_number.*subframe_number.*6+cc) = plot3([cbu_table_x1(cc,aa) cbu_table_x1(cc,aa)],[0 max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10],[0,0], 'Color', [0.6, 0.6, 0.6],'LineStyle','-.','LineWidth', 0.5);
            
        end
        
    end
    
end

xlim([min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10 max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10]);
set(gca,'xTick',min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10:((max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10)-(min(subframe_location_table_hor_x1(:))-(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10))/5:max(subframe_location_table_hor_x1(:))+(max(subframe_location_table_hor_x1(:))-min(subframe_location_table_hor_x1(:)))/10);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));
xtickangle(-108);

ylim([0 max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10]);
set(gca,'yTick',0:((max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10)-(0))/5:max(subframe_time_table(:))+(max(subframe_time_table(:))-min(subframe_time_table(:)))/10);  
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));
ytickangle(-2);

zlim([0 max(subframe_luminance_intensity_table(:))+(max(subframe_luminance_intensity_table(:))-min(subframe_luminance_intensity_table(:)))/10]);
set(gca,'zTick', 0:((max(subframe_luminance_intensity_table(:))+(max(subframe_luminance_intensity_table(:))-min(subframe_luminance_intensity_table(:)))/10)-(0))/5:max(subframe_luminance_intensity_table(:))+(max(subframe_luminance_intensity_table(:))-min(subframe_luminance_intensity_table(:)))/10);
set(gca, 'ZTicklabel', num2str(get(gca, 'ZTick')', '%.2f'));
ztickangle(-2);

title ('SPACE-TIME-INTENSITY PLOT (2)', 'Fontsize', 16)
text('Units', 'normalized', 'Position', [-0.05 -0.11], 'String', '>> X-AXIS: HORIZONTAL POSITION [DEG] // Y-AXIS: PRESENTATION TIME [SEC] // Z-AXIS: LUMINANCE [CD/M2] <<','Fontsize',14);

if subframe_number == 2
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+7) (frame_number.*subframe_number.*6+1)]),['3D subframe body (all SF, all FC)'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D subframe body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5);
    
elseif subframe_number == 3
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+7) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+13) (frame_number.*subframe_number.*6+1)]),['3D subframe body (all SF, all FC)'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D subframe body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5);
    
elseif subframe_number == 4
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+7) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+13) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+19) (frame_number.*subframe_number.*6+1)]),['3D subframe body (all SF, all FC)'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D subframe body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5);
    
elseif subframe_number == 5
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+7) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+13) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+19) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+25) (frame_number.*subframe_number.*6+1)]),['3D subframe body (all SF, all FC)'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D subframe body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5);
    
elseif subframe_number == 6
    
    legend (p([1 ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+1) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+7) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+13) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+19) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+25) ((eye_velo_hor_variable_degs_table_index-1).*subframe_number.*6+31) (frame_number.*subframe_number.*6+1)]),['3D subframe body (all SF, all FC)'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D subframe body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D subframe body (SF' num2str(subframe_number-(subframe_number-6)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],'Fontsize', 7.5);
    
end

set(gcf,'WindowStyle','docked');

view (5,40); 

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];

hold off;


f3 = figure (3); hold on; set(gca,'FontSize',11); % INTENSITY PROFILE - 2D plot refering to the perceived intensity of stimulation after temporal summation (intensity over time) of all subframes within the frame cycle of choice (x-axis: horizontal retinal position in [DEG] // y-axis: perceived intensity)

aa = eye_velo_hor_variable_degs_table_index; % loop for plotting of every single subframe's perceived intensity after temporal summation (refering to the chosen frame cycle with highest eye movement velocity in [DEG/SEC] and therefore the largest cbu effect)

for bb = 0:1:subframe_number-1
    
    x = subframe_location_table_hor_sorted1(1+(bb*4):1:4+(bb*4),aa);
    x(5,1) = x(1,1);
    
    y = subframe_integrated_luminance_table(1+(bb*4):1:4+(bb*4),aa);
    y(5,1) = y(1,1);
    
    p(bb+1) = plot(x, y, ':o','LineWidth', 0.5, 'Color', [sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerSize',4.0);
    
end

p(subframe_number+1) = plot(xyL_table(1:end,1), xyL_table(1:end,end), '-', 'Color', [0.85 0.3 0.3], 'LineWidth', 2.0); % plotting the line for the resulting perceived intensity (after temporal summation of single subframes) a.k.a. intensity profile

for cc = 1:1:size(cbu_table_x1,1) % loop for plotting of color sector threshold lines
    
    p(subframe_number+1+cc) = line([cbu_table_x1(cc,aa) cbu_table_x1(cc,aa)],[0 max(xyL_table(1:end,end))+(max(xyL_table(1:end,end))-min(xyL_table(1:end,end)))/10], 'Color', [0.6, 0.6, 0.6],'LineStyle','-.','LineWidth', 0.5);
    
end

location_table_short_sorted = location_table_short(:);

xlim([min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10 max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))+(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10]);
set(gca,'xTick',min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10:((max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))+(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10)-(min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10))/5:max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))+(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10);
xlabel('HORIZONTAL POSITION [DEG]','Fontsize',14)
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));
set(gca,'TickDir','out');
xtickangle(-90);

ylim([0 max(xyL_table(1:end,end))+(max(xyL_table(1:end,end))-min(xyL_table(1:end,end)))/10]);
set(gca,'yTick',0:((max(xyL_table(1:end,end))+(max(xyL_table(1:end,end))-min(xyL_table(1:end,end)))/10)-(0))/5:max(xyL_table(1:end,end))+(max(xyL_table(1:end,end))-min(xyL_table(1:end,end)))/10);
ylabel('PERCEIVED INTENSITY','Fontsize',14)
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));
ytickangle(0);

if subframe_number == 2
    
    legend ({['2D intensity profile (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')']},'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 3
    
    legend ({['2D intensity profile (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')']},'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 4
    
    legend ({['2D intensity profile (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')']},'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 5
    
    legend ({['2D intensity profile (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')']},'Fontsize', 7.5,'Location','northeastoutside');
    
elseif subframe_number == 6
    
    legend ({['2D intensity profile (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (SF' num2str(subframe_number-(subframe_number-6)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['2D intensity profile (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')']},'Fontsize', 7.5,'Location','northeastoutside');
    
end

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];

title ('INTENSITY PROFILE', 'Fontsize', 16)
set(gcf,'WindowStyle','docked');

hold off;


f4 = figure(4); f4 = plot3(0,0,0); hold on; set(gca,'FontSize',11); % INTENSITY MAP - 3D plot refering to the perceived intensity of stimulation after temporal summation (intensity over time) of all subframes within the frame cycle of choice (x-axis: horizontal retinal position in [DEG] // y-axis: vertical retinal position in [DEG] // z-axis: perceived intensity)

aa = eye_velo_hor_variable_degs_table_index; % loop for plotting of every single subframe's perceived intensity after temporal summation (refering to the chosen frame cycle with highest eye movement velocity in [DEG/SEC] and therefore the largest cbu effect)
 
for bb = 0:1:subframe_number-1
    
    for cc = 1:1:size(subframe_location_table_ver_x1,1)
        
        x = subframe_location_table_hor_sorted1(1+(bb*4):1:4+(bb*4),aa);
        x(5,1) = x(1,1);
        
        y = subframe_location_table_ver_x1(cc,1).*ones(5,1);
        
        z = subframe_integrated_luminance_table(1+(bb*4):1:4+(bb*4),aa);
        z(5,1) = z(1,1);
        
        p(bb.*6+cc) = plot3(x, y, z, ':o','LineWidth', 0.5, 'Color', [sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)], 'MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)],'MarkerSize',4.0);
        
    end
    
    p(bb.*6+cc+1) = plot3([x(1,1) x(1,1)],[subframe_location_table_ver_x1(1,1) subframe_location_table_ver_x1(2,1)], [z(1,1) z(1,1)], ':','LineWidth', 0.5, 'Color', [sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
    p(bb.*6+cc+2) = plot3([x(2,1) x(2,1)],[subframe_location_table_ver_x1(1,1) subframe_location_table_ver_x1(2,1)], [z(2,1) z(2,1)], ':','LineWidth', 0.5, 'Color', [sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
    p(bb.*6+cc+3) = plot3([x(3,1) x(3,1)],[subframe_location_table_ver_x1(1,1) subframe_location_table_ver_x1(2,1)], [z(3,1) z(3,1)], ':','LineWidth', 0.5, 'Color', [sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
    p(bb.*6+cc+4) = plot3([x(4,1) x(4,1)],[subframe_location_table_ver_x1(1,1) subframe_location_table_ver_x1(2,1)], [z(4,1) z(4,1)], ':','LineWidth', 0.5, 'Color', [sf_RGB(bb+1,1),sf_RGB(bb+1,2),sf_RGB(bb+1,3)]);
    
end

aa = eye_velo_hor_variable_degs_table_index; % loop for 3D plotting of resulting perceived intensity after temporal summation of all subframes of the chosen frame cycle with highest eye movement velocity in [DEG/SEC] and therefore the largest cbu effect
 
for dd = 1:1:size(subframe_location_table_ver_x1,1)
    
    x = xyL_table(1:end,1);
    x(end+1,1) = x(1,1);
    
    y = subframe_location_table_ver_x1(dd,1).*ones(size(xyL_table,1)+1,1);
    
    z = xyL_table(1:end,end);
    z(end+1,1) = z(1,1);
    
    p(bb.*6+cc+4+dd) = plot3(x, y, z, '-', 'Color', [0.85 0.3 0.3], 'LineWidth', 2.0);
    
end

for ee = 1:1:(integrated_subframes.*4)
    
    ff = find(xyL_table == subframe_location_table_hor_sorted1(ee,eye_velo_hor_variable_degs_table_index));
    
    p(bb.*6+cc+4+dd+ee) = plot3([subframe_location_table_hor_sorted1(ee,eye_velo_hor_variable_degs_table_index) subframe_location_table_hor_sorted1(ee,eye_velo_hor_variable_degs_table_index)],[subframe_location_table_ver_x1(1,1) subframe_location_table_ver_x1(2,1)], [xyL_table(ff,end) xyL_table(ff,end)], '-', 'Color', [0.85 0.3 0.3], 'LineWidth', 1.5);
    
end

p(bb.*6+cc+4+dd+ee+1) = plot3([0 0],[min(subframe_location_table_ver_x1)-(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10 max(subframe_location_table_ver_x1)+(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10], [0 0],'LineStyle','-.', 'Color', [0.65, 0.65, 0.65],'LineWidth', 1.0); % zero lines of retinal positioning for horizontal and vertical direction (intersection point represents fovea centralis)
p(bb.*6+cc+4+dd+ee+2) = plot3([min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10 max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))+(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10], [0 0], [0 0] , 'LineStyle','-.', 'Color', [0.65, 0.65, 0.65],'LineWidth', 1.0);

for r1 = 0.25:0.25:2 % loop for plotting circle grid (0.25DEG grid)
    
n = 1000;
t = linspace(0,2*pi,n);
x = r1*sin(t);
y = r1*cos(t);
z = zeros(1,size(x,2));
p(bb.*6+cc+4+dd+ee+2+r1/0.25) = line(x,y,z,'LineStyle','--','Color', [0.85, 0.85, 0.85],'LineWidth', 0.5);

end

for r2 = 0.5:0.5:2 % loop for plotting circle grid (0.5DEG grid)
    
n = 1000;
t = linspace(0,2*pi,n);
x = r2*sin(t);
y = r2*cos(t);
z = zeros(1,size(x,2));
p(bb.*6+cc+4+dd+ee+2+r1/0.25+r2/0.5) = line(x,y,z, 'LineStyle','--', 'Color', [0.75, 0.75, 0.75],'LineWidth', 0.5);

end

xlim([min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10 max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))+(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10]);
set(gca,'xTick',min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10:((max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))+(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10)-(min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10))/5:max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))+(max((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-min((location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))/10);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));
% xlabel('HORIZONTAL POSITION [DEG]', 'Position', [0.65, 0.0, -0.09], 'Rotation', 5, 'Fontsize',14);
xtickangle(-45);

ylim([min(subframe_location_table_ver_x1)-(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10 max(subframe_location_table_ver_x1)+(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10]);
set(gca,'yTick',min(subframe_location_table_ver_x1)-(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10:((max(subframe_location_table_ver_x1)+(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10)-(min(subframe_location_table_ver_x1)-(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10))/5:max(subframe_location_table_ver_x1)+(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));
% ylabel('VERTICAL POSITION [DEG]', 'Position', [-1.93, -0.085, 0.0], 'Rotation', -45, 'Fontsize',14);
ytickangle(+10);

zlim([0 max(xyL_table(1:end,end))+(max(xyL_table(1:end,end))-min(xyL_table(1:end,end)))/10]);
set(gca,'zTick',0:((max(xyL_table(1:end,end))+(max(xyL_table(1:end,end))-min(xyL_table(1:end,end)))/10)-(0))/5:max(xyL_table(1:end,end))+(max(xyL_table(1:end,end))-min(xyL_table(1:end,end)))/10);
set(gca, 'ZTicklabel', num2str(get(gca, 'ZTick')', '%.2f'));
% zlabel('PERCEIVED INTENSITY', 'Position', [-2.24, 0.0, 0.17], 'Fontsize',14);
ztickangle(-45);

text('Units', 'normalized', 'Position', [-0.03 -0.11], 'String', '>> X-AXIS: HORIZONTAL POSITION [DEG] // Y-AXIS: VERTICAL POSITION [DEG] // Z-AXIS: PERCEIVED INTENSITY <<','Fontsize',14);

title ('INTENSITY MAP', 'Fontsize', 16)

if subframe_number == 2
    
    legend (p([1 ((subframe_number-1).*6+1) (subframe_number.*6+1) (subframe_number.*6+3+subframe_number.*4) (subframe_number.*6+3+subframe_number.*4+2) (subframe_number.*6+3+subframe_number.*4+2+r1/0.25)]),['3D intensity body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D intensity body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['retinal zero lines'],['retinal 0.25deg circle grid'],['retinal 0.5deg circle grid'],'Fontsize', 7.5);
    
elseif subframe_number == 3
    
    legend (p([1 ((subframe_number-2).*6+1) ((subframe_number-1).*6+1) (subframe_number.*6+1) (subframe_number.*6+3+subframe_number.*4) (subframe_number.*6+3+subframe_number.*4+2) (subframe_number.*6+3+subframe_number.*4+2+r1/0.25)]),['3D intensity body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D intensity body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['retinal zero lines'],['retinal 0.25deg circle grid'],['retinal 0.5deg circle grid'],'Fontsize', 7.5);
    
elseif subframe_number == 4
    
    legend (p([1 ((subframe_number-3).*6+1) ((subframe_number-2).*6+1) ((subframe_number-1).*6+1) (subframe_number.*6+1) (subframe_number.*6+3+subframe_number.*4) (subframe_number.*6+3+subframe_number.*4+2) (subframe_number.*6+3+subframe_number.*4+2+r1/0.25)]),['3D intensity body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D intensity body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['retinal zero lines'],['retinal 0.25deg circle grid'],['retinal 0.5deg circle grid'],'Fontsize', 7.5);
    
elseif subframe_number == 5
    
    legend (p([1 ((subframe_number-4).*6+1) ((subframe_number-3).*6+1) ((subframe_number-2).*6+1) ((subframe_number-1).*6+1) (subframe_number.*6+1) (subframe_number.*6+3+subframe_number.*4) (subframe_number.*6+3+subframe_number.*4+2) (subframe_number.*6+3+subframe_number.*4+2+r1/0.25)]),['3D intensity body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D intensity body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['retinal zero lines'],['retinal 0.25deg circle grid'],['retinal 0.5deg circle grid'],'Fontsize', 7.5);
    
elseif subframe_number == 6
    
    legend (p([1 ((subframe_number-5).*6+1) ((subframe_number-4).*6+1) ((subframe_number-3).*6+1) ((subframe_number-2).*6+1) ((subframe_number-1).*6+1) (subframe_number.*6+1) (subframe_number.*6+3+subframe_number.*4) (subframe_number.*6+3+subframe_number.*4+2) (subframe_number.*6+3+subframe_number.*4+2+r1/0.25)]),['3D intensity body (SF' num2str(subframe_number-(subframe_number-1)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'], ['3D intensity body (SF' num2str(subframe_number-(subframe_number-2)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-3)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-4)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-5)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (SF' num2str(subframe_number-(subframe_number-6)) ', FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['3D intensity body (all SF, FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['retinal zero lines'],['retinal 0.25deg circle grid'],['retinal 0.5deg circle grid'],'Fontsize', 7.5);
    
end

view (-19,17); 

h = get(gca,'DataAspectRatio'); % equalize x- and y-axis by leaving z-axis out

if h(3)==1 
      set(gca,'DataAspectRatio',[1 1 2/max(h(1:2))]);
else  
    set(gca,'DataAspectRatio',[1 1 2*h(3)]);
end

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];

set(gcf,'WindowStyle','docked');
hold off;


f5 = figure (5); hold on; set(gca,'FontSize',11); % COLOR PROFILE - plotting of CIE xy chromaticity diagram with color gamut from applied subframes and horizontal CBU color profile

% CIE standard observer color matching function (CIE Lab 1931, 2.0deg field of view, delta_lambda = 5.0nm) 
% column 1 = wavelength (400-700nm), column 2 = x(lambda), column 3 = y(lambda), column 4 = z(lambda)
% absolute values in cie_1931_absolute, relative values in cie_1931_relative >> x(lambda) + y(lambda) + z(lambda) = 1
% from Shaw, M., & Fairchild, M. (2002). Evaluating the 1931 CIE color?matching functions. Color Research & Application, 27(5), 316-329 >> APPENDIX (!) >> further color spaces available (rgb, lms)

cie_1931_absolute = [400 0.0143000000000000 0.000400000000000000 0.0679000000000000;405 0.0232000000000000 0.000600000000000000 0.110000000000000;410 0.0435000000000000 0.00120000000000000 0.207000000000000;415 0.0776000000000000 0.00220000000000000 0.371000000000000;420 0.134000000000000 0.00400000000000000 0.646000000000000;425 0.215000000000000 0.00730000000000000 1.04000000000000;430 0.284000000000000 0.0116000000000000 1.39000000000000;435 0.329000000000000 0.0168000000000000 1.62000000000000;440 0.348000000000000 0.0230000000000000 1.75000000000000;445 0.348000000000000 0.0298000000000000 1.78000000000000;450 0.336000000000000 0.0380000000000000 1.77000000000000;455 0.319000000000000 0.0480000000000000 1.74000000000000;460 0.291000000000000 0.0600000000000000 1.67000000000000;465 0.251000000000000 0.0739000000000000 1.53000000000000;470 0.195000000000000 0.0910000000000000 1.29000000000000;475 0.142000000000000 0.113000000000000 1.04000000000000;480 0.0956000000000000 0.139000000000000 0.813000000000000;485 0.0580000000000000 0.169000000000000 0.616000000000000;490 0.0320000000000000 0.208000000000000 0.465000000000000;495 0.0147000000000000 0.259000000000000 0.353000000000000;500 0.00490000000000000 0.323000000000000 0.272000000000000;505 0.00240000000000000 0.407000000000000 0.212000000000000;510 0.00930000000000000 0.503000000000000 0.158000000000000;515 0.0291000000000000 0.608000000000000 0.112000000000000;520 0.0633000000000000 0.710000000000000 0.0782000000000000;525 0.110000000000000 0.793000000000000 0.0573000000000000;530 0.166000000000000 0.862000000000000 0.0422000000000000;535 0.226000000000000 0.915000000000000 0.0298000000000000;540 0.290000000000000 0.954000000000000 0.0203000000000000;545 0.360000000000000 0.980000000000000 0.0134000000000000;550 0.433000000000000 0.995000000000000 0.00870000000000000;555 0.512000000000000 1 0.00570000000000000;560 0.595000000000000 0.995000000000000 0.00390000000000000;565 0.678000000000000 0.979000000000000 0.00270000000000000;570 0.762000000000000 0.952000000000000 0.00210000000000000;575 0.843000000000000 0.915000000000000 0.00180000000000000;580 0.916000000000000 0.870000000000000 0.00170000000000000;585 0.979000000000000 0.816000000000000 0.00140000000000000;590 1.03000000000000 0.757000000000000 0.00110000000000000;595 1.06000000000000 0.695000000000000 0.00100000000000000;600 1.06000000000000 0.631000000000000 0.000800000000000000;605 1.05000000000000 0.567000000000000 0.000600000000000000;610 1 0.503000000000000 0.000300000000000000;615 0.938000000000000 0.441000000000000 0.000200000000000000;620 0.854000000000000 0.381000000000000 0.000200000000000000;625 0.751000000000000 0.321000000000000 0.000100000000000000;630 0.642000000000000 0.265000000000000 0;635 0.542000000000000 0.217000000000000 0;640 0.448000000000000 0.175000000000000 0;645 0.361000000000000 0.138000000000000 0;650 0.284000000000000 0.107000000000000 0;655 0.219000000000000 0.0816000000000000 0;660 0.165000000000000 0.0610000000000000 0;665 0.121000000000000 0.0446000000000000 0;670 0.0874000000000000 0.0320000000000000 0;675 0.0636000000000000 0.0232000000000000 0;680 0.0468000000000000 0.0170000000000000 0;685 0.0329000000000000 0.0119000000000000 0;690 0.0227000000000000 0.00820000000000000 0;695 0.0158000000000000 0.00570000000000000 0;700 0.0114000000000000 0.00410000000000000 0;400 0.0143000000000000 0.000400000000000000 0.0679000000000000];    
cie_1931_relative = [400 0.173120000000000 0.00484000000000000 0.822030000000000;405 0.173390000000000 0.00448000000000000 0.822120000000000;410 0.172820000000000 0.00477000000000000 0.822410000000000;415 0.172140000000000 0.00488000000000000 0.822980000000000;420 0.170920000000000 0.00510000000000000 0.823980000000000;425 0.170320000000000 0.00578000000000000 0.823890000000000;430 0.168490000000000 0.00688000000000000 0.824630000000000;435 0.167360000000000 0.00855000000000000 0.824090000000000;440 0.164070000000000 0.0108400000000000 0.825080000000000;445 0.161280000000000 0.0138100000000000 0.824910000000000;450 0.156720000000000 0.0177200000000000 0.825560000000000;455 0.151400000000000 0.0227800000000000 0.825820000000000;460 0.143990000000000 0.0296900000000000 0.826320000000000;465 0.135320000000000 0.0398400000000000 0.824840000000000;470 0.123730000000000 0.0577400000000000 0.818530000000000;475 0.109650000000000 0.0872600000000000 0.803090000000000;480 0.0912600000000000 0.132680000000000 0.776060000000000;485 0.0688000000000000 0.200470000000000 0.730720000000000;490 0.0453900000000000 0.295040000000000 0.659570000000000;495 0.0234600000000000 0.413280000000000 0.563270000000000;500 0.00817000000000000 0.538420000000000 0.453410000000000;505 0.00386000000000000 0.654970000000000 0.341170000000000;510 0.0138700000000000 0.750410000000000 0.235720000000000;515 0.0388500000000000 0.811640000000000 0.149510000000000;520 0.0743400000000000 0.833820000000000 0.0918400000000000;525 0.114550000000000 0.825780000000000 0.0596700000000000;530 0.155110000000000 0.805460000000000 0.0394300000000000;535 0.193030000000000 0.781520000000000 0.0254500000000000;540 0.229380000000000 0.754570000000000 0.0160600000000000;545 0.266000000000000 0.724100000000000 0.00990000000000000;550 0.301390000000000 0.692560000000000 0.00606000000000000;555 0.337350000000000 0.658890000000000 0.00376000000000000;560 0.373300000000000 0.624250000000000 0.00245000000000000;565 0.408510000000000 0.589870000000000 0.00163000000000000;570 0.444030000000000 0.554750000000000 0.00122000000000000;575 0.479030000000000 0.519950000000000 0.00102000000000000;580 0.512390000000000 0.486660000000000 0.000950000000000000;585 0.544980000000000 0.454240000000000 0.000780000000000000;590 0.576030000000000 0.423350000000000 0.000620000000000000;595 0.603640000000000 0.395790000000000 0.000570000000000000;600 0.626550000000000 0.372980000000000 0.000470000000000000;605 0.649110000000000 0.350520000000000 0.000370000000000000;610 0.665200000000000 0.334600000000000 0.000200000000000000;615 0.680100000000000 0.319750000000000 0.000150000000000000;620 0.691390000000000 0.308450000000000 0.000160000000000000;625 0.700490000000000 0.299410000000000 9.00000000000000e-05;630 0.707830000000000 0.292170000000000 0;635 0.714100000000000 0.285900000000000 0;640 0.719100000000000 0.280900000000000 0;645 0.723450000000000 0.276550000000000 0;650 0.726340000000000 0.273660000000000 0;655 0.728540000000000 0.271460000000000 0;660 0.730090000000000 0.269910000000000 0;665 0.730680000000000 0.269320000000000 0;670 0.731990000000000 0.268010000000000 0;675 0.732720000000000 0.267280000000000 0;680 0.733540000000000 0.266460000000000 0;685 0.734380000000000 0.265630000000000 0;690 0.734630000000000 0.265370000000000 0;695 0.734880000000000 0.265120000000000 0;700 0.735480000000000 0.264520000000000 0;400 0.173120000000000 0.00484000000000000 0.822030000000000];

% CIE standard observer color matching function (CIE Lab 1964, 10.0deg field of view, delta_lambda = 5.0nm) 
% column 1 = wavelength (400-700nm), column 2 = x(lambda), column 3 = y(lambda), column 4 = z(lambda)
% absolute values in cie_1964_absolute, relative values in cie_1964_relative >> x(lambda) + y(lambda) + z(lambda) = 1
% from Shaw, M., & Fairchild, M. (2002). Evaluating the 1931 CIE color?matching functions. Color Research & Application, 27(5), 316-329 >> APPENDIX (!) >> further color spaces available (rgb, lms)

cie_1964_absolute = [400 0.0191000000000000 0.00200000000000000 0.0860000000000000;405 0.0434000000000000 0.00450000000000000 0.197000000000000;410 0.0847000000000000 0.00880000000000000 0.389000000000000;415 0.141000000000000 0.0145000000000000 0.657000000000000;420 0.205000000000000 0.0214000000000000 0.973000000000000;425 0.265000000000000 0.0295000000000000 1.28000000000000;430 0.315000000000000 0.0387000000000000 1.55000000000000;435 0.358000000000000 0.0496000000000000 1.80000000000000;440 0.384000000000000 0.0621000000000000 1.97000000000000;445 0.387000000000000 0.0747000000000000 2.03000000000000;450 0.371000000000000 0.0895000000000000 1.99000000000000;455 0.343000000000000 0.106000000000000 1.90000000000000;460 0.302000000000000 0.128000000000000 1.75000000000000;465 0.254000000000000 0.153000000000000 1.55000000000000;470 0.196000000000000 0.185000000000000 1.32000000000000;475 0.132000000000000 0.220000000000000 1.03000000000000;480 0.0805000000000000 0.254000000000000 0.772000000000000;485 0.0411000000000000 0.298000000000000 0.570000000000000;490 0.0162000000000000 0.339000000000000 0.415000000000000;495 0.00510000000000000 0.395000000000000 0.302000000000000;500 0.00380000000000000 0.461000000000000 0.219000000000000;505 0.0154000000000000 0.531000000000000 0.159000000000000;510 0.0375000000000000 0.607000000000000 0.112000000000000;515 0.0714000000000000 0.686000000000000 0.0822000000000000;520 0.118000000000000 0.762000000000000 0.0607000000000000;525 0.173000000000000 0.823000000000000 0.0431000000000000;530 0.237000000000000 0.875000000000000 0.0305000000000000;535 0.304000000000000 0.924000000000000 0.0206000000000000;540 0.377000000000000 0.962000000000000 0.0137000000000000;545 0.452000000000000 0.982000000000000 0.00790000000000000;550 0.530000000000000 0.992000000000000 0.00400000000000000;555 0.616000000000000 0.999000000000000 0.00110000000000000;560 0.705000000000000 0.997000000000000 0;565 0.794000000000000 0.982000000000000 0;570 0.879000000000000 0.956000000000000 0;575 0.951000000000000 0.915000000000000 0;580 1.01000000000000 0.869000000000000 0;585 1.07000000000000 0.826000000000000 0;590 1.12000000000000 0.777000000000000 0;595 1.13000000000000 0.720000000000000 0;600 1.12000000000000 0.658000000000000 0;605 1.09000000000000 0.594000000000000 0;610 1.03000000000000 0.528000000000000 0;615 0.951000000000000 0.462000000000000 0;620 0.856000000000000 0.398000000000000 0;625 0.755000000000000 0.340000000000000 0;630 0.648000000000000 0.284000000000000 0;635 0.535000000000000 0.228000000000000 0;640 0.432000000000000 0.180000000000000 0;645 0.344000000000000 0.140000000000000 0;650 0.268000000000000 0.108000000000000 0;655 0.204000000000000 0.0812000000000000 0;660 0.153000000000000 0.0603000000000000 0;665 0.112000000000000 0.0441000000000000 0;670 0.0813000000000000 0.0318000000000000 0;675 0.0579000000000000 0.0226000000000000 0;680 0.0409000000000000 0.0159000000000000 0;685 0.0286000000000000 0.0111000000000000 0;690 0.0199000000000000 0.00770000000000000 0;695 0.0136000000000000 0.00540000000000000 0;700 0.00960000000000000 0.00370000000000000 0;400 0.0191000000000000 0.00200000000000000 0.0860000000000000];
cie_1964_relative = [400 0.178338000000000 0.0186740000000000 0.802988000000000;405 0.177215000000000 0.0183750000000000 0.804410000000000;410 0.175544000000000 0.0182380000000000 0.806218000000000;415 0.173538000000000 0.0178460000000000 0.808615000000000;420 0.170919000000000 0.0178420000000000 0.811239000000000;425 0.168307000000000 0.0187360000000000 0.812956000000000;430 0.165467000000000 0.0203290000000000 0.814204000000000;435 0.162167000000000 0.0224680000000000 0.815365000000000;440 0.158934000000000 0.0257030000000000 0.815364000000000;445 0.155316000000000 0.0299800000000000 0.814705000000000;450 0.151398000000000 0.0365230000000000 0.812079000000000;455 0.146020000000000 0.0451260000000000 0.808855000000000;460 0.138532000000000 0.0587160000000000 0.802752000000000;465 0.129790000000000 0.0781810000000000 0.792029000000000;470 0.115226000000000 0.108760000000000 0.776014000000000;475 0.0955140000000000 0.159190000000000 0.745297000000000;480 0.0727520000000000 0.229553000000000 0.697695000000000;485 0.0452100000000000 0.327797000000000 0.626994000000000;490 0.0210330000000000 0.440145000000000 0.538821000000000;495 0.00726400000000000 0.562598000000000 0.430138000000000;500 0.00555700000000000 0.674174000000000 0.320269000000000;505 0.0218320000000000 0.752764000000000 0.225404000000000;510 0.0495700000000000 0.802379000000000 0.148050000000000;515 0.0850400000000000 0.817056000000000 0.0979040000000000;520 0.125439000000000 0.810035000000000 0.0645260000000000;525 0.166490000000000 0.792032000000000 0.0414780000000000;530 0.207440000000000 0.765864000000000 0.0266960000000000;535 0.243473000000000 0.740029000000000 0.0164980000000000;540 0.278702000000000 0.711170000000000 0.0101280000000000;545 0.313475000000000 0.681046000000000 0.00547900000000000;550 0.347313000000000 0.650066000000000 0.00262100000000000;555 0.381165000000000 0.618155000000000 0.000681000000000000;560 0.414219000000000 0.585781000000000 0;565 0.447072000000000 0.552928000000000 0;570 0.479019000000000 0.520981000000000 0;575 0.509646000000000 0.490354000000000 0;580 0.537520000000000 0.462480000000000 0;585 0.564346000000000 0.435654000000000 0;590 0.590406000000000 0.409594000000000 0;595 0.610811000000000 0.389189000000000 0;600 0.629921000000000 0.370079000000000 0;605 0.647268000000000 0.352732000000000 0;610 0.661104000000000 0.338896000000000 0;615 0.673036000000000 0.326964000000000 0;620 0.682616000000000 0.317384000000000 0;625 0.689498000000000 0.310502000000000 0;630 0.695279000000000 0.304721000000000 0;635 0.701180000000000 0.298820000000000 0;640 0.705882000000000 0.294118000000000 0;645 0.710744000000000 0.289256000000000 0;650 0.712766000000000 0.287234000000000 0;655 0.715288000000000 0.284712000000000 0;660 0.717300000000000 0.282700000000000 0;665 0.717489000000000 0.282511000000000 0;670 0.718833000000000 0.281167000000000 0;675 0.719255000000000 0.280745000000000 0;680 0.720070000000000 0.279930000000000 0;685 0.720403000000000 0.279597000000000 0;690 0.721014000000000 0.278986000000000 0;695 0.715789000000000 0.284211000000000 0;700 0.721805000000000 0.278195000000000 0;400 0.178338000000000 0.0186740000000000 0.802988000000000];
   
% outer curved boundary of the color space perceivable by human eye after CIE Lab 1931 and CIE Lab 1964
p(1) = plot(cie_1931_relative(1:end,2), cie_1931_relative(1:end,3), '-o','LineWidth', 1.0, 'Color', [0.7 0.7 0.7],'MarkerEdgeColor','k','MarkerFaceColor',[0.8, 0.8, 0.8],'MarkerSize',3.0);
p(2) = plot(cie_1964_relative(1:end,2), cie_1964_relative(1:end,3), '-o','LineWidth', 1.0, 'Color', [0.7 0.7 0.7],'MarkerEdgeColor','k','MarkerFaceColor',[0.95, 0.95, 0.95],'MarkerSize',3.0);

p(3) = plot([sf_xcoordinate;sf_xcoordinate(1,1)], [sf_ycoordinate;sf_ycoordinate(1,1)], '-k', 'LineWidth', 1.0); % color gamut border lines

for aa = 1:1:size(sf_xcoordinate,1)
    p(4) = plot(sf_xcoordinate(aa,1), sf_ycoordinate(aa,1), 'o','MarkerEdgeColor','k','MarkerFaceColor',[sf_RGB(aa,1),sf_RGB(aa,2),sf_RGB(aa,3)],'MarkerSize',7.0); % SF spots
end

p(5) = plot(sf_xcoordinate_mix, sf_ycoordinate_mix, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', [mix_RGB(1,1),mix_RGB(1,2),mix_RGB(1,3)],'MarkerSize',7.0); % mixed color resulting from blending single SF, chromaticity coordinates calculated in DATA ENTRY MASK
p(6) = plot(0.3067,0.3180, 'x','LineWidth', 2.5,'MarkerEdgeColor','k','MarkerFaceColor', 'k','MarkerSize',6.0); % modified D65 (see confluence)

p(7) = plot(xyL_table(1:end,end-2), xyL_table(1:end,end-1), '--', 'Color', [0.85 0.3 0.3], 'LineWidth', 2.0); % plotting horizontal CBU color profile

for aa = 1:1:size(sf_xcoordinate,1)
    p(8) = plot([sf_xcoordinate_mix,sf_xcoordinate(aa,1)],[sf_ycoordinate_mix,sf_ycoordinate(aa,1)], '--k', 'LineWidth', 0.5); % connecting lines between SF spots and mixed color spots
end

%{ 
p(9) = plot(0.4476,0.4074, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % CIE-A >> % following some further standard lighting values are attached
p(10) = plot(0.3484,0.3516, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % CIE-B
p(11) = plot(0.3101,0.3162, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % CIE-C
p(12) = plot(0.3333,0.3333, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % CIE-E
p(13) = plot(0.3457,0.3585, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % D50
p(14) = plot(0.3324,0.3474, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % D55
p(15) = plot(0.3127,0.3290, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % D65
p(16) = plot(0.2990,0.3149, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % D75
p(17) = plot(0.2848,0.2932, 'o','LineWidth', 1.5,'MarkerEdgeColor','k','MarkerFaceColor', 'y','MarkerSize',4.0); % D93
%}

axis equal;

xlim([0 0.8]);
set(gca,'xTick',0:0.1:0.8);
xlabel('X COORDINATE','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));
set(gca,'TickDir','out');
xtickangle(-90);

ylim([0 0.9]);
set(gca,'yTick',0:0.1:0.9);
ylabel('Y COORDINATE','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));
ytickangle(0);

grid on;
ax = gca;
ax.GridColor = [0.7, 0.7, 0.7];

title ('COLOR PROFILE', 'Fontsize', 16);

if subframe_number == 2
    legend ({'spectral locus (cie lab 1931, 2.0deg)','spectral locus (cie lab 1964, 10.0deg)','color gamut','first primary color (SF1)','second primary color (SF2)','mixed color', 'modified D65', 'horizontal color profile'},'Fontsize', 7.5,'Location','northeastoutside');
elseif subframe_number == 3
    legend ({'spectral locus (cie lab 1931, 2.0deg)','spectral locus (cie lab 1964, 10.0deg)','color gamut','first primary color (SF1)','second primary color (SF2)','third pimary color (SF3)','mixed color', 'modified D65', 'horizontal color profile'},'Fontsize', 7.5,'Location','northeastoutside');
elseif subframe_number == 4
    legend ({'spectral locus (cie lab 1931, 2.0deg)','spectral locus (cie lab 1964, 10.0deg)','color gamut','first primary color (SF1)','second primary color (SF2)','third pimary color (SF3)','fourth primary color (SF4)','mixed color', 'modified D65', 'horizontal color profile'},'Fontsize', 7.5,'Location','northeastoutside');
elseif subframe_number == 5
    legend ({'spectral locus (cie lab 1931, 2.0deg)','spectral locus (cie lab 1964, 10.0deg)','color gamut','first primary color (SF1)','second primary color (SF2)','third pimary color (SF3)','fourth primary color (SF4)','fifth primary color (SF5)','mixed color', 'modified D65', 'horizontal color profile'},'Fontsize', 7.5,'Location','northeastoutside');
elseif subframe_number == 6
    legend ({'spectral locus (cie lab 1931, 2.0deg)','spectral locus (cie lab 1964, 10.0deg)','color gamut','first primary color (SF1)','second primary color (SF2)','third pimary color (SF3)','fourth primary color (SF4)','fifth primary color (SF5)','sixth primary color (SF6)','mixed color', 'modified D65', 'horizontal color profile'},'Fontsize', 7.5,'Location','northeastoutside');
end

set(gcf,'WindowStyle','docked');

hold off;


f6 = figure (6); hold on; set(gca,'FontSize',11); % COLOR MAP (OVERVIEW) - plotting colored illustration of perceived CBU profile on retina for the single frame cycle with largest cbu effect (global map)

RGB_table = xyz2rgb(XYZ_table(1:end,end-2:end),'OutputType','uint8'); % transfering XYZ tristimulus values to RGB for colored illustration of CBU profile

for aa = 2:1:size(RGB_table,1) % loop for plotting colored scanning lines for every retinal spot (X_VALUES) representing the perceived color perception at the current retinal position
    
p(aa) = line([xyL_table(aa,1) xyL_table(aa,1)],[(subframe_location_table_ver_x1(2,1)) (subframe_location_table_ver_x1(1,1))] ,'Color',[RGB_table(aa,1),RGB_table(aa,2),RGB_table(aa,3)],'LineStyle','-','LineWidth', 0.5);

end

location_table_short_sorted = location_table_short(:);

for xx = 1:1:integrated_subframes.*2 % loop for plotting of color sector threshold lines
    
    p(aa+xx) = plot([location_table_short_sorted((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+xx,1) location_table_short_sorted((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+xx,1)],[(subframe_location_table_ver_x1(2,1))-1 (subframe_location_table_ver_x1(1,1))+1],'Color','k','LineStyle','--','LineWidth', 1.0);
    
end

for r1 = 2.5:2.5:90 % loop for plotting circle grid (2.5DEG pale grid)
    
n = 1000;
t = linspace(0,2*pi,n);
x = r1*sin(t);
y = r1*cos(t);
p(aa+xx+r1/2.5) = line(x,y,'LineStyle','--','Color', [0.85, 0.85, 0.85],'LineWidth', 0.5);

end

for r2 = 5:5:90 % loop for plotting circle grid (5.0DEG grid)
    
n = 1000;
t = linspace(0,2*pi,n);
x = r2*sin(t);
y = r2*cos(t);
p(aa+xx+r1/2.5+r2/5) = line(x,y, 'LineStyle','--', 'Color', [0.75, 0.75, 0.75],'LineWidth', 0.5);

end

p(aa+xx+r1/2.5+r2/5+1) = line([0 0],[(-viewing_angle_ver_deg) (viewing_angle_ver_deg)],'LineStyle','-.', 'Color', [0.65, 0.65, 0.65],'LineWidth', 1.0); % zero lines of retinal positioning for horizontal and vertical direction (intersection point represents fovea centralis)
p(aa+xx+r1/2.5+r2/5+2) = line([(-viewing_angle_hor_deg) (viewing_angle_hor_deg)],[0 0],'LineStyle','-.', 'Color', [0.65, 0.65, 0.65],'LineWidth', 1.0);

p(aa+xx+r1/2.5+r2/5+3) = plot(min(subframe_location_table_hor_x1(:)),subframe_location_table_ver_x1(1,1)-content_height_deg./2,'x','LineWidth', 2.0,'MarkerEdgeColor','r','MarkerSize',8.0); % plots the spatial limit of retinal stimulation in horizontal direction during the whole sequence (mathematically smallest position value)
p(aa+xx+r1/2.5+r2/5+4) = plot(max(subframe_location_table_hor_x1(:)),subframe_location_table_ver_x1(1,1)-content_height_deg./2,'x','LineWidth', 2.0,'MarkerEdgeColor','r','MarkerSize',8.0); % plots the spatial limit of retinal stimulation in horizontal direction during the whole sequence (mathematically highest position value)

axis equal;

xlim([(-viewing_angle_hor_deg) (viewing_angle_hor_deg)]);
set(gca,'xTick',(-viewing_angle_hor_deg):(((viewing_angle_hor_deg)-(-viewing_angle_hor_deg))./10):(viewing_angle_hor_deg));
xlabel('HORZIONTAL POSITION [DEG]','Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));
set(gca,'TickDir','out');
xtickangle(-90);

ylim([(-viewing_angle_ver_deg) (viewing_angle_ver_deg)]);
set(gca,'yTick',(-viewing_angle_ver_deg):(((viewing_angle_ver_deg)-(-viewing_angle_ver_deg))./10):(viewing_angle_ver_deg));
ylabel('VERTICAL POSITION [DEG]','Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));
ytickangle(0);

legend (p([(aa+1) (aa+xx+1) (aa+xx+r1/2.5) (aa+xx+r1/2.5+r2/5+1) (aa+xx+r1/2.5+r2/5+3)]),['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['retinal 2.5deg circle grid'],['retinal 5.0deg circle grid'],['retinal zero lines'],['horizontal spatial limits of retinal stimulation (all FC)'],'Fontsize', 7.5,'Location','northeastoutside');

title ('COLOR MAP (OVERVIEW)', 'Fontsize', 16);

set(gcf,'WindowStyle','docked');

hold off;


f7 = figure(7); hold on; set(gca,'FontSize',11); % COLOR MAP (DETAIL) - plotting colored illustration of perceived CBU profile on retina for the single frame cycle with largest cbu effect (detail map)

RGB_table = xyz2rgb(XYZ_table(1:end,end-2:end),'OutputType','uint8'); % transfering XYZ tristimulus values to RGB for colored illustration of CBU profile

for aa = 2:1:size(RGB_table,1) % loop for plotting colored scanning lines for every retinal spot (X_VALUES) representing the perceived color perception at the current retinal position
    
p(aa) = line([xyL_table(aa,1) xyL_table(aa,1)],[(subframe_location_table_ver_x1(2,1)) (subframe_location_table_ver_x1(1,1))] ,'Color',[RGB_table(aa,1),RGB_table(aa,2),RGB_table(aa,3)],'LineStyle','-','LineWidth', 0.5);

end

location_table_short_sorted = location_table_short(:); 

for xx = 1:1:integrated_subframes.*2 % loop for plotting of color sector threshold lines
   
    p(aa+xx) = plot([location_table_short_sorted((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+xx,1) location_table_short_sorted((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+xx,1)],[(subframe_location_table_ver_x1(2,1))-1 (subframe_location_table_ver_x1(1,1))+1],'Color','k','LineStyle','--','LineWidth', 1.0);

end

for r1 = 0.25:0.25:2 % loop for plotting circle grid (0.25DEG grid)
    
n = 1000;
t = linspace(0,2*pi,n);
x = r1*sin(t);
y = r1*cos(t);
p(aa+xx+r1/0.25) = line(x,y,'LineStyle','--','Color', [0.65, 0.65, 0.65],'LineWidth', 0.5);

end

for r2 = 0.5:0.5:2 % loop for plotting circle grid (0.5DEG grid)
    
n = 1000;
t = linspace(0,2*pi,n);
x = r2*sin(t);
y = r2*cos(t);
p(aa+xx+r1/0.25+r2/0.5) = line(x,y, 'LineStyle','--', 'Color', [0.75, 0.75, 0.75],'LineWidth', 0.5);

end

p(aa+xx+r1/0.25+r2/0.5+1) = line([0 0],[min(subframe_location_table_ver_x1)-(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10 max(subframe_location_table_ver_x1)+(max(subframe_location_table_ver_x1)-min(subframe_location_table_ver_x1))/10],'LineStyle','-.', 'Color', [0.85, 0.85, 0.85],'LineWidth', 1.0); % zero lines of retinal positioning for horizontal and vertical direction (intersection point represents fovea centralis)
p(aa+xx+r1/0.25+r2/0.5+2) = line([min(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))-(max(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))-min(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))/10 max(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))+(max(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))-min(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))/10], [0 0], 'LineStyle','-.', 'Color', [0.85, 0.85, 0.85],'LineWidth', 1.0);

axis equal;

xlim([(min(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))) (max(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))]);
set(gca,'xTick',(min(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))):(((max(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1)))-(min(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))))./10):(max(location_table_short_sorted(((2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+1:(2.*integrated_subframes.*(eye_velo_hor_variable_degs_table_index-1))+2.*integrated_subframes),1))));
xlabel('HORIZONTAL POSITION [DEG]', 'Fontsize',14);
set(gca, 'XTicklabel', num2str(get(gca, 'XTick')', '%.2f'));
set(gca,'TickDir','out');
xtickangle(-90);

ylim([(subframe_location_table_ver_x1(2,1)) (subframe_location_table_ver_x1(1,1))]);
set(gca,'yTick',(subframe_location_table_ver_x1(2,1)):(((subframe_location_table_ver_x1(1,1))-(subframe_location_table_ver_x1(2,1)))./10):(subframe_location_table_ver_x1(1,1)));
ylabel('VERTICAL POSITION [DEG]', 'Fontsize',14);
set(gca, 'YTicklabel', num2str(get(gca, 'YTick')', '%.2f'));
ytickangle(0);

legend (p([(aa+1) (aa+xx+1) (aa+xx+r1/0.25) (aa+xx+r1/0.25+r2/0.5+1)]),['color sector threshold lines (FC' num2str(eye_velo_hor_variable_degs_table_index) ')'],['retinal 0.25deg circle grid'],['retinal 0.5deg circle grid'],['retinal zero lines'],'Fontsize', 7.5,'Location','northeastoutside');

title ('COLOR MAP (DETAIL)', 'Fontsize', 16);

set(gcf,'WindowStyle','docked');

hold off;

end