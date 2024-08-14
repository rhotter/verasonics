% Notice:
% This file is provided by Verasonics to end users as a programming example for the Verasonics Vantage Research
% Ultrasound System. Verasonics makes no claims as to the functionality or intended application of this program and
% the user assumes all responsibility for its use.
%
% File name: SetUpL11_5vWideBeam_Annotation.m - for L11_5v Linear array wide beam transmit. This file is provided to
% with annotation regarding the reference in the programming manual, tutorial, or training video.
%
% Description:
% The transmit aperture size can be set by the variable - numTx. The transmit aperture is translated across the array
% with the number of rays set by the variable - numRays. The receive aperture always covers the full 128 element
% aperture.
%
% Before going through this annotated script, it is highly recommended to go through our tutorial video on our web
% portal.
%
% This script generates all the sequence objects needed to acquire and process ultrasound images produced with an
% unfocused transmit pulse from a typical linear array transducer, in this case, a Verasonics L11-5v transducer. Each
% image consists of overlapping transmit beams, Wide Beam Scan, with a focus below bottom of field. Wide Beams and
% Regions overlap as they are scanned across field. On receive, all 128 channels are used to acquire the echo data from
% the full field of view. Intensities of reconstructed pixels are normalized by transmit beam intensity.
%
% This sequence illustrates the use of asynchronous acquisition and processing, where the frame rate of acquisition
% doesn't necessarily match the frame rate of processing. The acquisition of frames into memory can be at a fairly high
% rate, while the rate of processing can be at something reasonable for real-time viewing. A cineloop for both acquired
% RF data and processed image data will be created, either of which can be reviewed in a playback mode.
%
% update history:
% 10-01-2018 First Annotated script

% Clear all is required to remove items from workspace, freeing up system memory.
% More useful information can be found: https://www.mathworks.com/help/matlab/ref/clear.html
clear all

%% 1. Define system parameters (Section 2.1 in Programming Manual)
% P structure is the abbreviation of "PreSet". Parameters that are prefaced with a P are saved along with other values
% when a Preset is saved. If the user-defined variables will be modified by the callbacks of any UI controls, these
% variables should be placed in the P structure, that will be saved in the preSet file, to make preSetTool function
% correctly. It also helps user to manage the predefined variables.
P.numTx = 32;   % no. of elements in TX aperture.
P.numRays = 48; % no. of rays in frame, it's also no. regions of a image
P.startDepth = 5; % wavelength, used in Receive
P.endDepth = 192; % wavelength, used in Receive

% Define the number of transmitters available in the system. For a 128 channel Vantage system, or a 256 channel system
% using a single connector, this number is 128.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.

% Define the number of receive channels available in the system. For a 128 Vantage system, or a 256 channel Vantage
% system using one connector, this number is 128.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.

% In the example text above, we are setting the parameters for a 128 channel Vantage Unit, which has 128 transmitters
% and 128 receive channels. A 64 channel Vantage Unit would specify Resource.Parameters.numTransmit = 64 and
% Resource.Parameters.numRcvChannels = 64, instead of 128. For a Vantage 64LE, only Resource.Parameters.numRcvChannels
% would be set to 64.

% Define the speed of sound in the media to be imaged in meters/sec. For medical imaging, the speed of sound is an
% average of the various tissues, typically 1540 m/s. This parameter, if defined differently from the default of 1540
% m/s, should be defined before computing the transducer characteristics in Trans.
Resource.Parameters.speedOfSound = 1540; % speed of sound in m/sec

% The verbose attribute sets the level at which VSX and other functions report errors, warnings, and status messages on
% the Matlab command line. The supported values are:
% 0: Error messages only are displayed.
% 1: Error and warning messages are displayed.
% 2: Error, warning, and status messages are displayed. (default)
% 3. Error, warning, status, and debug status messages are displayed.
Resource.Parameters.verbose = 2;

% This attribute, when set to 1, causes the software to exit after initialization. Initialization adds all the missing
% default attributes and validates structures. Additional hidden structures such as DMAControl for specifying transfers
% of data to host memory are also created. This is useful for debugging a script that crashes Matlab when run. Once
% initialized, one can examine all the structures for added attributes and correct sequencing. It has been incorporated
% wit the EventAnalysisTool. If the EventAnalysisTool is used right after running the SetUp script (NOT VSX),
% EventAnalysisTool will use initialized only feature to check possible errors in the script.
Resource.Parameters.initializeOnly = 0;

% The simulateMode attribute specifies whether you want to run your script in simulate mode, or with acquired data from
% the Vantage hardware. A value of 0 means that you would like to use the Vantage hardware for acquisition. If the
% hardware is not detected, VSX will automatically switch to simulate mode, which is set with a value of 1. In simulate
% mode, you can choose whether to simulate the RF data from a Media model which you can define (simulateMode = 1), or
% whether to simply process whatever data is already in the RcvBuffer memory (simulateMode = 2). VSX creates a GUI
% toggle control to switch to simulateMode = 2 (named RcvDataLoop) from either simulateMode = 0 or 1.
Resource.Parameters.simulateMode = 0;
%  Resource.Parameters.simulateMode = 1 forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 stops sequence and processes RcvData continuously.

%% 2. Define the Transducer characteristics. (Section 2.2 in Programming Manual)
Trans.name = 'L11-5v';
Trans.units = 'mm';
Trans = computeTrans(Trans);

% Please first go through the Training Video 3_Defining the Transducer Parameters
%
% The Trans structure array defines the characteristics of the transducer attached to the system or used for simulation.
% If the transducer is one that is known by the software, which is the case with the L11-5v, it's attributes can be
% defined using the computeTrans.m function in the Utilities folder, and the full Trans structure will be created in the
% Matlab workspace, which can then be viewed with the Matlab workspace variable viewer or by typing "Trans" at the
% Matlab command prompt.
%
% For better understanding, some of the added attributes are briefly described below. Please refer to Section 2.2,
% Transducer Object, regarding detailed description, like transducer geometry, element position definition, and UTA
% type.
%
% Trans.frequency is the center frequency of the transducer, Fc. It specifies the conversion of wavelengths to time and
% or distance (given the speed of sound). Since time and distance attributes are specified in wavelengths, it is
% important to provide an appropriate center frequency. One can specify a different center frequency from the default
% frequency if desired, but this attribute should be set BEFORE calling computeTrans, since calculations of wavelengths
% will depend on this value.
%
% The elementWidth, ElementPos, and lensCorrection attributes are in units specified by Trans.units.
%
% Trans.maxHighVoltage sets the maximum high voltage that the system will allow to be applied to the transducer (+/-
% maxHighVoltage). Trans.maxHighVoltage can be set by the user before or after calling computeTrans, since computeTrans
% will not overwrite this attribute. CAUTION: The maximum hardware high voltage limit is capable of damaging some
% transducers and/or producing harmful acoustic power levels in biological tissue. To prevent damage, do not modify the
% high voltage limits in computeTrans. In scanning, use high voltage levels only as high as needed for a specific
% application.
%
% Trans.type determines the transducer geometry. A transducer type of 0 specifies a linear array. Other types supported
% are 1 (curved linear), and 2 (2D array). A linear array has elements whose elements can be specified with a single x
% coordinate value. Please refer to Section 2.2, Transducer Object regarding the transducer geometry and element
% position.
%
% Trans.connType identifies the connector type required by the user script. it is an integer value used as an index to
% identify one of the connector types supported by the Universal Transducer Adapter system.
%
% The Trans.spacing is sometimes known as pitch. The element spacing is provided by computeTrans.m in both mm and
% wavelengths.
%
% The position of all transducer elements in Trans.ElementPos. This array has a row for each element with four values:
% x, y, z, and alpha. The angle alpha is in radians, and specifies the angle of the normal of the element with the z
% axis. This parameter is typically only used for curved linear arrays. The x coordinate value is the only parameter we
% need to set for a linear array, and it is set to the position of the center of the element. The z axis passes through
% the center of the transducer array, between the central elements.
%
% It should be noted that for image reconstruction, the position of transducer elements is solely determined by the
% Trans.ElementPos array. Other attributes, such as Trans.spacing(Mm) and radius(Mm) can be defined for convenience in
% calculations, but are not used by the image reconstruction software.
%
% Trans.ElementSens defines the element sensitivity curve, sometimes known as a directivity pattern is computed using
% the above formula. This curve represents the fall off in intensity of an element's response to an echo coming in at an
% angle theta from the normal of the element (and also the fall off in intensity of the transmit intensity with angle).
% In this case, we are using a formula that calculates the curve for a narrow bar type element. If we knew the exact
% curve, as could be measured with a hydrophone in a water tank, we could enter these values directly.
%
% Trans.impedance is used to compute power limits. Verasonics provides an option for impedance measurement. Please refer
% to the Application Note in the Tools/ImpedanceMeasure folder for more details.
%
% Trans.Connector is an dependent attribute generated by VSX. It defines the composite "Element to Channel" mapping. The
% value of the first entry defines the channel number connected to the first element.

%% 3. Define the Pixel Data region for reconstruction. (Section 5.3 in Programming Manual)
% Please first go through the Training Video 4_Defining the Imaging Region
%
% The PData object (short for "pixel data") describes the pixel region(s) to be processed by the image reconstruction
% software. If no image reconstruction is required, this object can be omitted from a setup script. When image
% reconstruction is performed, the PData object specifies one or more pixel regions in a larger image space, which is
% defined relative to the transducer coordinate system. Image reconstruction is only performed at the pixel or voxel
% locations included in the regions specified. The resolution is determined by multiple factors, such as the transducer
% frequency, the size of the transducer aperture, and the length of the transmit burst. Typically,to sample the best
% axial and lateral resolution of the majority of ultrasound transducers with our pixel array, it is usually sufficient
% to use a pixel spacing in the lateral direction of one wavelength (which is typically close to the element spacing),
% and in the axial direction of half a wavelength. Note: The pixel density for the PData pixel array is not the same as
% the pixel density used for displaying the ultrasound image. Most displays have a very high resolution of pixels and
% thus require scaling up the pixel density for display.
%
% The rectangular coordinate PData specification for linear array is shown in Fig. 2.4.1 in the programming manual and
% The structure definition is as follows:

% The PData.PDelta specifies the spacing between pixels (in wavelengths) in the X, Y and Z dimensions. Set the
% PDelta(pdeltaX,0,pdeltaZ) attributes to a wavelength increment per pixel that is small enough to adequately sample the
% resolution of the image. We use separate pdeltaX and pdeltaZ values for the L11-5v linear array, since the depth
% resolution is typically much better than the lateral resolution. Note that the smaller the values, the longer the
% image or signal reconstruction will take for a given region.
PData(1).PDelta = [Trans.spacing/2, 0, 0.5];

% The PData.Size attribute is defined in rows, columns, and sections (for 3D volumes). NOT in X, Y and Z dimensions. In
% polar coordinates, the rows correspond to the R dimension, while the columns correspond to the Theta dimension. Please
% see programming manual, section 2.4 regarding the polar coordinate. Only rectangular coordinate will be described
% here.
PData(1).Size(1,1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % range, origin and pdelta set PData(1).Size.
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(1,3) = 1;      % single image page

% The PData.Origin attribute defines the location of the PData region in the coordinate system of the transducer. The
% PData.Origin is location of the upper left corner of the region. In this case we set the x dimension of the Origin to
% the center of the first transducer element. The z dimension is set to the start of the scan, as given by P.startDepth.
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.

% While it is possible to specify a PData array that selects a specific local region of the field of view for
% processing, it is preferable to use the PData.Region definition for this purpose. A number of geometric shapes can be
% defined, using a "Shape" structure. Please refer to section 2.4.1 PData.Region Objects in the programming manual for
% supported regions. The illustration of a Rectangle Shape Region within the PData area is shown in Fig. 12.3 in
% tutorial.

% Define P.numRays rectangular regions centered on TX beam origins.
PData(1).Region = repmat(struct('Shape',struct('Name','Rectangle',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',Trans.spacing*(P.numTx-12),...
                                               'height',P.endDepth-P.startDepth)),1,P.numRays);

% Multiple regions are required for one image, so the PData.Region(n).Shape.Position(1) is different for each region.
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end

% If the Region structure contains a supported Shape attribute, the utility function "computeRegions.m" can be used to
% compute the numPixels and PixelsLA attributes for the Region. Note: If no PData.Region is specified, the
% computeRegions utility function will create a default Region, which is the same size as the PData array.
% After running the SetUp script, please use showTXPD to evalute the region definition by clicking the "showRegion".
PData(1).Region = computeRegions(PData(1));

%% 4. Define a media model to use when in simulate mode. (Section 2.5 in Programming Manual)
% To run in simulate mode, we need to define a Media model to generate ultrasound echoes. The Media model is simply a
% set of point targets that have a defined spatial location and reflectivity. The spatial location is defined in a three
% dimensional coordinate system tied to the transducer - in the case of the linear array, the transducer elements are
% located along the x axis, centered around x=0, and the scan depth dimension is along the z axis. The y axis represents
% the elevational dimension, and is not used for a 2D linear scan.

% In this script, we call a separate Matlab script to define the media model, named 'pt1'. This script defines the Model
% name, the location of the point targets along with their reflectivity, and the number of point targets.
pt1;

% The Media.attenuation defines the the attenuation of the media in dB/cm/ MHz. A typical value is -0.5 dB/cm/MHz. Note
% that the attenuation value is used by both the simulation and the image reconstruction software to compensate for the
% attenuation loss on transmit only. The attenuation loss on receive is compensated by the Time Gain Control curve
% (TGC). With Doppler script, sometimes, the noise signal is higher than expected. Lower the Media.attenuation will be
% one possible solution.
Media.attenuation = -0.5;

% The Media.function attribute is used for dynamic media models, where the targets change their positions or
% characteristics between successive acquisitions. When present, and the Receive.callMediaFunc attribute is true (=1),
% the Matlab function (m file) specified by the string is called prior to simulation of the received data. For instance,
% 'movePoints' is used to move the all the points in the x direction after each frame. When to call this function is
% specified in the Receive structures to be defined later.
Media.function = 'movePoints';

%% 5. Define resources used by the sequence. (Section 2.6 in Programming Manual)
% Before describing the Resource attributes more completely, a few words about the data flow in the Verasonics system
% are in order. The Verasonics system can best be understood as implementing data transformations from one memory buffer
% type to another. The order of data flow with the various processes is as follows:
%
% Acquisition - RF data acquired in Vantage hardware modules' local channel memories are transferred to a RcvBuffer in
% the computer host's memory. The Vantage system provides large amounts of local memory (typically 48 MBytes or more for
% each channel)
%
% Pixel reconstruction - Processing takes input RF data from a RcvBuffer and output goes to InterBuffer (for I,Q pixels)
% and/or to ImageBuffer (for intensity pixels).
%
% Processing - From one InterBuffer or ImageBuffer to another (or perhaps the same buffer).
%
% Display - From one or more ImageBuffers to the DisplayWindow buffer. The DisplayWindow typically has more pixels than
% specified in the corresponding ImageBuffer, since the pixels often need to be interpolated up to a higher density to
% provide a useful image size on a high resolution display.
%
% The various Resource buffer definitions create actual Matlab arrays which hold the respective data. The corresponding
% Matlab arrays created have the following names:
%
% RcvBuffer(:) - Matlab cell array RcvData{n}(i,j,k), where n is the buffer number, i is the row index, j is the column
% index, and k is the frame no.
%
% InterBuffer(:) - Matlab cell array IQData{n}(i,j,k,l,m), where n is the buffer number, i is the row index, j is the
% column index, k is the section number, l is the page number, and m is the frame number.
%
% ImageBuffer(:) - Matlab cell arrays ImgData{n}(i,j,k,m) and ImgDataP{n}(i,j,k,m), where n is the buffer number, i is
% the row index, j is the column index, k is the section number, and m is the frame number.

% Set the datatype of the buffer to 16 bit signed integers. This is currently the only choice.
Resource.RcvBuffer(1).datatype = 'int16';

% Set the rowsPerFrame attribute to something larger than the sum of all the channel acquisitions in the frame, since
% the acquisitions will be stacked down a column. In this case, there are P.numRays acquisitions per frame. Considering
% that a maximum depth per acquisition would typically be less than 512 wavelengths (1024 round trip), or 4096 samples
% per acquisition (at 4 samples per wavelength), we can set a maximum rowsPerFrame of P.numRays*4096 samples.
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*4096;

% Set the colsPerFrame attribute to the number of channels available, in this case, 128.
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;

% Set the numFrames attribute to the number of frames to be captured in the RcvBuffer. The total number of frames should
% be greater than the number of frames that can be acquired during the processing time for a frame. In this way, we
% won't overwrite the data in a frame that is being processed. For example, if acquisition is set to acquire 100 frames
% per second, and we are able to process 50 frames per second, acquisition will acquire two frames while we are
% processing a frame. In this case, we should have at least two frames in our RcvBuffer. The organization of the buffer
% is shown in Fig. 2.6.1.1 in the programming manual. A RcvBuffer can contain one frame or an EVEN number of identical
% frames.
Resource.RcvBuffer(1).numFrames = 30;    % 30 frames stored in RcvBuffer.

% There is no need for an InterBuffer for storing the complex signal pixel reconstruction data if we are not using
% synthetic aperture reconstructions. We can reconstruct directly to intensity data in an ImageBuffer. The 10 frame
% ImageBuffer will receive the output of the image reconstruction, but won't have any of the additional processing
% needed for display, so please note that the ImgData is not the image shown on the displayWindow. Here, multiple
% acquisitions are required for one frame, so one inter buffer is necessary. If we don't define rows and cols, the
% InterBuffer and ImageBuffer will be sized the same as PData.Size. If the user would like to have access to the IQ data
% of each acquisition, or each frame, please refer to another example script, SetUpL11_5vFlashAngles_IQall.m, located in
% Specialty_Applications/ReconAll_IQData directory.
Resource.InterBuffer(1).numFrames = 1;   % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;

% The next Resource to define is a DisplayWindow for displaying the processed echo intensity image. This buffer will
% define the attributes of the Matlab display window used to show the image. The DisplayWindow typically has a higher
% pixel density than the ImageBuffer, as the ImageBuffer must only adequately sample the resolution of the ultrasound
% image, where the DisplayWindow must contain enough pixels to present a reasonable sized image on a high definition
% display. The pixels of the ImageBuffer will be interpolated up to the pixel density of the DisplayWindow. The
% DisplayWindow is created as a bit-mapped window with 8 bit pixels. A colormap is defined to set the mapping of echo
% intensity to gray scale values.

% Set the title that will appear in the frame of the display window.
Resource.DisplayWindow(1).Title = 'L11-5vWideBeam';

% Set the pixel density of the display window. This will determine the scale of the image data within the DisplayWindow.
% The larger the pdelta value, the smaller the image that will be rendered within the DisplayWindow.
Resource.DisplayWindow(1).pdelta = 0.35;

% Set the position of the Display window on the computer's display. To center the displayWindow vertically, we obtain
% the ScreenSize array, which provides the width and height of the computer display. The position parameters are
% [leftEdge, bottomEdge, width, height]. In this case, we set the size of the DisplayWindow (when converted to
% wavelengths) to match the size of the PData region (in wavelengths). We also set the DisplayWindow x and y reference
% point to match the PData.Origin x and z reference point. If we wanted a black border around our PData region on the
% display, we would set the DisplayWindow size somewhat bigger and the reference point slightly above and to the left of
% the PData.Origin.
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];

% The position relative to the transducer is given by the attribute, Resource.DisplayWindow.ReferencePt, whose values
% are in wavelengths, and specify the location of the upper left corner of the DisplayWindow in the transducer
% coordinate system. The Resource.DisplayWindow.mode attribute will be needed for when a 3D display routine is
% integrated. Please refer to section 2.6.4 in the programming manual and the example scripts in
% Specialty_Applications/Customized2DArray regarding the 3D display.
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane

% The additional display type was prompted by deficiencies in the Matlab display window for Matlab releases with 2015b
% or later, which include an update limit of 20 frames per second and a noticeable lag between the time of update and
% the appearance of the image on the display. The Verasonics viewer resolves these issues with a visible frame rate
% limit of 60 frames per second and a minimal delay in rendering the image to the display. The default display window is
% still the Matlab window, but the Verasonics display window can be specified by setting the Resource.DisplayWindow.Type
% attribute to 'Verasonics'.
Resource.DisplayWindow(1).Type = 'Verasonics';

% The numFrames attribute specifies the number of frames in the DisplayWindow cineloop buffer. If the user would like to
% have access to the display window, please refer to the nested function cineSave_Callback defined in the vsx_gui.m. It
% describe how the cineDisplay can be used for retrieve a particular frame from the cineloop buffer and how the frame
% can be saved with either Verasonics or Matlab viewer.
Resource.DisplayWindow(1).numFrames = 20;

% Set the type of units for the axes of the DisplayWindow. The default is 'wavelengths', but we can select 'mm' for
% millimeters.
Resource.DisplayWindow(1).AxesUnits = 'mm';

% Set the colormap for the DisplayWindow. The attribute here is a function that returns an array of 256x3 double values,
% with each row an R,G,B setting.
Resource.DisplayWindow(1).Colormap = gray(256);

% This completes the definition of the non-event structures of our sequence, that is, those that are not involved
% directly in specifying acquisition and processing actions. The remaining structures to be defined will be referenced
% directly or indirectly in our Event sequence by their index numbers. Obviously, one must already have a very good idea
% of what events will compose the sequence to define the referenced structures, so in some cases, the user may want to
% rough out a sequence of events before attempting to define the referenced items.

%% 6. Define the Transmit Waveform structure, TW. (Section 3.2.1 in Programming Manual)
% Please first go through the Training Video 5_Defining Transmit Parameters about the TW and TX definition
%
% We can now define the attributes of the objects that will be used in our sequence of events. The TW structure defines
% the characteristics of the transmit pulse. Each TW structure defines a single transmit waveform, so if we have several
% different waveforms, we have to create a separate TW structure for each one. However, for many scans, especially 2D
% echo scans such as defined by this script, it is sufficient to define a single waveform that is used for each transmit
% event.

% Set the type of transmit waveform definition to use. A simple transmit waveform that can be generated by the Vantage
% hardware can defined the type as 'parametric'. Other more complex transmit waveform types are 'envelope' and 'pulse
% code'. If one is writing a script for simulation only, additional types named 'function' and 'sampled' are available.
% Please check programming manual for more details. If the user is interested in our ArbWave capability, please refer to
% the ArbWaveToolbox in Tools folder and check our website: http://verasonics.com/arbitrary-waveform-generation-capability/

TW(1).type = 'parametric';

% Set the parameters for a parametric waveform. The values represent: A, the transmit frequency in MHz (which should
% have a period that can be implemented with an even number of 250MHz clock cycles); B, the duty cycle of a half cycle
% period, expressed as a fraction of 1.0 (used for transmit apodization); C, the number of half cycle periods in the
% transmit burst; D, the polarity of the first half cycle. In this instance, we want to set a frequency as close as
% possible to our transducer center frequency of 6.25 MHz, which happens to have a period of exactly 40 250MHz clock
% cycles. We set the half cycle on time to some fraction of the half cycle period to better approximate a sine wave, in
% this case 0.67, with that the third harmonic will be zero. For this transmit waveform, we want a short burst of 1
% cycle, so we set the number of half cycle periods to 2. Finally, we set the initial polarity of the first half cycle
% to 1, which is positive. In this case, we want all transmitters to use the same waveform, so we only have to specify a
% single set of parameters.
TW(1).Parameters = [Trans.frequency,0.67,1,1];

%% 7. Define the Transmit events structure, TX. (Section 3.2.2 in Programming Manual)
% Next we will define a transmit action using the TX structure. We need to define a TX structure for each unique
% transmit action in our sequence. Multiple transmits that have the same set of attributes can be specified with a
% single TX definition. The transmit waveform to use should be previously defined, and is referenced by the TX
% structure.
%
% The main attributes of the TX objects are shown below.

TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', 3*P.endDepth, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [], ...
                   'peakCutOff', 1.0, ...
                   'peakBLMax', 4.0), 1, P.numRays);

% TX.Origin, TX.Apod, TX.Delay and TX.TXPD need to be defined for each transmit beam with 'wavelength' unit
scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
end

% - Set event specific TX attributes.
for n = 1:P.numRays
    TX(n).Origin(1) = TxOrgX(n);
    % Compute transmit aperture apodization
    TX(n).Apod = +(((scaleToWvl*Trans.ElementPos(:,1))>(TxOrgX(n)-Trans.spacing*P.numTx/2))& ...
                 ((scaleToWvl*Trans.ElementPos(:,1))<(TxOrgX(n)+Trans.spacing*P.numTx/2)))';
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
end

% The TX.waveform is the index into an array of TW objects (see previous section) that define the transmit waveforms to
% be generated by each of the active transmitters.
%
% The TX.Origin Set the origin point for the transmit beam. In the case of a flat wavefront transmit, this is not
% particularly meaningful, and we just set the origin as the center of the aperture. The Origin is required for focused
% beams (TX.focus > 0), and virtual apex beams (TX.focus < 0) and specifies the point on the transducer element surface
% that the beam appears to originate from.
%
% The TX.Apod specifies the transmit apodization function to use for this transmit. For transmit events that use the
% Vantage hardware, the TX.Apod values can range from -1 to 1. If not all transmitters are needed for a transmit event,
% the TX.Apod array can be used to select a subset of the available transmitters, by specifying a zero value for the
% transmitters not used. The B parameter of TW.Parameters in the TW waveform structure reference by this TX will be
% modified by the TX.Apod values (the duty cycle fraction of the half cycle period will be multiplied by the TX.Apod
% value for each transmitter).
%
% TX.focus specifies the transmit focal point. The focus distance is the distance (in wavelengths) from the Origin point
% on the transducer to where the beam comes to a focus. A positive focus value defines the point of convergence for a
% multi-element transmit beam. A focus value of 0 specifies a flat wavefront across the transmit aperture, while a
% negative focus value specifies a diverging beam whose effective source is along the -z axis at the negative focal
% value (for no steering).
%
% TX.Steer specifies the beam steering to use for this transmit. The steering can be specified with two angles: the
% first entry specifies the angle of the beam projection into the x,z plane from the positive z axis (azimuth), and the
% second entry specifies the angle of the beam with respect to the x,z plane (elevation).
%
% If TX.Apod, TX.focus and TX.Steer are given (and also TX.Origin for TX.focus~=0), the method 'computeTXDelays' can be
% called to compute the individual transmit delays. This function also requires access to Trans, the global object
% specifying transducer geometry (defined earlier). The method computes TX.Delay values for all elements with a non-zero
% TX.Apod value. These delay times are computed in wavelengths of the Trans.frequency attribute, which translate to
% times based on the speed of sound defined in the Resource.Parameters.speedOfSound attribute.
%
% For some applications, it may be more convenient to define the transmit focus at a fixed location. In these cases, the
% TX.FocalPt attribute can be used instead of TX.Focus, TX.Steer, and TX.Origin when calling the computeTXDelays
% function. The TX.FocalPt attribute specifies a fixed coordinate, (x,y,z) where the transmit energy will be focused. If
% the TX.FocalPt attribute is missing or empty, the TX.focus attributes will be used, but if present, the TX.focus,
% TX.Steer and TX.Origin attributes will be ignored. If mm is desired, please use TX.FocalPtMm.
%
% The TXPD attribute of the TX object holds information on the simulated transmit field at the pixel locations defined
% in a PData array. A function is provided for computing the field parameters called computeTXPD. If added to a TX
% structure, a reconstruction that references the TX structure (through the ReconInfo.txnum attribute) will use the TXPD
% data in processing the reconstruction, and perform the follow:
%
% 1) The reconstruction output intensity will be normalized at each pixel location according to the peak transmit
% intensity at the pixel.
%
% 2) The time of travel from the origin of the transmit burst to the time of peak intensity at the pixel location will
% use the TXPD calculation instead of an algorithmic approximation.
%
% 3) Only pixels with a peak intensity greater than TX.peakCutOff (in transmit equivalents) and a burst length at the
% pixel location less than TX.peakBLMax will be reconstructed. These new TX attributes default to 1.25 and (the minimum
% transmit burst length + 0.25 (wavelengths)), respectively, if not provided.
%
% The TXPD attributes can be added to facilitate image reconstruction for certain types of scanning sequences,
% particularly those that involve multiple overlapping transmit beams. In these sequences, the pixels within the
% overlapping beams may be reconstructed multiple times, namely once for each transmit/receive event that insonifies the
% pixel.
%
% How to visulaize the transmit beam and inspect how the TX overlaps the PData? --> showTXPD
% Just type showTXPD in the command window

%% 8. Define the Time Gain Control waveform structure, TGC. (Section 3.3.3 in Programming Manual)
% We now can address our next desired action - "Receive echo signals from plane wave transmit and store in local
% memory." For this action we need to program the analog signal amplifiers on each receive channel and specify how the
% signals are digitized, filtered and stored. The TGC object defines the time-gain-compensation curve for the receive
% portion of the acquisition event.
%
% The TGC.CntrlPts array specifies the TGC level at 8 equally spaced points over the course of the Receive acquisition
% interval, starting at zero and ending at a depth set by TGC.rangeMax, in units of wavelengths of Trans.frequency.
% Therefore TGC.rangeMax is typically set equal to Receive.endDepth. The various control values will be able to be
% changed during run time by means of a set of TGC sliders on the user interface window, so the values given here are
% used to set the starting position of the sliders. The utility function, 'computeTGCWaveform(TGC)', computes the actual
% TGC curve which will be sent to the hardware. The TGC curve is re- computed whenever a control slider changes, and the
% new curve is sent to the hardware.

TGC.CntrlPts = [51,245,410,475,549,659,736,791];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% 9. Define the receive events structure, Receive. (Section 3.3 in Programming Manual)
% Please first go through the Training Video 6_Defining Receive Parameters Now we need to define a Receive operation
% that will specify other attributes of the input signal processing. The Receive structure defines the receiver
% characteristics of each acquisition event. Please refer to Fig. 3.3.1 in the programming manual for better
% understanding the role of each parameter in the receive signal path.
%
% There should be a separate Receive specification for each acquisition event in the sequence that goes to a unique
% location in the RcvBuffer. Since there are potentially a large number of Receive structures, it is best to predefine
% the structures with default attributes, then over-write the default attributes with the desired ones in a loop. The
% definition of the default set of attributes is as follows:

maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Resource.Parameters.numTransmit), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numRays*Resource.RcvBuffer(1).numFrames);

% The Receive.Apod attribute sets the individual receive channel gain. A value of 0 turns a receive channel off, writing
% all zeros to the channels memory. The value of one is set here for all channels, meaning we want full amplitude on
% channel outputs. The number of values in the Receive.Apod array is typically equal to the number of elements in the
% array for non-multiplexed transducers. The Receive.Apod array values are used to set the scaling multiplier after the
% Input Filter and the values can range from -4.0 to 4.0. Note that the length of Receive.Apod should be equal to
% Resource.Parameters.numTransmit and the length of non-zero values is dependent on the system. For instance, only 64
% channels are available for the 64LE, only 64 indices of the Receive.Apod array can be set to values other than zero,
% but the length of Receive.Apod is 128.
%
% The Receive.startDepth attribute determines the range at which A/D sampling begins, and is set to the same wavelength
% value set by P.startDepth.
%
% The Receive.endDepth attribute determines the range at which A/D sampling ends. For this attribute, we need to set a
% larger value than found in P.endDepth, since to reconstruct pixels at P.endDepth, we need to acquire RF data for the
% longer path lengths of elements that are not directly over the reconstruction pixel. One can compute the worst case
% longest path length for the L11-5v as the square root of the sum of the squares of P.endDepth and the transducer
% aperture, like the maxAcqLength calculation list above.
%
% The Receive.TGC attribute is the number of a Time Gain Control waveform defined in a TGC structure (see the previous
% object defined). This curve defines how the receiver gain increases with time over the depth of the scan.
%
% The next three attributes, 'bufnum', 'framenum', and 'acqnum' define where the acquisition data for a Receive goes in
% the host memory RcvBuffer. Recall that the RcvBuffer (in this case, RcvBuffer number one) was defined with multiple
% frames (Fig. 2.6.1.1 in the programming manual). Within each frame, each column contains all the acquisition data for
% a corresponding channel, in separate acqnum segments, which are packed consecutively down the column. Note that the
% acquisition data are first stored in local memory on the VDAS modules, so these attributes specify the target location
% for the acquisition data when they are transferred to host memory.
%
% The Receive.sampleMode attribute is used to specify the method of sampling for receive data stored in the RcvBuffer.
% Here, it is set to 'NS200BW', which will translate to a sample rate for the RcvBuffer of four times Trans.frequency
% (set to the nearest realizable sample rate based on the 250MHz master clock).  To illustrate how these different
% sampling schemes are implemented, the samples for a cosine wave at 90% of the realizable center frequency (effective
% center frequency or demodFrequency) are shown in Fig. 3.3.1.1 in the programming manual.
%
% The Receive.mode attribute specifies how a receive channel is to store its RF data. For Receive.mode = 0, the RF data
% replaces data in the local memory. For Receive.mode = 1, the acquired data are added to the data already in local
% memory. This mode allows accumulating a number of acquisitions in local memory before transferring the accumulated sum
% to the RcvBuffer (for usage, see coding examples section).
%
% The Receive.LowPassCoef and Receive.InputFilter attributes are missing in our specification, meaning that we will use
% the default values that the system provides. Please use filterTool to visualize the default filter settings and design
% your own filter if needed. Just type filterTool in the command window or check the tutorial video on the website:
% http://verasonics.com/tools-demonstration/
%
% The last attribute, 'callMediaFunction', is used only in simulate mode. When set to 1 for a Receive structure, the
% Media.function name provided in the Media structure is called before acquiring simulated data. This allows moving the
% media points between frames so that the simulation can process moving targets.
%
% The following loop sets the attributes of individual Receive structures:

for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end
end

% In this loop, which modifies the receive structure associated with each acquisition frame, the parameter that is being
% set is: Receive.framenum. The frame number to be used in the RcvBuffer for the acquisition data.
% How to visualize the memory allocation and check possible index error? --> EventAnalysisTool

%% 10. Define the reconstruction structure, Recon. (Section 3.4.1 in Programming Manual)
% With the pixels that we want to process defined by the PData.Region structures, we now need to specify the attributes
% of the reconstruction process, using the Recon and ReconInfo structures. The Recon structure will define the general
% characteristics of the reconstruction for a PData object, including the source and destination buffers, and which
% ReconInfo structures to apply to the reconstruction process. The input data for a Recon structure is always RF data in
% a RcvBuffer, and output data is directed to an InterBuffer for complex signal data and/or an ImageBuffer for magnitude
% data. There can be multiple parts to a reconstruction; for example, in a synthetic aperture reconstruction that
% combines the complex signal data from separate acquisitions. The multiple parts of a reconstruction are defined using
% ReconInfo structures, which will be described below. Since we will need an Image buffer for the output of the
% reconstruction, we must first define the buffer with a Resource specification (should be defined in
% Resource.ImageBuffer already)

Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'RINums',(1:P.numRays));

% Recon.senscutoff is a value from 0.0 to 1.0 that sets the threshold of sensitivity below which an element is excluded
% from the reconstruction summation for a given pixel. A typical value would be 0.6, corresponding to a 4.44 dB loss in
% signal. The Trans.ElementSens function is used to determine an element's relative sensitivity to the pixel point being
% reconstructed if the sensitivity is above the senscutoff value, the element's signal is included in
% reconstruction, otherwise not.
%
% Recon.pdatanum specifies the number of the PData structure that defines the pixel locations of the reconstructed data.
%
% Recon.rcvBufFrame is an optional attribute that when provided, over-rides the frame number specified by the ReconInfo
% structures. If set to a non-zero, positive value, the value is the frame number that will be used for reconstruction.
% If set to -1, and running with the Vantage hardware, the processing determines the last acquisition frame transferred
% into the RcvBuffer used by this Recon structure, and processes this frame.
%
% Recon.IntBufDest and Recon.ImgBufDest each specify the destination buffer and frame that will receive the
% reconstructed output, using a two element array. (The RcvBuffer used for source RF data is provided in the Receive
% object referenced in the ReconInfo objects.) Depending on the reconstruction mode specified in the ReconInfo object,
% one or both buffer destinations may be required. For the frame number value in Recon.ImgBufDest, a value of -1 can be
% used, indicating that the next available frame number should be used for output.
%
% The Recon.RINums attribute is a row vector that specifies the ReconInfo structure indices associated with this
% reconstruction. Each ReconInfo object contains information on how to reconstruct pixels for a specific pixel data
% region.

%% 11. Define the ReconInfo structures. (Section 3.4.2 in Programming Manual)
% The ReconInfo structure(s) define the details of each step in a reconstruction process. Each ReconInfo structure's
% index is referenced in a Recon structure and we must define a ReconInfo structure for each reference. Multiple Recon
% structures must not reference the same ReconInfo structures. The ReconInfo structures reference the transmit and
% receive structures used to acquire the acquisition data, and specify the region of the PData array to reconstruct.
% They also define a mode of reconstruction, that determines whether to output complex signal data or magnitude data,
% and whether to replace the data in the output buffer or accumulate with data already in the buffer.

ReconInfo = repmat(struct('mode', 'accumIQ', ...
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 0, ...
                   'threadSync', 1), 1, P.numRays);

% - Set specific ReconInfo attributes.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

% The ReconInfo.mode attribute defines the output of the reconstruction and how to process it into the output buffer.
% See section 3.4.2 for a complete definition of the different modes. Here, multiple acquisitions are used for one
% frame, so the IQ data of each acquisition will be 'accumulated' and the mode will be 'accumIQ'.
%
% The Pre attribute is used to specify a function that executes prior to the reconstruction action, while the Post
% attribute specifies a function that executes afterwards. The 'clearInterBuf' and 'clearImageBuf' Pre functions clear
% the entire PData pixel region for the destination frame. These functions are required if the reconstruction output
% goes only to a Region that is a subset of the pixels in the entire PData array, since pixels not written to may
% contain data from previous reconstructions. If the reconstruction output goes to the entire PData array, the
% 'replaceIQ' or 'replaceIntensity' reconstruction modes can be used instead.
%
% The Post functions, 'IQ2IntensityImageBuf', 'IQ2IntensityImageBufAdd' and 'IQ2IntensityImageBufMul' are used when
% multiple IQ reconstructions over partial Regions of the PData array are combined in the InterBuffer, as in the case of
% multiple overlapping wide beams.
%
% ReconInfo.txnum specifies the index of the TX object used for transmit. This object is needed to determine the origin
% and characteristics of the transmit beam.
%
% ReconInfo.rcvnum specifies the index of the Receive object used for acquisition. The Receive object contains the
% information on how the RF data was acquired, as well as the storage location of the data. The rcvnum attribute is set
% to 1, but in the Recon structure we have specified rcvBufFrame = -1, meaning that the 1 will get replaced during run
% time with the Receive index of the most recent frame acquired.
%
% ReconInfo.scaleFactor can be provided to scale the output of any reconstruction mode, since it is applied to the I,Q
% reconstruction output before any additional output processing (such as computing the magnitude or accumulating to an
% InterBuffer). It's main use is for weighting multiple synthetic aperture acquisitions which are to be combined in an
% InterBuffer.
%
% ReconInfo.regionnum is the number of the region that is to be reconstructed with this ReconInfo object. Each ReconInfo
% object specifies the reconstruction parameters for one and only one region of the PData structure specified in
% Recon.pdatanum.
%
% The 'threadSync' attribute, when set to 1 (true), specifies that the reconstruction threads should be synchronized
% after the reconstruction of the region specified in the ReconInfo. Thread synchronization is required when a
% reconstruction has multiple ReconInfo structures with unique regionnums to be processed and the pixels in the
% ReconInfo region overlap with pixels in one or more ReconInfo regions yet to be processed. Without thread
% synchronization, different threads might be processing different regions with shared pixels at the same time, writing
% to the same pixels in a race condition.
%
% The rule of thumb:
% 1) One Recon structure is required per Image frame (from one or more acquisitions)
% 2) One ReconInfo structure is required per acquisition (one Tx/Rx event)
% 3) The set of ReconInfo objects should be unique for each Recon
%
% How to inspect my Recon/ReconInfo setting? --> Run EventAnalsysTool and select the Recon number and the tables will be
% displayed

%% 12. Define the Process structure(s). (Section 3.5 in Programming Manual)
% After reconstruction, our echo intensity output data will be found in an ImageBuffer frame, at the same number of
% points (pixels) as the PData array. To put a high quality image up on the display screen of the computer, we typically
% must perform additional image processing operations, such as scaling, compression, reject, and persistence processing.
% Image processing flow is shown in Fig. 3.5.1.1 in the programming manual. Processing operations in the Verasonics are
% defined by a Process object, which defines a classname for the type of processing, a method to use, and list of
% parameters that help define additional elements of the processing. For processing image data, we use the classname
% 'Image', and the method 'imageDisplay'. The main attributes of this Process structure are shown below.

pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...            % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% The first two parameters, imgbufnum and framenum defined the frame in the ImageBuffer to be use as input. Here again a
% value of -1 for the frame number indicates that the most recent frame reconstructed should be processed.
%
% The next parameter, pdatanum, is the index of the PData structure that was used for reconstruction. This structure is
% needed to get the image size (rows and columns) as well as the pixel delta values.
%
% The pgain parameter (short for processing gain) specifies a gain factor to apply to the reconstruction intensity
% values as they are processed. A value of 1.0 leaves the data unchanged. The reconstruction intensity values are
% normalized by the number of channels contributing to the pixel reconstruction, but if some or all of the channels have
% low amplitude signals, the reconstruction values will need to be amplified to produce an adequate display intensity
% value.
%
% The reject value cuts off low level intensities before mapping intensities to the display. The value is the
% percentage of the lower quartile of the intensity range to reject.
%
% Set the persistence method to 'simple', and the level to 20, which means that the intensity value output consists of
% 20% of the previous frame's intensity plus 80% of the new intensity value.
%
% The interpMethod sets the interpolation method. In this case, '4pt' means to use a 4 point bi-linear interpolation
% method. This is currently the only method supported.
%
% The next attribute, grainRemoval, can be used to enable a 3x3 matrix filter that aims to eliminate single pixels that
% differ significantly from their neighbors. This filter has three settings - 'low', 'medium' and 'high' which offer
% varying amounts of filtering. An additional 3x3 matrix filter, processMethod, can also be invoked, which aims at
% reducing variation in line structures detected within the filter kernel.
%
% The averageMethod attribute allows implementing a running average of the current and previous frames. The available
% choices are 'none', 'runAverage2', and 'runAverage3', of which the latter two select the running average of the
% current frame with the previous one or two frames. This allows combining frames with different processing to achieve
% special effects, such as spatial compounding.
%
% Set the 'compressMethod' and 'compressFactor'. The 'compressMethod' options are 'power' or 'log'. The 'power' option
% raises the intensity data to the power n, where n is a fraction set by the 'compressFactor' value (40 equates to 0.5).
% Additional compression or expansion can accomplished by modifying the display windows colormap.
%
% The mappingMethod parameter specifies the portion of the colormap to be used by the intensity values. The choices are
% lowerHalf, upperHalf, and full. The lowerHalf and upperHalf choices are used when the color palette needs to be split
% between two ultrasound modes, such as 2D and color Doppler. For normal 2D processing, we use the full setting.
%
% Setting the 'display' attribute to 1 means that we want to display the image after processing. There are some
% conditions where we don't want the display to appear until additional processing has been performed.  The
% 'displayWindow' attribute sets the index of the Resource.DisplayWindow that we want the image to appear in.
%
% The other two objects supported in the Process.classname are Doppler and External. Please refer to other example
% scripts and programming manual regarding the use of them.

%% 13. Define the SeqControl and Event structures. (Section 3.6 and 3.7 in Programming Manual)
% We are now at the point where we are ready to define the individual sequence events for our sequence. Each event can
% reference a transmit action, a receive action, a reconstruction action, and a processing action. Finally, a SeqControl
% action can be set that controls factors such as DMA transfers and event flow control. For actions that are not needed,
% a zero index is set, which indicates no action for that item. In the case of tx and rcv actions, these must both be
% set or only the tx index set (for a transmit only event). A receive only Event can have a tx reference to a TX
% structure with all transmitters disabled (TX.Apod = zeros(128)).
%
% It is often easier to define the SeqControl structures at the point where they are needed in the sequence of Events.
% The exception to this approach would be 'transferToHost' which needs be unique in the event sequence.
%
% It is important to keep in mind that the hardware and software sequencer are running with different speed. Some
% actions are hardware only, like tx/rcv, and some actions are software only, like recon and process. Please refer to
% the Table 8.1 in the programming tutorial about Hardware and Software Sequencer Actions, and Section 3.6 in the
% programming manual about other useful seqControl commands, like trigger and sync.

% The 'jump' command will have the event sequence jump to the event defined in the argument.
SeqControl(1).command = 'jump';
SeqControl(1).argument = 1;

% The 'timeToNextAcq' command is used to set a precise time between acquisition events. It can also be used to set the
% duration between transmit only Events (Event.rcv = 0). The time duration is set in the argument field in microseconds
% (10 - 4190000).
SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 220;  % 220 usec
SeqControl(3).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(3).argument = 30000;  % 30000 usec = 30msec time between frames

% When 'returnToMatlab' command is encountered, the processing function, runAcq, suspends execution and returns to the
% Matlab environment.  The next time runAcq is called (with no startEvent specified), execution resumes at the next
% Event in the sequence. The callback of any UI control will not be executed until the sequencer reaches
% 'returnToMatlab'.
SeqControl(4).command = 'returnToMatlab';

nsc = 5; % Define a counter variable, nsc, to keep track of new SeqControl indices.

% Specify Event structure arrays.
n = 1;

% The acquisition and reconstruction for all the frames in the RcvBuffer are defined in a loop. This makes it easy to
% change the number of frames acquired in the RF cineloop by simply changing the Resource.RcvBuffer.numFrames attribute.
for i = 1:Resource.RcvBuffer(1).numFrames

    % With this script, multiple acquisitions are required for one frame, so a for loop is needed for each rcv frame.
    % This for loop is used for data acquisition only (hardware sequencer)
    for j = 1:P.numRays

        % Event.info only provides information as reference
        Event(n).info = ['Acq ',num2str(j),' for frame ', num2str(i)];

        % Event.tx is the index of TX. For instance, if Event(1).tx = 1, then TX(1) will be used in Event(1)
        Event(n).tx = j;

        % Event.rcv is the index of Receive. For instance, if Event(1).rcv = 1, then Receive(1) will be used in Event(1)
        % Each acquisition requires a unique Receive definition
        Event(n).rcv = P.numRays*(i-1)+j;

        % Recon and Process are not required for data acquisition
        Event(n).recon = 0;
        Event(n).process = 0;

        % SeqControl(2) defines the interval between acquisition
        Event(n).seqControl = 2;
        n = n+1;
    end

    % The last acquisition Event's seqControl needs to be modified for data transfer and frame rate control
    Event(n-1).seqControl = [3,nsc];

    % Again, 'transferToHost' can't be predefined! In addition, use the EventAnalysisTool and select the seqControl
    % number for transferToHost, the estimated DMA information will be displayed
      SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
      nsc = nsc+1;

    % Specify a reconstruction after each acquired frame. The hardware sequencer ignores these events and continues to
    % acquire new frames while the software is doing reconstruction processing. The software processing ignores tx and
    % rcv actions (unless in simulate mode), and with Recon.rcvBufFrame = -1, will try and process the most recently
    % transferred frame of data. If the RcvDataLoop control is toggled, the reconstruction will process the next frame
    % in the RcvBuffer from the last one processed, and will sequentially process every frame in the buffer.
    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;

    % Since the RF cineloop might be quite long, resulting in a lot of frames reconstructed before reaching the end of the
    % sequence, we insert a 'returnToMatlab' at every 5th frame reconstructed. This allows for more rapid response to GUI
    % controls, which can only happen while back in the Matlab environment.
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame reconstructed
        Event(n).seqControl = 4;
    end
    n = n+1;
end

% Jump back to the first event will 'returnToMatlab' automatically. Please refer to the programming manual for more
% information
Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

%% 14. User specified UI Controls (Section 4.1 in Programming Manual)
% In the Verasonics software, we have provided a way to specify user GUI controls and callback functions within the
% SetUp scripts, thus allowing the controls to be kept with the script that uses them. This mechanism incidentally also
% provides a way to execute Matlab commands just before our script runs, which can be useful for scripts that need to
% initialize something just before running.
%
% To facilitate creating user controls on the Verasonics GUI, newer software releases provide some built-in UI functions
% that can be positioned at certain predefined locations. The currently available VsXXX controls that can be selected
% are: VsSlider, VsPushButton, VsToggleButton, and VsButtonGroup. Please refer to the section 4.1.1 for more details.

% User specified UI Control Elements - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};

% To make it easier to create and interpret long callback functions, a utility function is provided that will translate
% text between two delimiters into cell array strings, called 'text2cell' whose calling inputs are the filename to read
% from (optional; if missing defaults to calling script file) and the delimiter string. This allows writing normal
% Matlab code at the end of the Setup file that can be automatically interpreted into the UI.Callback cell array.
%
% For instance, here, text2Cell will look for the delimiter, %SensCutOffCallback, and convert all codes in between into cells.
% VSX will then create a associated callback function named UserB7Callback located in the tempdir.
UI(1).Callback = text2cell('%SensCutOffCallback');

% - Range Change
MinMaxVal = [64,300,P.endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% - Peak CutOff
UI(3).Control = {'UserB2','Style','VsSlider','Label','Peak Cutoff',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakCutOff],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(3).Callback = text2cell('%PeakCutOffCallback');

% - Max. Burst Length
UI(4).Control = {'UserB1','Style','VsSlider','Label','Max. BL',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakBLMax],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(4).Callback = text2cell('%MaxBLCallback');

%% 15. frame rate estimation and save into a mat file for VSX
% The frameRateFactor variable is used to account for the fact that, in the Event sequence, a script may process more
% than one frame before it returns to matlab (typically when it jumps back to event #1 to start over). For example, if
% the Event Sequence processed four frames between two returnToMatlab events, the tElapsed measured by VSX would
% represent the total time used to process all four frames, so the frame rate implied by the moving average
% sequencePeriod value would have to be multiplied by four, to provide an estimate of the actual frame rate of the
% software process, not the hardware frame rate.
frameRateFactor = 3;

% Save all structures into a .mat file. If the variable filename is defined with the name of the Setup script matfile,
% VSX will use it, and not ask for the "Name of .mat file to process". For instance, the user can uncomment the
% following commands to skip the query of "Name of .mat file to process" after running VSX.
% filename = 'L11-5vWideBeam';
% save(['MatFiles/',filename]);
save('MatFiles/L11-5vWideBeam');

% The EventAnalysisTool will use the initializedOnly feature to find possible errors and display many useful
% information. It's recommended to run EventAnalysiTool before running the VSX, so we list it here to remind the user to
% use it.
EventAnalysisTool
return

%% 16. Callback routines to be converted by text2cell function.
% Some papameters will be changed by the callback function and the software uses the Control structure to update parameters
% Please refer to section 4.2 in the programming manual for details.

%SensCutOffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutOffCallback

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);

PData = evalin('base','PData');
PData.Size(1) = ceil((P.endDepth-P.startDepth)/PData.PDelta(3));
PData.Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing*(P.numTx-12),...
                    'height',P.endDepth-P.startDepth)),1,P.numRays);
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end
PData.Region = computeRegions(PData);
assignin('base','PData',PData);
% Update TXPD data of TX structures.
TX = evalin('base','TX');
for i = 1:size(TX,2)
    TX(i).TXPD = computeTXPD(TX(i),PData);
end
assignin('base','TX',TX);
% Update Receive structures.
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData.Size(1)*PData.PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','TX','Receive','Recon','TGC','DisplayWindow'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%PeakCutOffCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakCutOff = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%PeakCutOffCallback

%MaxBLCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakBLMax = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%MaxBLCallback
