package nanoj.squirrel.java.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.gui._BaseDialog_;
import nanoj.kernels.Kernel_VoronoiImage;
import nanoj.squirrel.java.gui.tools.SetMaximumStackSize_;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.awt.*;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.*;
import static nanoj.core.java.array.ArrayMath.calculateMSE;
import static nanoj.core.java.array.ArrayMath.calculatePPMCC;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.MAX_FITTING;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.getShiftFromCrossCorrelationPeak;
import static nanoj.core.java.image.transform.CrossCorrelationMap.calculateCrossCorrelationMap;
import static nanoj.core.java.image.transform.CrossCorrelationMap.cropCCM;
import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_Errors;
import static nanoj.squirrel.java.gui.ImagesHelper.magnify;

/**
 * Created by sculley on 09/01/2019.
 */
public class ErrorMapV2Blocked_ extends _BaseDialog_ {
    private static final float ROOT2 = (float) Math.sqrt(2);
    private static final String OVERBLUR_MESSAGE = "RSF is unrealistic size... Optimising to PSF estimate";

    String titleRefImage = "", titleRSFImage = "", titleSRImage = "", noRSFString = "-- RSF unknown, estimate via optimisation --";
    boolean showAdvancedSettings, _showAdvancedSettings = false;
    boolean doSlidingWindow;

    protected ErrorMap_ExtraSettings_ errorMap_ExtraSettings = new ErrorMap_ExtraSettings_();
    protected ErrorMap_ExtraSettings_ _errorMap_ExtraSettings;

    protected SetMaximumStackSize_ maximumStackSize = new SetMaximumStackSize_();
    ArrayList<String> titles = new ArrayList<String>();
    String[] imageTitles;

    int maxSigma;
    boolean showPlot, framePurge, borderControl;
    boolean doRegistration;
    int maxExpectedMisalignment;
    int maxSRStackSize;
    boolean showIntensityNormalised, showConvolved, showRSF, showPositiveNegative;

    //Convolve convolve = new Convolve();

    double sigmaGuess;
    ImagePlus impRef, impSR, impRSF;
    ImageStack imsRef, imsSR, imsRSF;
    int nSlicesRef, nSlicesSR, nSlicesRSF;
    int w_SR, h_SR, w_Ref, h_Ref;
    int magnification, magnification2;
    boolean noCrop = true;
    private boolean borderCrop = false;
    private boolean registrationCrop = false;
    private int maxMag, blocksPerXAxis, blocksPerYAxis;

    private final static Kernel_VoronoiImage VI = new Kernel_VoronoiImage();

    DecimalFormat df = new DecimalFormat("00.00");


    @Override
    public boolean beforeSetupDialog(String s) {
        autoOpenImp = false;
        useSettingsObserver = true;
        maxSRStackSize = (int) maximumStackSize.getPrefs("maxSRImages",50);
        int nImages = WindowManager.getImageCount();
        if(nImages<2){
            log.error("At least 2 images are required to run this code!");
            return false;
        }

        imageTitles = WindowManager.getImageTitles();
        for(int n=0; n<nImages; n++){
            ImagePlus thisImp = WindowManager.getImage(imageTitles[n]);
            if(thisImp.getStackSize()<=maxSRStackSize){
                titles.add(thisImp.getTitle());
            }
        }

        imageTitles = new String[titles.size()];
        for(int n=0; n<titles.size(); n++){
            imageTitles[n] = titles.get(n);
        }

        return true;
    }

    @Override
    public void setupDialog() {
        // create string array with 'estimate RSF' as option
        String[] imageTitlesRSF = new String[imageTitles.length + 1];
        imageTitlesRSF[0] = noRSFString;
        for(int i=0; i<imageTitles.length;i++) imageTitlesRSF[i+1] = imageTitles[i];

        gd = new NonBlockingGenericDialog("Calculate Error Map, RSE and RSP");
        gd.addHelp("https://bitbucket.org/rhenriqueslab/nanoj-squirrel");

        if (titleRefImage.equals("")) {
            titleRefImage = getPrefs("titleRefImage", "");
            titleSRImage = getPrefs("titleSRImage", "");
            titleRSFImage = getPrefs("titleRSFImage", "");
        }

        if (!titleRefImage.equals("") && contains(titleRefImage, imageTitles))
            gd.addChoice("Reference Image", imageTitles, titleRefImage);
        else
            gd.addChoice("Reference Image", imageTitles, imageTitles[0]);

        if (!titleSRImage.equals("") && contains(titleSRImage, imageTitles))
            gd.addChoice("Super-Resolution Reconstruction(s)", imageTitles, titleSRImage);
        else
            gd.addChoice("Super-Resolution Reconstruction(s)", imageTitles, imageTitles[0]);

        if (!titleRSFImage.equals("") && contains(titleRSFImage, imageTitles))
            gd.addChoice("RSF Estimate Image(s)", imageTitlesRSF, titleRSFImage);
        else
            gd.addChoice("RSF Estimate Image(s)", imageTitlesRSF, imageTitlesRSF[0]);

        gd.addNumericField("Max. Mag. in Optimization (default: 5)", getPrefs("maxMag", 5), 0);
        gd.addMessage("The higher the maximum magnification, the more precise (marginally) but slower the algorithm will be...", new Font("Arial", Font.ITALIC, 12));

        gd.addCheckbox("Show_Advanced_Settings", false);
        gd.addCheckbox("Do sliding window", getPrefs("doSlidingWindow", false));
    }

    @Override
    public boolean loadSettings() {
        titleRefImage = gd.getNextChoice();
        titleSRImage = gd.getNextChoice();
        titleRSFImage = gd.getNextChoice();
        maxMag = (int) gd.getNextNumber();

        showAdvancedSettings = gd.getNextBoolean();
        if (!_showAdvancedSettings && showAdvancedSettings && _errorMap_ExtraSettings == null) {
            _errorMap_ExtraSettings = new ErrorMap_ExtraSettings_();
            _errorMap_ExtraSettings.start();
            _showAdvancedSettings = true;
        }
        else if (_errorMap_ExtraSettings != null && _errorMap_ExtraSettings.gd != null) {
            _errorMap_ExtraSettings.gd.dispose();
            _errorMap_ExtraSettings = null;
        }
        if (showAdvancedSettings) _showAdvancedSettings = true;
        else _showAdvancedSettings = false;

        doSlidingWindow = gd.getNextBoolean();

        setPrefs("titleRefImage", titleRefImage);
        setPrefs("titleSRImage", titleSRImage);
        setPrefs("titleRSFImage", titleRSFImage);
        setPrefs("maxMag", maxMag);
        setPrefs("doSlidingWindow", doSlidingWindow);
        if (maxMag < 1) return false;

        prefs.savePreferences();

        impRef = WindowManager.getImage(titleRefImage);
        impSR = WindowManager.getImage(titleSRImage);
        if(titleRSFImage!=noRSFString) impRSF = WindowManager.getImage(titleRSFImage);
        else sigmaGuess = 0;

        imsRef = impRef.getImageStack().convertToFloat();
        imsSR = impSR.getImageStack().convertToFloat();
        if (impRSF!=null) imsRSF = impRSF.getImageStack().convertToFloat();

        w_SR = imsSR.getWidth();
        h_SR = imsSR.getHeight();

        magnification = w_SR / imsRef.getWidth();
        magnification2 = magnification * magnification;

        nSlicesSR = impSR.getStackSize();
        nSlicesRef = impRef.getStackSize();
        if(impRSF!=null) {nSlicesRSF = impRSF.getStackSize();}

        // check if SR size is integer multiple
        double _magnification = ((double) w_SR) / imsRef.getWidth();
        if (magnification != _magnification) {
            log.error("SR image size needs to be an integer multiple of the reference image size!");
            return false;
        }

        // check number of frames in Reference
        if (nSlicesRef>1) {
            log.error("Not supported for multiple reference frames, please check reference image selection!");
            return false;
        }

        // check number of frames in RSF
        if (impRSF!=null && imsRSF.getSize() != 1 && nSlicesRSF != nSlicesSR){
            IJ.error("Number of slices in RSF image(s) needs to be 1 or the same number of SR reconstructions.");
            return false;
        }

        return true;
    }

    // Import settings from advanced dialog
    public boolean loadSettingsFromPrefs() {

        maxSigma = errorMap_ExtraSettings.getPrefs("maxSigma", 200);
        showPlot = errorMap_ExtraSettings.getPrefs("showPlot", false);
        framePurge = errorMap_ExtraSettings.getPrefs("framePurge", false);
        borderControl = errorMap_ExtraSettings.getPrefs("borderControl", true);

        doRegistration = errorMap_ExtraSettings.getPrefs("doRegistration", true);
        maxExpectedMisalignment = errorMap_ExtraSettings.getPrefs("maxExpectedMisalignment", 0);

        showIntensityNormalised = errorMap_ExtraSettings.getPrefs("showIntensityNormalised", true);
        showConvolved = errorMap_ExtraSettings.getPrefs("showConvolved", true);
        showRSF = errorMap_ExtraSettings.getPrefs("showRSF", false);
        showPositiveNegative = errorMap_ExtraSettings.getPrefs("showPositiveNegative", false);

        prefs.savePreferences();

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        long startTime = System.currentTimeMillis();


        // Check and purge homogeneous/empty slices
        if(framePurge) {
            ArrayList<Integer> slicesToDelete = new ArrayList<Integer>();

            for (int n = 1; n <= nSlicesSR; n++) {

                log.status("Checking stack, frame " + n);
                log.progress(n, nSlicesSR);

                ImageProcessor ip = imsSR.getProcessor(n);
                FloatProcessor fp = imsSR.getProcessor(n).convertToFloatProcessor();
                float[] pixels = (float[]) fp.getPixels();

                boolean deleteFlag = true;

                float firstVal = pixels[0];
                for (int p = 1; p < pixels.length; p++) {
                    if (pixels[p] != firstVal) {
                        deleteFlag = false;
                        break;
                    }
                }

                if (deleteFlag) slicesToDelete.add(n);
            }

            log.msg("Empty Frame Check:");

            if (!slicesToDelete.isEmpty()) {
                int nSlicesToDelete = slicesToDelete.size();
                int[] slicesToDeleteArray = new int[nSlicesToDelete];
                for (int i = 0; i < nSlicesToDelete; i++) {
                    slicesToDeleteArray[i] = slicesToDelete.get(i);
                }

                log.msg("\t Frames purged: " + Arrays.toString(slicesToDeleteArray));

                for (int i = nSlicesToDelete - 1; i >= 0; i--) {
                    imsSR.deleteSlice(slicesToDeleteArray[i]);
                }

                // update number of slices in SR stack
                nSlicesSR = imsSR.getSize();
            } else {
                log.msg("\t All frames contain data!");
            }
        }

        // Crop black borders from images
        if(borderControl){
            checkAndCropBorders();

            // update width and height
            w_SR = imsSR.getWidth();
            h_SR = imsSR.getHeight();
        }

        // Run cross-correlation realignment and crop realignment borders
        if(doRegistration) {
            log.status("Registering and cropping images...");
            realignImages();
        }

        // update width and height
        w_SR = imsSR.getWidth();
        h_SR = imsSR.getHeight();
        int nPixelsSR = w_SR*h_SR;

        // Set up reference float processors
        FloatProcessor fpRef = imsRef.getProcessor(1).convertToFloatProcessor();
        fpRef.resetRoi();

        w_Ref = fpRef.getWidth();
        h_Ref = fpRef.getHeight();

        float[] pixelsRef = (float[]) fpRef.getPixels();

        FloatProcessor fpRefScaledToSR = (FloatProcessor) fpRef.duplicate();
        fpRefScaledToSR.setInterpolationMethod(ImageProcessor.BICUBIC);
        fpRefScaledToSR = (FloatProcessor) fpRefScaledToSR.resize(w_SR, h_SR);
        float[] pixelsRefScaledToSR = (float[]) fpRefScaledToSR.getPixels();

        // output stacks
        ImageStack imsSRIntensityScaled = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsSRConvolved = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsEMap = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsAlphaMaps = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsBetaMaps = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsSigmaMaps = new ImageStack(w_SR, h_SR, nSlicesSR);

        ResultsTable rt = new ResultsTable();

        // handle sigma boundary

        float maxSigmaBoundary;

        // Adjust maxSigmaBoundary if there is calibration data
        String pixelUnitRef = impRef.getCalibration().getUnit();
        double pixelSizeNm;
        if(pixelUnitRef.equals("nm")){
            log.msg("Detected nm!");
            pixelSizeNm = impRef.getCalibration().pixelWidth;
            maxSigmaBoundary = (float) (maxSigma/pixelSizeNm)*magnification;
            log.msg("Maximum sensible RSF sigma = "+maxSigma+"nm ("+maxSigmaBoundary+" SR pixels)");
        }
        else if(pixelUnitRef.equals("micron")|| pixelUnitRef.equals("microns") || pixelUnitRef.equals("µm")){
            log.msg("Detected µm!");
            pixelSizeNm = impRef.getCalibration().pixelWidth * 1000;
            maxSigmaBoundary = (float) (maxSigma/pixelSizeNm)*magnification;
            log.msg("Maximum sensible RSF sigma = "+maxSigma+"nm ("+maxSigmaBoundary+" SR pixels)");
        }
        else{
            //assign arbitrary pixel size of 100nm
            pixelSizeNm = 100;
            maxSigmaBoundary = (float) (maxSigma/pixelSizeNm)*magnification;
            log.msg("Maximum sensible RSF sigma = "+maxSigmaBoundary+" SR pixels");
        }


        // blocking
        int blockSize = (int) (10*maxSigmaBoundary);
        while(blockSize%magnification!=0) blockSize++;
        log.msg("Block size is "+blockSize);

        blocksPerXAxis = (int) floor(w_SR/blockSize);
        blocksPerYAxis = (int) floor(h_SR/blockSize);

        double inc = 1;
        if(doSlidingWindow) inc = 0.5;
        double totalBlocks = (blocksPerYAxis/inc)*(blocksPerXAxis*inc);


        long loopStart = System.nanoTime();

        for(int s=0; s<nSlicesSR; s++) {

            log.msg("-----------------------------------");
            log.msg("Processing super-resolution frame "+(s+1));
            log.msg("-----------------------------------");

            FloatProcessor fpSR = imsSR.getProcessor(s+1).convertToFloatProcessor();
            float[] pixelsSR = (float[]) fpSR.getPixels();

            // create arraylists for parameter maps
            ArrayList<Float> xPositionsList = new ArrayList<Float>();
            ArrayList<Float> yPositionsList = new ArrayList<Float>();
            ArrayList<Float> alphaList = new ArrayList<Float>();
            ArrayList<Float> betaList = new ArrayList<Float>();
            ArrayList<Float> sigmaList = new ArrayList<Float>();


            for (double nYB = 0; nYB < blocksPerYAxis; nYB+=inc) {
                for (double nXB = 0; nXB < blocksPerXAxis; nXB+=inc) {

                    double thisBlock = nYB*blocksPerXAxis + nXB;
                    log.progress(thisBlock/totalBlocks);

                    boolean localOverblurFlag = false;
                    FloatProcessor fpSRBlock = getBlockFp(fpSR, nYB, nXB);
                    float[] pixelsRefBlock = getBlockPixels(fpRef, nYB, nXB);

                    int blockWidthSR = fpSRBlock.getWidth();
                    int blockHeightSR = fpSRBlock.getHeight();
                    int xStartSR = (int) nXB*blockWidthSR;
                    int yStartSR = (int) nYB*blockHeightSR;
                    int nPixelsSRBlock = blockWidthSR*blockHeightSR;
                    int blockWidthRef = blockWidthSR/magnification;
                    int blockHeightRef = blockHeightSR/magnification;

                    if(nPixelsSRBlock!=(pixelsRefBlock.length*magnification2)) continue;

                    float[] ones = new float[nPixelsSRBlock];
                    for(int i=0; i<nPixelsSRBlock; i++){ones[i] = 1;}


                    // UNIVARIATE OPTIMIZER - LINEAR MATCHING
                    /// setup optimizer
                    sigmaOptimiseFunction f =  new ErrorMapV2Blocked_.sigmaOptimiseFunction(fpSRBlock, pixelsRefBlock, ones);
                    UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);

                    /// run optimizer
                    UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                            new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, blockWidthSR/2)); //limit to block width
                    float sigma_linear = (float) result.getPoint();

                    // GET ALPHA AND BETA
                    FloatProcessor blurredFp = (FloatProcessor) fpSRBlock.duplicate();
                    FloatProcessor blurredOnes = new FloatProcessor(blockWidthSR, blockHeightSR, ones);
                    blurredFp.blurGaussian(sigma_linear);
                    blurredOnes.blurGaussian(sigma_linear);

                    blurredFp = (FloatProcessor) blurredFp.resize(blockWidthRef, blockHeightRef);
                    blurredOnes = (FloatProcessor) blurredOnes.resize(blockWidthRef, blockHeightRef);

                    float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRefBlock, (float[]) blurredOnes.getPixels());
                    float alpha = aB[0];
                    float beta = aB[1];

                    // check if sigma hit the boundary

                    float alphaBoundary = alpha, betaBoundary = beta, sigmaBoundary = sigma_linear;
                    if(sigma_linear>maxSigmaBoundary){
                        sigmaBoundary = maxSigmaBoundary;
                        localOverblurFlag = true;
                        //log.msg(OVERBLUR_MESSAGE);
                        // calculate alpha and beta with maxSigmaBoundary as blur
                        blurredFp = (FloatProcessor) fpSRBlock.duplicate();
                        blurredOnes = new FloatProcessor(blockWidthSR, blockHeightSR, ones);
                        blurredFp.blurGaussian(maxSigmaBoundary);
                        blurredOnes.blurGaussian(maxSigmaBoundary);

                        blurredFp = (FloatProcessor) blurredFp.resize(blockWidthRef, blockHeightRef);
                        blurredOnes = (FloatProcessor) blurredOnes.resize(blockWidthRef, blockHeightRef);

                        aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRefBlock, (float[]) blurredOnes.getPixels());
                        alphaBoundary = aB[0];
                        betaBoundary = aB[1];
                    }

                    if(alphaBoundary<0) continue;

                    xPositionsList.add((float) (xStartSR + blockWidthSR/2));
                    yPositionsList.add((float) (yStartSR + blockHeightSR/2));
                    alphaList.add(alphaBoundary);
                    betaList.add(betaBoundary);
                    sigmaList.add(sigmaBoundary);

                }
            }

            // convert arraylists into float arrays and create voronoi maps
            float[] xPositions = arrayListToFloatArray(xPositionsList);
            float[] yPositions = arrayListToFloatArray(yPositionsList);
            float[] alphas = arrayListToFloatArray(alphaList);
            float[] betas = arrayListToFloatArray(betaList);
            float[] sigmas = arrayListToFloatArray(sigmaList);

            FloatProcessor fpAlphaMap = VI.calculateImage(w_SR, h_SR, xPositions, yPositions, alphas);
            FloatProcessor fpBetaMap = VI.calculateImage(w_SR, h_SR, xPositions, yPositions, betas);
            FloatProcessor fpSigmaMap = VI.calculateImage(w_SR, h_SR, xPositions, yPositions, sigmas);

            fpAlphaMap.blurGaussian(blockSize/2);
            fpBetaMap.blurGaussian(blockSize/2);
            fpSigmaMap.blurGaussian(blockSize/2);

            float[] pixelsAlphaMap = (float[]) fpAlphaMap.getPixels();
            float[] pixelsBetaMap = (float[]) fpBetaMap.getPixels();
            float[] pixelsSRIntensityScaled = new float[nPixelsSR];

            for(int i=0; i<nPixelsSR; i++){
                pixelsSRIntensityScaled[i] = pixelsSR[i]*pixelsAlphaMap[i]+pixelsBetaMap[i];
            }
            FloatProcessor fpSRIntensityScaled = new FloatProcessor(w_SR, h_SR, pixelsSRIntensityScaled);
            imsSRIntensityScaled.setProcessor(fpSRIntensityScaled, s+1);

            FloatProcessor fpSRConvolved = new FloatProcessor(w_SR, h_SR);
            // go back through blocks and convolve with local sigma
            for(int nYB=0; nYB<blocksPerYAxis; nYB++){
                for(int nXB=0; nXB<blocksPerXAxis; nXB++){
                    int blockWidth = w_SR/blocksPerXAxis;
                    int blockHeight = h_SR/blocksPerYAxis;

                    int xStartSR = nXB*blockWidth;
                    int yStartSR = nYB*blockHeight;
                    int xCentreSR = xStartSR+blockWidth/2;
                    int yCentreSR = yStartSR+blockHeight/2;

                    FloatProcessor fp = fpSRIntensityScaled.duplicate().convertToFloatProcessor();
                    float thisSigma = fpSigmaMap.getf(xCentreSR, yCentreSR);
                    fp.blurGaussian(thisSigma);

                    fp = getBlockFp(fp, nYB, nXB);

                    fpSRConvolved = setBlock(fpSRConvolved, fp, xStartSR, yStartSR, fp.getWidth(), fp.getHeight());
                }
            }
            imsSRConvolved.setProcessor(fpSRConvolved, s+1);
            imsAlphaMaps.setProcessor(fpAlphaMap, s+1);
            imsBetaMaps.setProcessor(fpBetaMap, s+1);
            imsSigmaMaps.setProcessor(fpSigmaMap, s+1);

            // CALCULATE METRICS AND MAP
            log.status("Calculating similarity...");

            /// metrics
            FloatProcessor fpSRConvolved_RefSize = (FloatProcessor) fpSRConvolved.resize(w_Ref, h_Ref);
            float[] pixelsSRConvolved_RefSize = (float[]) fpSRConvolved_RefSize.getPixels();
            double globalRMSE = sqrt(calculateMSE(pixelsSRConvolved_RefSize, pixelsRef));
            double globalPPMCC = calculatePPMCC(pixelsSRConvolved_RefSize, pixelsRef, true);

            /// error map
            float[] pixelsEMap = new float[nPixelsSR];
            float[] pixelsSRC = (float[]) fpSRConvolved.getPixels();

            float maxRef = -Float.MAX_VALUE;
            for(int p=0; p<nPixelsSR; p++){
                float vRef = pixelsRefScaledToSR[p];
                float vSRC = pixelsSRC[p];

                maxRef = max(maxRef, vRef); //why
                if(showPositiveNegative) pixelsEMap[p] = (vRef - vSRC);
                else pixelsEMap[p] = abs(vRef-vSRC);
            }
            imsEMap.setProcessor(new FloatProcessor(w_SR, h_SR, pixelsEMap), s+1);

            // CALCULATE METRICS AND MAP FOR BOUNDARY PROBLEM CASES
            double globalRMSEBoundary = globalRMSE, globalPPMCCBoundary = globalPPMCC;

            /// put error values into table
            rt.incrementCounter();
            rt.addValue("Frame", s+1);
            rt.addValue("RSP (Resolution Scaled Pearson-Correlation)", globalPPMCC);
            rt.addValue("RSE (Resolution Scaled Error)", globalRMSE);
//            rt.addValue("Overblur warning", String.valueOf(globalOverblurFlag));
//            if(!globalOverblurFlag){
//                rt.addValue("RSP - constrained", "N/A");
//                rt.addValue("RSE - constrained", "N/A");
//            }
//            else{
//                rt.addValue("RSP - constrained", globalPPMCCBoundary);
//                rt.addValue("RSE - constrained", globalRMSEBoundary);
//            }

            if(nSlicesSR>1) {
                double frameTime = ((System.nanoTime() - loopStart) / (s+1)) / 1e9;
                double remainingTime = frameTime * (nSlicesSR - (s+1));
                int _h = (int) (remainingTime / 3600);
                int _m = (int) (((remainingTime % 86400) % 3600) / 60);
                int _s = (int) (((remainingTime % 86400) % 3600) % 60);
                log.msg("Estimated time remaining to complete analysis:");
                log.msg("\t" + String.format("%02d:%02d:%02d", _h, _m, _s));
            }

        }

        /// BUILD RSF STACK
        //TODO: this.

        // Output images

        String borderString = new String();
        if(borderCrop) {borderString = "Border cropped ";}

        String registrationString = new String();
        if(registrationCrop) {
            if (borderCrop) registrationString = "and Registered";
            else registrationString = "Registered";
        }

        if(!noCrop) new ImagePlus(titleRefImage+" - " +borderString+" " +registrationString, imsRef).show();
        if(showIntensityNormalised || !noCrop){
            new ImagePlus(titleSRImage +" - " +borderString +registrationString+" - intensity-normalised", imsSRIntensityScaled).show();
        }

        if(showConvolved){
            new ImagePlus(titleSRImage+" - Convolved with RSF", imsSRConvolved).show();
        }
        //if(showRSF) new ImagePlus("Optimised RSF stack", imsRSF_rescaled).show();

        ImagePlus impEMap = new ImagePlus(titleSRImage+" - Resolution Scaled Error-Map", imsEMap);
        applyLUT_SQUIRREL_Errors(impEMap);
        impEMap.show();


        new ImagePlus("Alpha maps", imsAlphaMaps).show();
        new ImagePlus("Beta maps", imsBetaMaps).show();
        new ImagePlus("Sigma maps", imsSigmaMaps).show();

        rt.show("RSP and RSE values");

        IJ.run("Tile");

        log.msg("SQUIRREL analysis took "+(System.currentTimeMillis()-startTime)/1000+"s");

    }

    //helper function lifted from imageshelper
    protected static boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
    }

    private void checkAndCropBorders() {
        int w_Ref = imsRef.getWidth();
        int h_Ref = imsRef.getHeight();

        int maxTopBorder = 0, maxRightBorder = 0, maxBottomBorder = 0, maxLeftBorder = 0;

        for(int n=1; n<=nSlicesSR; n++) {

            log.status("Checking borders of frame "+n);
            log.progress(n, nSlicesSR);

            FloatProcessor fpSR = imsSR.getProcessor(n).convertToFloatProcessor();

            int topBorderCounter = 0, rightBorderCounter = 0, bottomBorderCounter = 0, leftBorderCounter = 0;
            float topBorderValue = 0, rightBorderValue = 0, bottomBorderValue = 0, leftBorderValue = 0;

            outerloop:
            while(topBorderValue==0){
                for(int x=0; x<w_SR; x++){
                    topBorderValue = fpSR.getPixelValue(x, topBorderCounter);
                    if(topBorderValue>0 && topBorderValue!=Float.NaN) break outerloop;
                }
                topBorderCounter++;
            }

            outerloop:
            while(rightBorderValue==0){
                for(int y=0; y<h_SR; y++){
                    rightBorderValue = fpSR.getPixelValue(w_SR - rightBorderCounter -1, y);
                    if(rightBorderValue>0 && rightBorderValue!=Float.NaN) break outerloop;
                }
                rightBorderCounter++;
            }

            outerloop:
            while(bottomBorderValue==0){
                for(int x=0; x<w_SR; x++){
                    bottomBorderValue = fpSR.getPixelValue(x, h_SR - bottomBorderCounter-1);
                    if(bottomBorderValue>0 && bottomBorderValue!=Float.NaN) break outerloop;
                }
                bottomBorderCounter++;
            }

            outerloop:
            while(leftBorderValue==0){
                for(int y=0; y<h_SR; y++){
                    leftBorderValue = fpSR.getPixelValue(leftBorderCounter, y);
                    if(leftBorderValue>0 && leftBorderValue!=Float.NaN) break outerloop;
                }
                leftBorderCounter++;
            }

            //Check if largest border in stacks
            maxTopBorder = max(topBorderCounter, maxTopBorder);
            maxRightBorder = max(rightBorderCounter, maxRightBorder);
            maxBottomBorder = max(bottomBorderCounter, maxBottomBorder);
            maxLeftBorder = max(leftBorderCounter, maxLeftBorder);

        }

        log.msg("Border Control:");

        if(maxTopBorder==0 && maxRightBorder==0 && maxBottomBorder==0 && maxLeftBorder==0){
            log.msg("\t No borders to crop!");
            borderCrop = false;
            return;
        }

        noCrop = false;
        borderCrop = true;

        // Convert borders to reference image scale
        int roundedTop = (int) (ceil(maxTopBorder/magnification));
        int roundedLeft = (int) (ceil(maxLeftBorder/magnification));
        int roundedBottom = (int) (ceil(maxBottomBorder/magnification));
        int roundedRight = (int) (ceil(maxRightBorder/magnification));


        log.msg("\t Rounded border widths (reference scale): Top = "+roundedTop+", Right = "+roundedRight+", Bottom = "+roundedBottom+", Left = "+roundedLeft);

        int newWidth = w_Ref - roundedRight - roundedLeft;
        int newHeight = h_Ref - roundedBottom - roundedTop;

        imsRef = imsRef.duplicate().crop(roundedLeft, roundedTop, 0, newWidth, newHeight, nSlicesRef);
        imsSR = imsSR.crop(roundedLeft*magnification, roundedTop*magnification, 0, newWidth*magnification, newHeight*magnification, nSlicesSR);

    }

    private void realignImages(){

        int w_Ref = imsRef.getWidth();
        int h_Ref = imsRef.getHeight();

        // Magnify ref image
        ImageStack imsRefResized = magnify(imsRef, w_SR, h_SR);
        FloatProcessor fpRef = imsRefResized.getProcessor(1).convertToFloatProcessor();

        // Cross-correlation realignment
        ImageStack imsSRTranslated = new ImageStack(w_SR, h_SR, nSlicesSR);
        ResultsTable rt = new ResultsTable();

        double maxShiftX = 0, maxShiftY = 0, absMaxShiftX = 0, absMaxShiftY = 0;

        for(int n=1; n<=nSlicesSR; n++) {
            log.status("Checking registration of slice "+n);
            log.progress(n, nSlicesSR);

            FloatProcessor fpSR = imsSR.getProcessor(n).convertToFloatProcessor();

            log.status("calculating cross-correlation...");
            FloatProcessor ipCCM = (FloatProcessor) calculateCrossCorrelationMap(fpRef, fpSR, true);

            if (maxExpectedMisalignment != 0) ipCCM = cropCCM(ipCCM, maxExpectedMisalignment);

            log.status("calculating cross-correlation peaks...");

            float[] shift = getShiftFromCrossCorrelationPeak(ipCCM, MAX_FITTING);
            float shiftX = shift[1];
            float shiftY = shift[2];
            rt.incrementCounter();
            rt.addValue("Frame", n);
            rt.addValue("x-shift (pixels)", shiftX);
            rt.addValue("y-shift (pixels)", shiftY);

            if(abs(shiftX)>absMaxShiftX) {maxShiftX = shiftX; absMaxShiftX = abs(shiftX);}
            if(abs(shiftY)>absMaxShiftY) {maxShiftY = shiftY; absMaxShiftY = abs(shiftY);}

            fpSR.setInterpolationMethod(fpSR.BICUBIC);
            fpSR.resetRoi();
            fpSR.translate(shiftX, shiftY);
            imsSRTranslated.setProcessor(fpSR, n);

        }

        //rt.show("Cross-correlation results");

        log.msg("Registration:");

        //log.msg("maxShiftX = "+maxShiftX+" maxShiftY = "+maxShiftY);
        if(absMaxShiftX<0.5) maxShiftX = 0;
        if(absMaxShiftY<0.5) maxShiftY = 0;

        if(maxShiftX==0 && maxShiftY==0){
            log.msg("\t Images are already registered");
            registrationCrop = false;
            return;
        }

        noCrop = false;
        registrationCrop = true;

        // Crop shifty borders from images
        int roundedMaxShiftX = (int) (ceil(absMaxShiftX/magnification)*signum(maxShiftX));
        int roundedMaxShiftY = (int) (ceil(absMaxShiftY/magnification)*signum(maxShiftY));

        log.msg("\t Maximum misalignment (on reference pixel scale): x = "+roundedMaxShiftX+", y = "+roundedMaxShiftY+". Cropping registration borders.");

        fpRef.resetRoi();

        imsRef = imsRef.duplicate().crop(max(roundedMaxShiftX, 0), max(roundedMaxShiftY, 0), 0,
                w_Ref - abs(roundedMaxShiftX), h_Ref - abs(roundedMaxShiftY), nSlicesRef);

        imsSR = imsSRTranslated.crop(max(roundedMaxShiftX*magnification, 0), max(roundedMaxShiftY*magnification, 0), 0,
                w_SR - abs(roundedMaxShiftX*magnification), h_SR - abs(roundedMaxShiftY*magnification), nSlicesSR);
    }

    float[] getBlockPixels(FloatProcessor fp, double nYB, double nXB){
        int w = fp.getWidth();
        int h = fp.getHeight();

        int blockWidth = w/blocksPerXAxis;
        int blockHeight = h/blocksPerYAxis;

        int xStart = (int) (nXB*blockWidth);
        int yStart = (int) (nYB*blockHeight);

        blockWidth = Math.min(blockWidth, w-xStart);
        blockHeight = Math.min(blockHeight, h-yStart);

        Rectangle rectangle = new Rectangle(xStart, yStart, blockWidth, blockHeight);
        fp.setRoi(rectangle);

        FloatProcessor fpCrop = fp.crop().convertToFloatProcessor();
        fp.resetRoi();

        return (float[]) fpCrop.getPixels();

    }

    FloatProcessor getBlockFp(FloatProcessor fp, double nYB, double nXB){
        int w = fp.getWidth();
        int h = fp.getHeight();

        int blockWidth = w/blocksPerXAxis;
        int blockHeight = h/blocksPerYAxis;

        int xStart = (int) (nXB*blockWidth);
        int yStart = (int) (nYB*blockHeight);

        blockWidth = Math.min(blockWidth, w-xStart);
        blockHeight = Math.min(blockHeight, h-yStart);

        Rectangle rectangle = new Rectangle(xStart, yStart, blockWidth, blockHeight);
        fp.setRoi(rectangle);

        FloatProcessor fpCrop = fp.crop().convertToFloatProcessor();
        fp.resetRoi();

        return fpCrop;
    }

    private class sigmaOptimiseFunction implements UnivariateFunction {

        FloatProcessor fpSRBlock;
        float[] pixelsRefBlock, ones;
        int blockWidthSR, blockHeightSR, blockWidthRef, blockHeightRef;
        ArrayList<Float> sigmaList = new ArrayList<Float>();
        ArrayList<Float> errorList = new ArrayList<Float>();

        public sigmaOptimiseFunction(FloatProcessor fpSRBlock, float[] pixelsRefBlock, float[] ones){
            this.fpSRBlock = fpSRBlock;
            this.pixelsRefBlock = pixelsRefBlock;
            this.ones = ones;
            this.blockWidthSR = fpSRBlock.getWidth();
            this.blockHeightSR = fpSRBlock.getHeight();
            this.blockWidthRef = blockWidthSR/magnification;
            this.blockHeightRef = blockHeightSR/magnification;
        }

        public double value(double sigma) {
            FloatProcessor blurredFp = (FloatProcessor) fpSRBlock.duplicate();
            FloatProcessor blurredOnes = new FloatProcessor(blockWidthSR, blockHeightSR, ones);
            blurredFp.blurGaussian(sigma);
            blurredOnes.blurGaussian(sigma);

            blurredFp = (FloatProcessor) blurredFp.resize(blockWidthRef, blockHeightRef);
            blurredOnes = (FloatProcessor) blurredOnes.resize(blockWidthRef, blockHeightRef);

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRefBlock, (float[]) blurredOnes.getPixels());

            FloatProcessor finalFpSR = (FloatProcessor) fpSRBlock.duplicate();
            finalFpSR.multiply(aB[0]);
            finalFpSR.add(aB[1]);
            finalFpSR.blurGaussian(sigma);
            FloatProcessor finalFpSRResized = (FloatProcessor) finalFpSR.resize(blockWidthRef, blockHeightRef);

            double error = calculateRMSE(pixelsRefBlock, (float[]) finalFpSRResized.getPixels());
            //log.status("Optimising... sigma="+df.format(sigma)+", alpha="+df.format(aB[0])+", beta="+df.format(aB[1])+". Error="+df.format(error));
            sigmaList.add((float) sigma);
            errorList.add((float) error);
            return error;
        }

        public ArrayList<Float> getSigmaList() {
            return sigmaList;
        }

        public ArrayList<Float> getErrorList() {
            return errorList;
        }
    }

    private float[] calculateAlphaBeta(float[] xA, float[] y, float[] oneA){
         /*
            xA = scaled and translated super-resolution
            oneA = scaled ones
            y =  reference
            */

        float N = 0;
        for(int i=0; i<oneA.length; i++){
            N += oneA[i]*oneA[i];
        }

        float nPixels = xA.length;
        assert(nPixels==y.length);
        assert(nPixels==oneA.length);


        float xATxA = 0, xAT1A = 0, yTxA = 0, yT1A = 0;

        for(int i=0; i<nPixels; i++){
            yTxA += y[i]*xA[i];
            yT1A += y[i]*oneA[i];
            xAT1A += xA[i]*oneA[i];
            xATxA += xA[i]*xA[i];
        }

        float numerator = N*yTxA - yT1A*xAT1A;
        float denominator = N*xATxA  - xAT1A*xAT1A;
        float alphaHat = numerator/denominator;
        float betaHat = yT1A/N - alphaHat*(xAT1A/N);

        return new float[] {alphaHat, betaHat};
    }

    private float calculateRMSE(float[] array1, float[] array2){

        int N = array1.length;
        double MSE = 0;

        for(int i=0; i<N; i++){
            MSE += (array1[i]-array2[i])*(array1[i]-array2[i]);
        }
        MSE /= N;

        return (float) Math.sqrt(MSE);
    }

    private FloatProcessor setBlock(FloatProcessor fp, FloatProcessor fpBlock, int xStart, int yStart, int blockWidth, int blockHeight){

        for(int j=0; j<blockHeight; j++){
            for(int i=0; i<blockWidth; i++){
                int y = yStart+j;
                int x = xStart+i;
                fp.setf(x, y, fpBlock.getf(i,j));
            }
        }

        return fp;
    }

    float getIntegratedGaussian(float dx, float dy, float sigma2) {
        float Ex = 0.5f * (erf((dx + 0.5f) / sigma2) - erf((dx - 0.5f) / sigma2));
        float Ey = 0.5f * (erf((dy + 0.5f) / sigma2) - erf((dy - 0.5f) / sigma2));
        float vKernel = Ex * Ey;
        return vKernel;
    }

    float erf(float g) {
        float x = abs(g);
        if (x >= 4.0f)
            return (g > 0.0f) ? 1.0f : -1.0f;

        // constants
        float a1 =  0.254829592f;
        float a2 = -0.284496736f;
        float a3 =  1.421413741f;
        float a4 = -1.453152027f;
        float a5 =  1.061405429f;
        float p  =  0.3275911f;

        // A&S formula 7.1.26
        float t = 1.0f / (1.0f + p*x);
        float y = (float) (1.0f - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x));

        return (g > 0.0f) ? y : -y;
    }

    private FloatProcessor getRSF(float sigma) {
        // calculate the final RSF
        float sigma2 = (float) (ROOT2*abs(sigma));
        int radius = max(((int) sigma) * 3, 1);
        int size = radius * 2 + 1;
        float vKernelSum = 0;

        FloatProcessor fpRSF = new FloatProcessor(size, size);
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {
                float vKernel = getIntegratedGaussian(dx, dy, sigma2);
                vKernelSum+= vKernel;
                fpRSF.setf(dx+radius, dy+radius, vKernel);
            }
        }
        fpRSF.multiply(1./vKernelSum);
        return fpRSF;
    }

    private float[] arrayListToFloatArray(ArrayList<Float> arrayList){
        int nElements = arrayList.size();
        float[] array = new float[nElements];

        for(int n=0; n<nElements; n++) array[n] = arrayList.get(n);
        return array;
    }
}
