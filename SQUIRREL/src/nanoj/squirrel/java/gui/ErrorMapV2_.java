package nanoj.squirrel.java.gui;

import histogram2.HistogramMatcher;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.kernels.Kernel_SquirrelSwarmOptimizer;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;
import nanoj.squirrel.java.gui.tools.SetMaximumStackSize_;
import nanoj.squirrel.java.minimizers.GaussianFitMinimizer;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.*;
import static nanoj.core.java.array.ArrayCasting.toArray;
import static nanoj.core.java.array.ArrayMath.calculateMSE;
import static nanoj.core.java.array.ArrayMath.calculatePPMCC;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.MAX_FITTING;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.getShiftFromCrossCorrelationPeak;
import static nanoj.core.java.image.transform.CrossCorrelationMap.calculateCrossCorrelationMap;
import static nanoj.core.java.image.transform.CrossCorrelationMap.cropCCM;
import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_Errors;
import static nanoj.squirrel.java.gui.ImagesHelper.magnify;

/**
 * Created by sculley on 02/06/2017.
 */
public class ErrorMapV2_ extends _BaseSQUIRRELDialog_ {

    private static final float ROOT2 = (float) Math.sqrt(2);
    private static Kernel_SquirrelSwarmOptimizer kPSOErrorMap = new Kernel_SquirrelSwarmOptimizer();

    String titleRefImage = "", titleRSFImage = "", titleSRImage = "", noRSFString = "-- RSF unknown, estimate via optimisation --";
    boolean showAdvancedSettings, _showAdvancedSettings = false;

    protected ErrorMap_ExtraSettings_ errorMap_ExtraSettings = new ErrorMap_ExtraSettings_();
    protected ErrorMap_ExtraSettings_ _errorMap_ExtraSettings;

    protected SetMaximumStackSize_ maximumStackSize = new SetMaximumStackSize_();
    ArrayList<String> titles = new ArrayList<String>();
    String[] imageTitles;

    boolean framePurge, borderControl;
    boolean doRegistration;
    int maxExpectedMisalignment;
    int maxSRStackSize;
    boolean showIntensityNormalised, showConvolved, showRSF, showPositiveNegative;

    HistogramMatcher histogramMatcher = new HistogramMatcher();

    double sigmaGuess;
    ImagePlus impRef, impSR, impRSF;
    ImageStack imsRef, imsSR, imsRSF;
    int nSlicesRef, nSlicesSR, nSlicesRSF;
    int w_SR, h_SR, w_Ref, h_Ref;
    int magnification, magnification2;
    boolean noCrop = true;
    private boolean borderCrop = false;
    private boolean registrationCrop = false;
    private int maxMag;

    @Override
    public boolean beforeSetupDialog(String arg) {
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

        setPrefs("titleRefImage", titleRefImage);
        setPrefs("titleSRImage", titleSRImage);
        setPrefs("titleRSFImage", titleRSFImage);
        setPrefs("maxMag", maxMag);
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
    public void execute() {

        long startTime = System.currentTimeMillis();


        // Check and purge homogeneous/empty slices
        if(framePurge) {
            ArrayList<Integer> slicesToDelete = new ArrayList<Integer>();

            for (int n = 1; n <= nSlicesSR; n++) {

                log.status("Checking stack, frame " + n);
                log.progress(n, nSlicesSR);

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

        // Set up pixel arrays for optimization
        float[] pixelsRef = (float[]) fpRef.getPixels();
        float[] ones = new float[nPixelsSR];
        for(int i=0; i<nPixelsSR; i++){ones[i] = 1;}
        ImageProcessor ipRef = imsRef.getProcessor(1);
        int[] histRef = ipRef.getHistogram();

        // Set up pixel array for error map generation
        FloatProcessor fpRefScaledToSR = (FloatProcessor) fpRef.duplicate();
        fpRefScaledToSR.setInterpolationMethod(ImageProcessor.BICUBIC);
        fpRefScaledToSR = (FloatProcessor) fpRefScaledToSR.resize(w_SR, h_SR);
        float[] pixelsRefScaledToSR = (float[]) fpRefScaledToSR.getPixels();

        // Set up output
        int maxRSFWidth = 0;
        FloatProcessor[] fpRSFArray = new FloatProcessor[nSlicesSR];

        ImageStack imsSRConvolved = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsSRConvolved_histMatched = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsEMap = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsEMap_histMatched = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsSRNormalised = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsSRNormalised_histMatched = new ImageStack(w_SR, h_SR, nSlicesSR);
        ResultsTable rt = new ResultsTable();

        /*
        Nyquist calculations:
        Resolution of microscope = FWHM of PSF
        => If correctly sampled, resolution = 2 pixels. More likely to oversample than undersample so make generous
        assumption that resolution = 4 pixels
        => Sigma of Gaussian PSF would be 4/2.35482 pixels = 1.6986
        => Calculating sigma on upsampled grid therefore multiply by magnification
         */

        //float nyquistFactor = (4/2.35482f)*magnification;
        float nyquistFactor = 20;

        long loopStart = System.nanoTime();

        for(int n=1; n<=nSlicesSR; n++){

            log.msg("-----------------------------------");
            log.msg("Processing super-resolution frame "+n);
            log.msg("-----------------------------------");

            // Get SR FloatProcessor
            FloatProcessor fpSR = imsSR.getProcessor(n).convertToFloatProcessor();
            ImageProcessor ipSR = imsSR.getProcessor(n);

            log.msg("*****LINEAR RESCALING*****");

            // UNIVARIATE OPTIMIZER - LINEAR MATCHING
            /// setup optimizer
            sigmaOptimiseFunction f =  new ErrorMapV2_.sigmaOptimiseFunction(fpSR, pixelsRef, ones);
            UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
            /// run optimizer
            UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                    new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0,nyquistFactor)); //NYQUIST ASSUMED
            float sigma_linear = (float) result.getPoint();
            log.msg("Best sigma is: "+sigma_linear);
            log.msg("Best error is: "+result.getValue());

            float[] errorList = toArray(f.getErrorList(), 1.0f);
            float[] sigmaList = toArray(f.getSigmaList(), 1.0f);

            Plot plot = new Plot("Brent optimiser - frame "+n+": Error vs Sigma", "Sigma", "Error");
            plot.addPoints(sigmaList, errorList, Plot.CIRCLE);
            plot.show();

            // GET ALPHA AND BETA
            FloatProcessor blurredFp = (FloatProcessor) fpSR.duplicate();
            FloatProcessor blurredOnes = new FloatProcessor(w_SR, h_SR, ones);
            blurredFp.blurGaussian(sigma_linear);
            blurredOnes.blurGaussian(sigma_linear);

            blurredFp = (FloatProcessor) blurredFp.resize(w_Ref, h_Ref);
            blurredOnes = (FloatProcessor) blurredOnes.resize(w_Ref, h_Ref);

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRef, (float[]) blurredOnes.getPixels());
            float alpha = aB[0];
            float beta = aB[1];

            log.msg("Alpha is: "+alpha+", beta is: "+beta);

            log.msg("*****HISTOGRAM MATCHING*****");
            /////////////////////
            //Histogram Matching
            /// setup optimizer
            sigmaOptimiseMatchHistogram f_mh =  new ErrorMapV2_.sigmaOptimiseMatchHistogram(fpSR, pixelsRef, histRef);
            UnivariateOptimizer optimizer_mh = new BrentOptimizer(1e-10, 1e-14);
            /// run optimizer
            UnivariatePointValuePair result_mh = optimizer.optimize(new MaxEval(1000),
                    new UnivariateObjectiveFunction(f_mh), GoalType.MINIMIZE, new SearchInterval(0,nyquistFactor)); //NYQUIST ASSUMED
            float sigma_matched = (float) result_mh.getPoint();
            log.msg("Best sigma is: "+sigma_matched);
            log.msg("Best error is: "+result_mh.getValue());

            float[] errorList_mh = toArray(f_mh.getErrorList(), 1.0f);
            float[] sigmaList_mh = toArray(f_mh.getSigmaList(), 1.0f);

            Plot plot_mh = new Plot("Brent optimiser - histogram matching - frame "+n+": Error vs Sigma", "Sigma", "Error");
            plot_mh.addPoints(sigmaList_mh, errorList_mh, Plot.CIRCLE);
            plot_mh.show();

            // GET ALPHA AND BETA
            FloatProcessor blurredFp_mh = (FloatProcessor) fpSR.duplicate();
            FloatProcessor blurredOnes_mh = new FloatProcessor(w_SR, h_SR, ones);
            blurredFp_mh.blurGaussian(sigma_matched);
            blurredOnes_mh.blurGaussian(sigma_matched);

            blurredFp_mh = (FloatProcessor) blurredFp_mh.resize(w_Ref, h_Ref);
            blurredOnes_mh = (FloatProcessor) blurredOnes_mh.resize(w_Ref, h_Ref);

            float[] aB_mh = calculateAlphaBeta((float[]) blurredFp_mh.getPixels(), pixelsRef, (float[]) blurredOnes_mh.getPixels());
            float alpha_mh = aB_mh[0];
            float beta_mh = aB_mh[1];

            log.msg("Alpha is: "+alpha_mh+", beta is: "+beta_mh);


            // POPULATE OUTPUT STACKS

            /// intensity-scaled stack
            FloatProcessor fpSRIntensityScaled = (FloatProcessor) fpSR.duplicate();
            fpSRIntensityScaled.multiply(alpha);
            fpSRIntensityScaled.add(beta);
            imsSRNormalised.setProcessor(fpSRIntensityScaled, n);

            FloatProcessor fpSRIntensityScaled_mh = (FloatProcessor) fpSR.duplicate();
            fpSRIntensityScaled_mh.multiply(alpha);
            fpSRIntensityScaled_mh.add(beta);
            imsSRNormalised_histMatched.setProcessor(fpSRIntensityScaled_mh, n);

            /// intensity-scaled and convolved stack
            FloatProcessor fpSRIntensityScaledBlurred = (FloatProcessor) fpSRIntensityScaled.duplicate();
            fpSRIntensityScaledBlurred.blurGaussian(sigma_linear);
            imsSRConvolved.setProcessor(fpSRIntensityScaledBlurred, n);

            FloatProcessor fpSRIntensityScaledBlurred_mh = (FloatProcessor) fpSRIntensityScaled_mh.duplicate();
            fpSRIntensityScaledBlurred_mh.blurGaussian(sigma_matched);
            imsSRConvolved_histMatched.setProcessor(fpSRIntensityScaledBlurred_mh, n);

            /// Generate RSF for this image
            log.status("Generating RSF");
            FloatProcessor fpRSF = getRSF(sigma_linear/magnification);
            fpRSFArray[n-1] = fpRSF;
            maxRSFWidth = max(maxRSFWidth, fpRSF.getWidth());

            // CALCULATE METRICS AND MAP
            log.status("Calculating similarity...");

            /// metrics
            FloatProcessor fpIntensityScaledBlurred_RefSize = (FloatProcessor) fpSRIntensityScaledBlurred.resize(w_Ref, h_Ref);
            float[] pixelsIntensityScaledBlurred_RefSize = (float[]) fpIntensityScaledBlurred_RefSize.getPixels();
            double globalRMSE = sqrt(calculateMSE(pixelsIntensityScaledBlurred_RefSize, pixelsRef));
            double globalPPMCC = calculatePPMCC(pixelsIntensityScaledBlurred_RefSize, pixelsRef, true);

            /// error map
            float[] pixelsEMap = new float[nPixelsSR];
            float[] pixelsSRC = (float[]) fpSRIntensityScaledBlurred.getPixels();

            float maxRef = -Float.MAX_VALUE;
            for(int p=0; p<nPixelsSR; p++){
                float vRef = pixelsRefScaledToSR[p];
                float vSRC = pixelsSRC[p];

                maxRef = max(maxRef, vRef); //why
                if(showPositiveNegative) pixelsEMap[p] = (vRef - vSRC);
                else pixelsEMap[p] = abs(vRef-vSRC);
            }
            imsEMap.setProcessor(new FloatProcessor(w_SR, h_SR, pixelsEMap), n);

            /// put error values into table
            rt.incrementCounter();
            rt.addValue("Frame", n);
            rt.addValue("RSP (Resolution Scaled Pearson-Correlation)", globalPPMCC);
            rt.addValue("RSE (Resolution Scaled Error)", globalRMSE);

            if(nSlicesSR>1) {
                double frameTime = ((System.nanoTime() - loopStart) / n) / 1e9;
                double remainingTime = frameTime * (nSlicesSR - n);
                int _h = (int) (remainingTime / 3600);
                int _m = (int) (((remainingTime % 86400) % 3600) / 60);
                int _s = (int) (((remainingTime % 86400) % 3600) % 60);
                log.msg("Estimated time remaining to complete analysis:");
                log.msg("\t" + String.format("%02d:%02d:%02d", _h, _m, _s));
            }
        }

        /// BUILD RSF STACK
        ImageStack imsRSF_rescaled = new ImageStack(maxRSFWidth, maxRSFWidth, nSlicesSR);
        int centreLargest = maxRSFWidth/2;

        for(int n=1; n<=nSlicesSR; n++){
            FloatProcessor fp = fpRSFArray[n-1];
            int centre = fp.getWidth()/2;
            FloatProcessor fpInsert = new FloatProcessor(maxRSFWidth, maxRSFWidth);
            fpInsert.insert(fp, centreLargest-centre, centreLargest-centre);
            imsRSF_rescaled.setProcessor(fpInsert, n);
        }

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
            new ImagePlus(titleSRImage +" - " +borderString +registrationString+" - intensity-normalised", imsSRNormalised).show();
            new ImagePlus(titleSRImage +" - " +borderString +registrationString+" - histomatch - intensity-normalised", imsSRNormalised_histMatched).show();
        }

        if(showConvolved){
            new ImagePlus(titleSRImage+" - Convolved with RSF", imsSRConvolved).show();
            new ImagePlus(titleSRImage+" - histomatch - Convolved with RSF", imsSRConvolved_histMatched).show();

        }
        if(showRSF) new ImagePlus("Optimised RSF stack", imsRSF_rescaled).show();

        rt.show("RSP and RSE values");

        ImagePlus impEMap = new ImagePlus(titleSRImage+" - Resolution Scaled Error-Map", imsEMap);
        applyLUT_SQUIRREL_Errors(impEMap);
        impEMap.show();

        log.msg("SQUIRREL analysis took "+(System.currentTimeMillis()-startTime)/1000+"s");

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

    private double[] floatArrayToDoubleArray(Object pixels) {
        float[] floatArray = (float[]) pixels;
        double[] doubleArray = new double[floatArray.length];
        for(int i=0; i<floatArray.length; i++){
            doubleArray[i] = (double) floatArray[i];
        }
        return doubleArray;
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

    //helper function lifted from imageshelper
    protected static boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
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

    private class sigmaOptimiseFunction implements UnivariateFunction{

        FloatProcessor fpSR;
        float[] pixelsRef, ones;
        ArrayList<Float> sigmaList = new ArrayList<Float>();
        ArrayList<Float> errorList = new ArrayList<Float>();

        public sigmaOptimiseFunction(FloatProcessor fpSR, float[] pixelsRef, float[] ones){
            this.fpSR = fpSR;
            this.pixelsRef = pixelsRef;
            this.ones = ones;
        }

        public double value(double sigma) {
            FloatProcessor blurredFp = (FloatProcessor) fpSR.duplicate();
            FloatProcessor blurredOnes = new FloatProcessor(w_SR, h_SR, ones);
            blurredFp.blurGaussian(sigma);
            blurredOnes.blurGaussian(sigma);

            blurredFp = (FloatProcessor) blurredFp.resize(w_Ref, h_Ref);
            blurredOnes = (FloatProcessor) blurredOnes.resize(w_Ref, h_Ref);

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRef, (float[]) blurredOnes.getPixels());

            FloatProcessor finalFpSR = (FloatProcessor) fpSR.duplicate();
            finalFpSR.multiply(aB[0]);
            finalFpSR.add(aB[1]);
            finalFpSR.blurGaussian(sigma);
            FloatProcessor finalFpSRResized = (FloatProcessor) finalFpSR.resize(w_Ref, h_Ref);

            double error = calculateRMSE(pixelsRef, (float[]) finalFpSRResized.getPixels());
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

    private class sigmaOptimiseMatchHistogram implements UnivariateFunction{

        ImageProcessor ipSR;
        float[] pixelsRef;
        int[] histRef;

        ArrayList<Float> sigmaList = new ArrayList<Float>();
        ArrayList<Float> errorList = new ArrayList<Float>();

        public sigmaOptimiseMatchHistogram(ImageProcessor ipSR, float[] pixelsRef, int[] histRef){
            this.ipSR = ipSR;
            this.pixelsRef = pixelsRef;
            this.histRef = histRef;
        }

        public double value(double sigma) {
            FloatProcessor blurredFp = (FloatProcessor) ipSR.duplicate();
            blurredFp.blurGaussian(sigma);
            blurredFp = (FloatProcessor) blurredFp.resize(w_Ref, h_Ref);
            int[] blurredSRHist = blurredFp.getHistogram();

            int[] matchedHistograms = histogramMatcher.matchHistograms(blurredSRHist, histRef);
            blurredFp.applyTable(matchedHistograms);

            double error = calculateRMSE(pixelsRef, (float[]) blurredFp.getPixels());
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
}
