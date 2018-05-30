package nanoj.squirrel.java.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.measure.CurveFitter;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.image.filtering.Convolve;
import nanoj.kernels.Kernel_SquirrelSwarmOptimizer;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;
import nanoj.squirrel.java.gui.tools.SetMaximumStackSize_;
import nanoj.squirrel.java.minimizers.GaussianFitMinimizer;
import org.w3c.dom.Window;

import java.awt.*;
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
import static nanoj.kernels.Kernel_SquirrelSwarmOptimizer.*;
import static nanoj.squirrel.java.gui.ImagesHelper.magnify;

/**
 * Created by sculley on 02/06/2017.
 */
public class ErrorMapSimple_ extends _BaseSQUIRRELDialog_ {

    private static Kernel_SquirrelSwarmOptimizer kPSOErrorMap = new Kernel_SquirrelSwarmOptimizer();

    String titleRefImage = "", titleRSFImage = "", titleSRImage = "", noRSFString = "-- RSF unknown, estimate via optimisation --";
    boolean showAdvancedSettings, _showAdvancedSettings = false;

    protected ErrorMap_ExtraSettings_ errorMap_ExtraSettings = new ErrorMap_ExtraSettings_();
    protected ErrorMap_ExtraSettings_ _errorMap_ExtraSettings;

    protected SetMaximumStackSize_ maximumStackSize = new SetMaximumStackSize_();
    ArrayList<String> titles = new ArrayList<String>();
    String[] imageTitles;

    boolean framePurge, borderControl;
    int maxExpectedMisalignment;
    int maxSRStackSize;
    boolean showIntensityNormalised, showConvolved, showRSF, showPositiveNegative;

    double sigmaGuess, alphaGuess, betaGuess;
    ImagePlus impRef, impSR, impRSF;
    ImageStack imsRef, imsSR, imsRSF;
    int nSlicesRef, nSlicesSR, nSlicesRSF;
    int w_SR, h_SR;
    int magnification, magnification2;
    boolean noCrop = true;
    private boolean visualiseParameterEvolution = false;
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
        log.status("Registering and cropping images...");
        realignImages();

        // update width and height
        w_SR = imsSR.getWidth();
        h_SR = imsSR.getHeight();

        // Set up reference float processors
        FloatProcessor fpRef = imsRef.getProcessor(1).convertToFloatProcessor();
        fpRef.resetRoi();

        FloatProcessor fpRef_scaledToSR = (FloatProcessor) fpRef.duplicate();
        fpRef_scaledToSR.setInterpolationMethod(ImageProcessor.BICUBIC);
        fpRef_scaledToSR = (FloatProcessor) fpRef_scaledToSR.resize(w_SR, h_SR);
        float[] pixelsRef_scaledToSR = (float[]) fpRef_scaledToSR.getPixels();

        // Set up output
        int maxRSFWidth = 0;
        FloatProcessor[] fpRSFArray = new FloatProcessor[nSlicesSR];

        ImageStack imsSRConvolved = new ImageStack(w_SR, h_SR, nSlicesSR);

        int nPixelsSR = w_SR*h_SR;
        ImageStack imsEMap = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsSRNormalised = new ImageStack(w_SR, h_SR, nSlicesSR);
        ResultsTable rt = new ResultsTable();

        long loopStart = System.nanoTime();

        ArrayList<float[]> errorEvolution = new ArrayList<float[]>(nSlicesSR);

        float minSigma, maxSigma, sigmaStep;

        if(sigmaGuess==0){
            minSigma = 0.1f;
            maxSigma = 10;
            sigmaStep = 0.1f;
        }
        else{
            minSigma = (float) (sigmaGuess*0.25);
            maxSigma = (float) (sigmaGuess*1.75);
            sigmaStep = (float) (sigmaGuess*0.05);
        }

        int nSigmaSteps = (int) ((maxSigma-minSigma)/sigmaStep);
        int nPixelsRef = fpRef.getWidth()*fpRef.getHeight();
        float[] pixelsRef = (float[]) fpRef.getPixels();

        float[] sigmas = new float[nSigmaSteps];
        float[] errors = new float[nSigmaSteps];

        for(int n=1; n<=nSlicesSR; n++) {

            log.msg("-----------------------------------");
            log.msg("Processing super-resolution frame " + n);
            log.msg("-----------------------------------");

            // Set up SR float processors
            FloatProcessor fpSR = imsSR.getProcessor(n).convertToFloatProcessor();
            fpSR.resetRoi();
            FloatProcessor fpSR_scaledToRef = (FloatProcessor) fpSR.duplicate();
            fpSR_scaledToRef = (FloatProcessor) fpSR_scaledToRef.resize(fpRef.getWidth(), fpRef.getHeight());

            for (int m = 0; m < nSigmaSteps; m++) {

                // get alpha and beta
                CurveFitter cf = new CurveFitter(floatArrayToDoubleArray(fpSR_scaledToRef.getPixels()), floatArrayToDoubleArray(fpRef.getPixels()));
                cf.doFit(cf.STRAIGHT_LINE);
                double beta = cf.getParams()[0];
                double alpha = cf.getParams()[1];

                // rescale image
                fpSR_scaledToRef = linearRescale(fpSR.duplicate().resize(fpRef.getWidth(), fpRef.getHeight()), alpha, beta);

                // convolve (nanoj)
                Convolve cv = new Convolve();
                FloatProcessor
                cv.

                // convolve with gaussian of size sigma (imagej)
                float thisSigma = m * sigmaStep + minSigma;
                fpSR_scaledToRef.blurGaussian(thisSigma);

                // calculate error
                float[] pixelsSR_scaledToRef = (float[]) fpSR_scaledToRef.getPixels();
                float error = 0;

                for (int e = 0; e < nPixelsRef; e++) {
                    error += pow(pixelsRef[e] - pixelsSR_scaledToRef[e], 2) / nPixelsRef;
                }

                sigmas[m] = thisSigma;
                errors[m] = error;
            }

            errorEvolution.add(errors);
        }

        Plot plot = new Plot("sigma vs error", "sigma", "error", sigmas, errors);

        plot.show();

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

    FloatProcessor linearRescale(ImageProcessor ip, double gradient, double offset){
        FloatProcessor fp = ip.convertToFloatProcessor();

        float[] pixels = (float[]) fp.getPixelsCopy();

        for(int i=0; i<pixels.length; i++){
            float v = pixels[i];
            pixels[i] = (float) (v*gradient + offset);
        }

        return new FloatProcessor(fp.getWidth(), fp.getHeight(), pixels);

    }


    //helper function lifted from imageshelper
    protected static boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
    }
}
