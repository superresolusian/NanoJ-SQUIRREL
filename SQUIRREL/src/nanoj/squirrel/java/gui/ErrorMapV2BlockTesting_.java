package nanoj.squirrel.java.gui;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.gui._BaseDialog_;
import org.apache.commons.math3.analysis.UnivariateFunction;

import java.io.IOException;
import java.util.ArrayList;

import static java.lang.Math.max;
import static java.lang.Math.min;

/**
 * Created by sculley on 20/09/2018.
 */
public class ErrorMapV2BlockTesting_ extends _BaseDialog_ {

    String titleRefImage = "", titleRSFImage = "", titleSRImage = "";
    String[] imageTitles;

    int blocksPerXAxis, blocksPerYAxis;
    double maxSigma;

    ImagePlus impRef, impSR;
    ImageStack imsRef, imsSR;
    int w_SR, h_SR, w_Ref, h_Ref;
    int magnification, magnification2;

    int blockWidthRef, blockHeightRef, blockWidthSR, blockHeightSR;
    int bufferSize;

    @Override
    public boolean beforeSetupDialog(String arg) {

        autoOpenImp = false;
        useSettingsObserver = true;
        int nImages = WindowManager.getImageCount();
        if(nImages<2){
            log.error("At least 2 images are required to run this code!");
            return false;
        }

        imageTitles = WindowManager.getImageTitles();

        return true;
    }

    @Override
    public void setupDialog() {

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

        gd.addNumericField("Blocks per axis (default: 10)", getPrefs("blocksPerAxis", 10), 0);
        gd.addNumericField("Max sigma (nm)", getPrefs("maxSigma", 200),0);

    }

    @Override
    public boolean loadSettings() {
        titleRefImage = gd.getNextChoice();
        titleSRImage = gd.getNextChoice();

        int blocksPerAxis = (int) max(gd.getNextNumber(),1);
        maxSigma = gd.getNextNumber();

        setPrefs("titleRefImage", titleRefImage);
        setPrefs("titleSRImage", titleSRImage);
        setPrefs("blocksPerAxis", blocksPerAxis);
        setPrefs("maxSigma", maxSigma);

        prefs.savePreferences();

        impRef = WindowManager.getImage(titleRefImage);
        impSR = WindowManager.getImage(titleSRImage);

        imsRef = impRef.getImageStack().convertToFloat();
        imsSR = impSR.getImageStack().convertToFloat();


        w_Ref = impRef.getWidth();
        h_Ref = impRef.getHeight();

        w_SR = imsSR.getWidth();
        h_SR = imsSR.getHeight();

        magnification = w_SR / imsRef.getWidth();
        magnification2 = magnification * magnification;

        double _magnification = ((double) w_SR) / imsRef.getWidth();
        if (magnification != _magnification) {
            log.error("SR image size needs to be an integer multiple of the reference image size!");
            return false;
        }

        blocksPerXAxis = blocksPerAxis;
        blocksPerYAxis = (int) (blocksPerXAxis * (((double) h_Ref)/w_Ref));

        blockWidthRef =  w_Ref/ blocksPerXAxis;
        blockHeightRef = h_SR / blocksPerYAxis;

        blockWidthSR = blockWidthRef*magnification;
        blockHeightSR = blockHeightRef*magnification;

        bufferSize = magnification;

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        // Set up reference float processors
        FloatProcessor fpRef = imsRef.getProcessor(1).convertToFloatProcessor();
        fpRef.resetRoi();

        w_Ref = fpRef.getWidth();
        h_Ref = fpRef.getHeight();

        int nPixelsSR = w_SR*h_SR;
        int nSlicesSR = 1;

        // Set up pixel arrays for optimization
        float[] pixelsRef = (float[]) fpRef.getPixels();
        float[] ones = new float[nPixelsSR];
        for(int i=0; i<nPixelsSR; i++){ones[i] = 1;}

        // Set up pixel array for error map generation
        FloatProcessor fpRefScaledToSR = (FloatProcessor) fpRef.duplicate();
        fpRefScaledToSR.setInterpolationMethod(ImageProcessor.BICUBIC);
        fpRefScaledToSR = (FloatProcessor) fpRefScaledToSR.resize(w_SR, h_SR);
        float[] pixelsRefScaledToSR = (float[]) fpRefScaledToSR.getPixels();

        // Set up output
        int maxRSFWidth = 0;
        FloatProcessor[] fpRSFArray = new FloatProcessor[nSlicesSR];

        ImageStack imsSRConvolved = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsEMap = new ImageStack(w_SR, h_SR, nSlicesSR);
        ImageStack imsSRNormalised = new ImageStack(w_SR, h_SR, nSlicesSR);
        ResultsTable rt = new ResultsTable();

        /*
        For maximum expected emission wavelength λem = 700nm, PSF FWHM = 305nm (Born & Wolf)
        => sigma of Gaussian PSF = 305/2.35482 nm = 130nm
        If no calibration data, assume arbitrary pixel size of 100nm
        => sigma of Gaussian PSF = 130/100 = 1.3 pixels
        Add 10% leeway => maximum sigma = 1.3*1.1 = 1.43 pixels
        Sigma is calcuated on upsampled grid therefore multiply by magnification

        However, this never seems to be quite enough blur. Instead, choose max sigma = 200nm...
         */

        float maxSigmaBoundary = 0;

        // Adjust maxSigmaBoundary if there is calibration data
        String pixelUnitRef = impRef.getCalibration().getUnit();
        double pixelSizeNm;
        if(pixelUnitRef.equals("nm")){
            log.msg("Detected nm!");
            pixelSizeNm = impRef.getCalibration().pixelWidth;
            maxSigmaBoundary = (float) (maxSigma/pixelSizeNm)*magnification;
            log.msg("Setting maximum sigma to "+maxSigma+"nm ("+maxSigmaBoundary+" pixels) to avoid overblurring");
        }
        else if(pixelUnitRef.equals("micron")|| pixelUnitRef.equals("microns") || pixelUnitRef.equals("µm")){
            log.msg("Detected µm!");
            pixelSizeNm = impRef.getCalibration().pixelWidth * 1000;
            maxSigmaBoundary = (float) (maxSigma/pixelSizeNm)*magnification;
            log.msg("Setting maximum sigma to "+maxSigma+"nm ("+maxSigmaBoundary+" pixels) to avoid overblurring");
        }
        else{
            //assign arbitrary pixel size of 100nm
            pixelSizeNm = 100;
            maxSigmaBoundary = (float) (maxSigma/pixelSizeNm)*magnification;
            log.msg("Setting maximum sigma to "+maxSigmaBoundary+" pixels to avoid overblurring");
        }


        for (int nYB = 0; nYB < blocksPerYAxis; nYB++) {
            for (int nXB = 0; nXB < blocksPerXAxis; nXB++) {

            }
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

    private class sigmaOptimiseFunctionBlocks implements UnivariateFunction {

        int xB, yB;
        FloatProcessor fpSR;
        float[] pixelsRef, ones;
        ArrayList<Float> sigmaList = new ArrayList<Float>();
        ArrayList<Float> errorList = new ArrayList<Float>();

        public sigmaOptimiseFunctionBlocks(int xB, int yB, FloatProcessor fpSR, float[] pixelsRef, float[] ones){
            this.xB = xB;
            this.yB = yB;
            this.fpSR = fpSR;
            this.pixelsRef = pixelsRef;
            this.ones = ones;
        }

        private FloatProcessor getROIWithBuffer(FloatProcessor fp, int x, int y, int w, int h, int bufferSize) {

            int x0 = max(0, x-bufferSize);
            int x1 = min(w_SR, x+w+bufferSize);
            int y0 = max(0, y-bufferSize);
            int y1 = min(h_SR, y+h+bufferSize);

            int w_ = x1-x0;
            int h_ = y1-y0;

            FloatProcessor fpCrop = new FloatProcessor(w_,h_);
            for (int j=0; j<h_; j++) {
                for (int i=0; i<w_; i++) {
                    fpCrop.setf(i, j, fp.getf(x0+i, y0+j));
                }
            }
            return fpCrop;
        }

        private FloatProcessor getROIWithoutBuffer(FloatProcessor fp, int x, int y, int w, int h, int bufferSize) {

            int x0 = min(x, bufferSize);
            int y0 = min(y, bufferSize);

            FloatProcessor fpCrop = new FloatProcessor(w,h);
            for (int j=0; j<h; j++) {
                for (int i=0; i<w; i++) {
                    fpCrop.setf(i, j, fp.getf(x0+i, y0+j));
                }
            }
            return fpCrop;
        }


        public double value(double sigma) {
            //grab SR block + 5 pixel buffer region from SR image.
            int xStartSR = xB * blockWidthSR;
            int yStartSR = yB * blockHeightSR;

            //account for blocks which extend past edges of image
            int thisBlockWidthSR = Math.min(blockWidthSR, w_SR-xStartSR);
            int thisBlockHeightSR = Math.min(blockHeightSR, h_SR-yStartSR);

            FloatProcessor blurredFpBlock = getROIWithBuffer(fpSR, xStartSR, yStartSR, thisBlockWidthSR, thisBlockHeightSR, bufferSize);
            float[] ones = new float[blurredFpBlock.getWidth()*blurredFpBlock.getHeight()];
            for(int i=0; i<ones.length; i++) ones[i] = 1;
            FloatProcessor blurredOnes = new FloatProcessor(blurredFpBlock.getWidth(), blurredFpBlock.getHeight(), ones);

            blurredFpBlock.blurGaussian(sigma);
            blurredOnes.blurGaussian(sigma);

            blurredFpBlock = getROIWithoutBuffer(blurredFpBlock, xStartSR, yStartSR, thisBlockWidthSR, thisBlockHeightSR, bufferSize);
            FloatProcessor blurredFp = (FloatProcessor) blurredFpBlock.resize(blockWidthRef, blockHeightRef);
            blurredOnes = (FloatProcessor) blurredOnes.resize(w_Ref, h_Ref);

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRef, (float[]) blurredOnes.getPixels());


            FloatProcessor finalFpSR = (FloatProcessor) fpSR.duplicate();
            finalFpSR.multiply(aB[0]);
            finalFpSR.add(aB[1]);
            finalFpSR.blurGaussian(sigma);
            FloatProcessor finalFpSRResized = (FloatProcessor) finalFpSR.resize(w_Ref, h_Ref);

            double error = calculateRMSE(pixelsRef, (float[]) finalFpSRResized.getPixels());
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

    //helper function lifted from imageshelper
    protected static boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
    }
}
