package nanoj.squirrel.java.gui.testing;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.gui._BaseDialog_;
import nanoj.core.java.image.filtering.Convolve;

import java.io.IOException;

import static java.lang.Math.max;
import static java.lang.Math.pow;

/**
 * Created by sculley on 12/11/2018.
 */
public class SRGammaTest_ extends _BaseDialog_ {
    private double pixelSizeGT, pixelSizeRef, pixelSizeSR;
    private double refFWHM;
    private double SRFWHM;
    private double minGamma, maxGamma, incGamma;

    Convolve cv = new Convolve();

    @Override
    public boolean beforeSetupDialog(String s) {
        autoOpenImp = true;
        useSettingsObserver = true;

        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Generate reference and different gamma images");
        gd.addNumericField("GT pixel size (nm)", getPrefs("pixelSizeGT", 10), 0);

        gd.addNumericField("Reference image pixel size (nm)", getPrefs("pixelSizeRef", 100), 0);
        gd.addNumericField("Reference PSF FWHM (nm)", getPrefs("refFWHM", 300), 0);

        gd.addNumericField("SR image pixel size (nm)", getPrefs("pixelSizeSR", 20), 0);
        gd.addNumericField("SR FWHM (nm)", getPrefs("SRFWHM", 50), 0);

        gd.addNumericField("Minimum gamma", getPrefs("minGamma", 0.5), 2);
        gd.addNumericField("Maximum gamma", getPrefs("maxGamma", 3), 2);
        gd.addNumericField("Gamma increment", getPrefs("incGamma", 0.5), 2);

    }

    @Override
    public boolean loadSettings() {
        pixelSizeGT = gd.getNextNumber();

        pixelSizeRef = gd.getNextNumber();
        refFWHM = gd.getNextNumber();

        pixelSizeSR = gd.getNextNumber();
        SRFWHM = gd.getNextNumber();

        minGamma = gd.getNextNumber();
        maxGamma = gd.getNextNumber();
        incGamma = gd.getNextNumber();

        setPrefs("pixelSizeGT", pixelSizeGT);
        setPrefs("pixelSizeRef", pixelSizeRef);
        setPrefs("refFWHM", refFWHM);
        setPrefs("pixelSizeSR", pixelSizeSR);
        setPrefs("SRFWHM", SRFWHM);
        setPrefs("minGamma", minGamma);
        setPrefs("maxGamma", maxGamma);
        setPrefs("incGamma", incGamma);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        FloatProcessor fpGT = (FloatProcessor) imp.getProcessor();
        int wGT = fpGT.getWidth();
        int hGT = fpGT.getHeight();
        float[] pixelsGT = (float[]) fpGT.getPixels();
        int nPixelsGT = pixelsGT.length;

        double refSigma = refFWHM/2.35482;
        double SRSigma = SRFWHM/2.35482;

        double refSigmaGTPx = refSigma/pixelSizeGT;

        // create reference image

        ///normalize so that max intensity is 1000 to prevent loss during convolution

        float[] pixelsRef = normalise(pixelsGT.clone(), 1000);

        /// blur
        FloatProcessor fpRef = new FloatProcessor(wGT, hGT, pixelsRef);
        fpRef.blurGaussian(refSigmaGTPx);

        /// resize
        int magnificationRef = (int) (pixelSizeRef/pixelSizeGT);
        int wRef = wGT/magnificationRef;
        int hRef = hGT/magnificationRef;
        fpRef = (FloatProcessor) fpRef.resize(wRef, hRef);

        /// display
        new ImagePlus("Reference image", fpRef).show();

        // create SR images
        int magnificationSR = (int) (pixelSizeSR/pixelSizeGT);
        int wSR = wGT/magnificationSR;
        int hSR = hGT/magnificationSR;

        double SRSigmaGTPx = SRSigma/pixelSizeGT;
        double SRFWHMGTPx = SRFWHM/pixelSizeGT;

        int nImages = 0;
        for(double g=minGamma; g<=maxGamma; g+=incGamma) nImages++;

        ImageStack imsSR = new ImageStack(wSR, hSR, nImages);

        int imCounter = 1;

        FloatProcessor fpSR = fpGT.duplicate().convertToFloatProcessor();
        fpSR.blurGaussian(SRSigmaGTPx);
        fpSR.setInterpolationMethod(ImageProcessor.NONE);
        fpSR = fpSR.resize(wSR, hSR).convertToFloatProcessor();

        for(double g=minGamma; g<=maxGamma; g+=incGamma){

            FloatProcessor fp = fpSR.duplicate().convertToFloatProcessor();

            float[] pixels = (float[]) fp.getPixels();
            for(int i=0; i<pixels.length; i++) pixels[i] = (float) pow(pixels[i], g);
            fp.setPixels(normalise(pixels, 100));

            imsSR.setProcessor(fp, imCounter);
            imsSR.setSliceLabel("g="+g, imCounter);
            imCounter++;
        }

        new ImagePlus("SR images", imsSR).show();

    }

    float[] normalise(float[] pixels, float targetMax){
        float maxVal = 0;
        for(int i=0; i<pixels.length; i++) maxVal = max(maxVal, pixels[i]);
        for(int i=0; i<pixels.length; i++) pixels[i] *= (targetMax/maxVal);
        return pixels;
    }
}
