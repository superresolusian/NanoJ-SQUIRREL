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

/**
 * Created by sculley on 12/11/2018.
 */
public class SRLinearityTest_ extends _BaseDialog_ {
    private double pixelSizeGT, pixelSizeRef, pixelSizeSR;
    private double refFWHM;
    private double SRFWHM;

    Convolve cv = new Convolve();

    @Override
    public boolean beforeSetupDialog(String s) {
        autoOpenImp = true;
        useSettingsObserver = true;

        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Generate reference and different linearity images");
        gd.addNumericField("GT pixel size (nm)", getPrefs("pixelSizeGT", 10), 0);

        gd.addNumericField("Reference image pixel size (nm)", getPrefs("pixelSizeRef", 100), 0);
        gd.addNumericField("Reference PSF FWHM (nm)", getPrefs("refFWHM", 300), 0);

        gd.addNumericField("SR image pixel size (nm)", getPrefs("pixelSizeSR", 20), 0);
        gd.addNumericField("SR FWHM (nm)", getPrefs("SRFWHM", 50), 0);

    }

    @Override
    public boolean loadSettings() {
        pixelSizeGT = gd.getNextNumber();

        pixelSizeRef = gd.getNextNumber();
        refFWHM = gd.getNextNumber();

        pixelSizeSR = gd.getNextNumber();
        SRFWHM = gd.getNextNumber();

        setPrefs("pixelSizeGT", pixelSizeGT);
        setPrefs("pixelSizeRef", pixelSizeRef);
        setPrefs("refFWHM", refFWHM);
        setPrefs("pixelSizeSR", pixelSizeSR);
        setPrefs("SRFWHM", SRFWHM);

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


        /// case 1 - step function
        //// create kernel
        int wStepKernel = (int) SRFWHMGTPx;
        if(wStepKernel%2==0) wStepKernel+=1;
        float[] stepKernel = new float[wStepKernel*wStepKernel];
        for(int i=0; i<wStepKernel*wStepKernel; i++) stepKernel[i] = 1f/(wStepKernel*wStepKernel);
        FloatProcessor fpStepKernel = new FloatProcessor(wStepKernel, wStepKernel, stepKernel);

        //// do convolution
        FloatProcessor fpStep = new FloatProcessor(wGT, hGT, pixelsGT.clone());
        fpStep = cv.convolve2D(fpStep, fpStepKernel);

        //// resize
        //fpStep = (FloatProcessor) fpStep.resize(wSR, hSR);

        //// normalise image to 1
        float[] pixelsStep = normalise((float[]) fpStep.getPixels(), 1);
        fpStep.setPixels(pixelsStep);

        /// case 2 - gaussian blur
        FloatProcessor fpGauss = new FloatProcessor(wGT, hGT, pixelsGT.clone());
        fpGauss.blurGaussian(SRSigmaGTPx);

        //// resize
        //fpGauss = (FloatProcessor) fpGauss.resize(wSR, hSR);

        //// normalise image to 1
        float[] pixelsGauss = normalise((float[]) fpGauss.getPixels(), 1);
        fpGauss.setPixels(pixelsGauss);

//
//        /// case 3 - triangle
//        //// create kernel
//        float maxKernelVal = 10;
//        int wTriangleKernel = (int) (2*SRFWHMGTPx)+1;
//        float[] triangleKernel = new float[wTriangleKernel];
//        triangleKernel[0] = 1;
//        triangleKernel[wTriangleKernel-1] = 1;
//        triangleKernel[(int) SRFWHMGTPx] = maxKernelVal;
//
//        float grad = (float) ((maxKernelVal-1)/(floor(wTriangleKernel/2)));
//        for(int i=1; i<(int) SRFWHMGTPx; i++) triangleKernel[i] = i*grad + 1;
//        for(int i=(int) SRFWHMGTPx + 1; i<wTriangleKernel; i++) triangleKernel[i] = -i*grad + 1;
//        //TODO: note - not normalised
//
//        FloatProcessor fpTriangleKernel = new FloatProcessor(triangleKernel.length, 1, triangleKernel);
//
//        //// horizontal convolution
//        FloatProcessor fpTriangle = new FloatProcessor(wGT, hGT, pixelsGT);
//        cv.convolve2D(fpTriangle, fpTriangleKernel);
//
//        //// vertical convolution
//        fpTriangle.rotateRight();
//        cv.convolve2D(fpTriangle, fpTriangleKernel);
//        fpTriangle.rotateLeft();
//
//        //// resize
//        //fpTriangle = (FloatProcessor) fpTriangle.resize(wSR, hSR);
//
//        //// normalise image to 1
//        float[] pixelsTriangle = normalise((float[]) fpTriangle.getPixels(), 1);
//        fpTriangle.setPixels(pixelsTriangle);

        /// case 3 - squared gaussian
        float[] pixelsGauss2 = new float[pixelsGauss.length];
        for(int i=0; i<pixelsGauss2.length; i++) pixelsGauss2[i] = pixelsGauss[i]*pixelsGauss[i];

        //// normalise image to 1
        pixelsGauss2 = normalise(pixelsGauss2, 1);
        FloatProcessor fpGauss2 = new FloatProcessor(wGT, hGT, pixelsGauss2);


        /// case 4 - shadowed gauss
        //// create blurred images
        FloatProcessor fpLargeGauss = new FloatProcessor(wGT, hGT, pixelsGT.clone());
        fpLargeGauss.blurGaussian(SRSigmaGTPx*2);
        FloatProcessor fpShadow = new FloatProcessor(wGT, hGT, pixelsGT.clone());
        fpShadow.blurGaussian(SRSigmaGTPx);

        //// perform subtraction
        float[] pixelsShadow = (float[]) fpShadow.getPixels();
        float[] pixelsLargeGauss = (float[]) fpLargeGauss.getPixels();
        for(int i=0; i<nPixelsGT; i++) pixelsShadow[i] -= pixelsLargeGauss[i];

        //// resize
        //fpShadow.setPixels(pixelsShadow);
        //fpShadow = (FloatProcessor) fpShadow.resize(wSR, hSR);

        //// normalise image to 1
        pixelsShadow = normalise((float[]) fpShadow.getPixels(), 1);
        fpShadow.setPixels(pixelsShadow);

        /// Create image stack

        ImageStack imsSR = new ImageStack(wSR, hSR, 5);
        fpGT.setInterpolationMethod(ImageProcessor.NONE);
        imsSR.setProcessor(fpGT.duplicate().resize(wSR, hSR, true), 1);
        imsSR.setSliceLabel("GT", 1);
        fpStep.setInterpolationMethod(ImageProcessor.NONE);
        imsSR.setProcessor(fpStep.resize(wSR, hSR, true), 2);
        imsSR.setSliceLabel("Square function", 2);
        fpGauss.setInterpolationMethod(ImageProcessor.NONE);
        imsSR.setProcessor(fpGauss.resize(wSR, hSR, true), 3);
        imsSR.setSliceLabel("Gaussian blur", 3);
        fpGauss2.setInterpolationMethod(ImageProcessor.NONE);
        imsSR.setProcessor(fpGauss2.resize(wSR, hSR, true), 4);
        imsSR.setSliceLabel("Squared Gaussian function", 4);
        fpShadow.setInterpolationMethod(ImageProcessor.NONE);
        imsSR.setProcessor(fpShadow.resize(wSR, hSR, true), 5);
        imsSR.setSliceLabel("Shadow Gauss", 5);

        new ImagePlus("SR images", imsSR).show();

    }

    float[] normalise(float[] pixels, float targetMax){
        float maxVal = 0;
        for(int i=0; i<pixels.length; i++) maxVal = max(maxVal, pixels[i]);
        for(int i=0; i<pixels.length; i++) pixels[i] *= (targetMax/maxVal);
        return pixels;
    }
}
