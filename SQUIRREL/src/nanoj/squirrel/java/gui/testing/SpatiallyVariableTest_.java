package nanoj.squirrel.java.gui.testing;

import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import nanoj.core.java.gui._BaseDialog_;

import java.awt.*;
import java.io.IOException;
import java.util.Random;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;

/**
 * Created by sculley on 22/01/2019.
 */
public class SpatiallyVariableTest_ extends _BaseDialog_ {

    String[] imageTitles;

    String titleGTImage="", titleAlphaMap="", titleBetaMap="", titleSigmaMap="";
    double meanAlpha, stdAlpha, meanBeta, stdBeta, meanSigma, stdSigma;
    int magnification;

    ImagePlus impGT, impAlpha, impBeta, impSigma;

    int blocksPerXAxis, blocksPerYAxis;

    Random random = new Random();

    @Override
    public boolean beforeSetupDialog(String s) {
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

        gd = new NonBlockingGenericDialog("Simulate reference with spatially inhomogenous alpha, beta, sigma");

        if (titleGTImage.equals("")) {
            titleGTImage = getPrefs("titleGTImage", "");
            titleAlphaMap = getPrefs("titleAlphaMap", "");
            titleBetaMap = getPrefs("titleBetaMap", "");
            titleSigmaMap = getPrefs("titleSigmaMap", "");
        }

        if (!titleGTImage.equals("") && contains(titleGTImage, imageTitles))
            gd.addChoice("Ground Truth Image", imageTitles, titleGTImage);
        else
            gd.addChoice("Ground Truth Image", imageTitles, imageTitles[0]);

        if (!titleAlphaMap.equals("") && contains(titleAlphaMap, imageTitles))
            gd.addChoice("Alpha map", imageTitles, titleAlphaMap);
        else
            gd.addChoice("Alpha map", imageTitles, imageTitles[0]);

        if (!titleBetaMap.equals("") && contains(titleBetaMap, imageTitles))
            gd.addChoice("Beta map", imageTitles, titleBetaMap);
        else
            gd.addChoice("Beta map", imageTitles, imageTitles[0]);

        if (!titleBetaMap.equals("") && contains(titleSigmaMap, imageTitles))
            gd.addChoice("Sigma map (for SR image)", imageTitles, titleSigmaMap);
        else
            gd.addChoice("Sigma map (for SR image)", imageTitles, imageTitles[0]);

        gd.addNumericField("Magnification", getPrefs("magnification", 5), 0);

    }

    @Override
    public boolean loadSettings() {
        titleGTImage = gd.getNextChoice();
        titleAlphaMap = gd.getNextChoice();
        titleBetaMap = gd.getNextChoice();
        titleSigmaMap = gd.getNextChoice();

        magnification = (int) gd.getNextNumber();

        setPrefs("titleGTImage", titleGTImage);
        setPrefs("titleAlphaMap", titleAlphaMap);
        setPrefs("titleBetaMap", titleBetaMap);
        setPrefs("titleSigmaMap", titleSigmaMap);

        setPrefs("magnification", magnification);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {
        impGT = WindowManager.getImage(titleGTImage);
        impAlpha = WindowManager.getImage(titleAlphaMap);
        impBeta = WindowManager.getImage(titleBetaMap);
        impSigma = WindowManager.getImage(titleSigmaMap);

        FloatProcessor fpGT = impGT.getProcessor().convertToFloatProcessor();
        FloatProcessor fpAlpha = impAlpha.getProcessor().convertToFloatProcessor();
        FloatProcessor fpBeta = impBeta.getProcessor().convertToFloatProcessor();
        FloatProcessor fpSigma = impSigma.getProcessor().convertToFloatProcessor();

        int wGT = fpGT.getWidth();
        int hGT = fpGT.getHeight();

        int wRef = wGT/magnification;
        int hRef = hGT/magnification;

        //prepare maps
        fpAlpha = fpAlpha.resize(wGT, hGT).convertToFloatProcessor();
        fpBeta = fpBeta.resize(wGT, hGT).convertToFloatProcessor();
        fpSigma = fpSigma.resize(wGT, hGT).convertToFloatProcessor();

        //generate SR image

        //block-wise convolve image
        int blockSize = 25;
        blocksPerXAxis = (int) floor(wGT / blockSize);
        blocksPerYAxis = (int) floor(hGT / blockSize);

        FloatProcessor fpSRConvolved = new FloatProcessor(wGT, hGT);

        for(int nYB=0; nYB<blocksPerYAxis; nYB++){
            for(int nXB=0; nXB<blocksPerXAxis; nXB++){

                int blockWidth = wGT/blocksPerXAxis;
                int blockHeight = hGT/blocksPerYAxis;

                int xStartSR = nXB*blockWidth;
                int yStartSR = nYB*blockHeight;
                int xCentreSR = xStartSR+blockWidth/2;
                int yCentreSR = yStartSR+blockHeight/2;

                FloatProcessor fp = fpGT.duplicate().convertToFloatProcessor();
                float thisSigma = fpSigma.getf(xCentreSR, yCentreSR);
                fp.blurGaussian(thisSigma);

                fp = getBlockFp(fp, nYB, nXB);

                fpSRConvolved = setBlock(fpSRConvolved, fp, xStartSR, yStartSR, fp.getWidth(), fp.getHeight());
            }
        }

        //intensity-scale image
        FloatProcessor fpIntensityScaled = fpGT.duplicate().convertToFloatProcessor();
        for(int y=0; y<hGT; y++){
            for(int x=0; x<wGT; x++){
                float val = fpIntensityScaled.getf(x, y);
                float alpha = fpAlpha.getf(x, y);
                float beta = fpBeta.getf(x, y);
                fpIntensityScaled.setf(x, y, val*alpha+beta);
            }
        }

        fpIntensityScaled.blurGaussian(7.5);

        FloatProcessor fpRef = fpIntensityScaled.resize(wRef, hRef).convertToFloatProcessor();

        new ImagePlus("Scaled alpha map", fpAlpha).show();
        new ImagePlus("Scaled beta map", fpBeta).show();
        new ImagePlus("Scaled sigma map", fpSigma).show();

        new ImagePlus("SR image", fpSRConvolved).show();

        new ImagePlus("Reference image", fpRef).show();

    }

    protected static boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
    }


    public FloatProcessor scaleMap(FloatProcessor fp, double mean, double sigma){

        float[] pixels = (float[]) fp.getPixels();
        int nPixels = pixels.length;

        float maxVal = Float.MIN_VALUE;
        float minVal = Float.MAX_VALUE;

        for(int i=0; i<nPixels; i++){
            maxVal = max(pixels[i], maxVal);
            minVal = min(pixels[i], minVal);
        }

        for(int i=0; i<nPixels; i++){
            float val = pixels[i];
            val = (val-minVal)/(maxVal-minVal);
            val = (float) (random.nextGaussian()*sigma + mean)*val;
            val = max(val, 0);
            pixels[i] = val;
        }

        fp.setPixels(pixels);

        return fp;

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
}
