package nanoj.squirrel.java.gui.testing;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import nanoj.core.java.gui._BaseDialog_;
import nanoj.core.java.threading.NanoJThreadExecutor;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by sculley on 18/01/2019.
 */
public class RollingBlockTest_ extends _BaseDialog_ {

    ArrayList<String> titles = new ArrayList<String>();
    String[] imageTitles;

    String titleRefImage = "", titleSRImage = "";

    ImagePlus impRef, impSR;
    ImageStack imsRef, imsSR;
    int nSlicesRef, nSlicesSR;
    int w_SR, h_SR, w_Ref, h_Ref;
    int magnification;

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
        for(int n=0; n<nImages; n++){
            ImagePlus thisImp = WindowManager.getImage(imageTitles[n]);
            titles.add(thisImp.getTitle());
        }

        imageTitles = new String[titles.size()];
        for(int n=0; n<titles.size(); n++){
            imageTitles[n] = titles.get(n);
        }

        return true;
    }

    @Override
    public void setupDialog() {

        gd = new NonBlockingGenericDialog("Test blocking");

        if (titleRefImage.equals("")) {
            titleRefImage = getPrefs("titleRefImage", "");
            titleSRImage = getPrefs("titleSRImage", "");
        }

        if (!titleRefImage.equals("") && contains(titleRefImage, imageTitles))
            gd.addChoice("Reference Image", imageTitles, titleRefImage);
        else
            gd.addChoice("Reference Image", imageTitles, imageTitles[0]);

        if (!titleSRImage.equals("") && contains(titleSRImage, imageTitles))
            gd.addChoice("Super-Resolution Reconstruction(s)", imageTitles, titleSRImage);
        else
            gd.addChoice("Super-Resolution Reconstruction(s)", imageTitles, imageTitles[0]);

    }

    @Override
    public boolean loadSettings() {

        titleRefImage = gd.getNextChoice();
        titleSRImage = gd.getNextChoice();

        setPrefs("titleRefImage", titleRefImage);
        setPrefs("titleSRImage", titleSRImage);

        impRef = WindowManager.getImage(titleRefImage);
        impSR = WindowManager.getImage(titleSRImage);

        imsSR = impSR.getImageStack();
        imsRef = impRef.getImageStack();

        w_SR = imsSR.getWidth();
        h_SR = imsSR.getHeight();

        w_Ref = imsRef.getWidth();
        h_Ref = imsRef.getHeight();

        magnification = w_SR / w_Ref;

        nSlicesSR = impSR.getStackSize();
        nSlicesRef = impRef.getStackSize();

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        FloatProcessor fpSR = imsSR.getProcessor(1).convertToFloatProcessor();
        FloatProcessor fpRef = imsRef.getProcessor(1).convertToFloatProcessor();

        float maxSigmaBoundary = 7.5f;

        int blockSize = (int) (5*maxSigmaBoundary);
        log.msg("blockSize is "+blockSize);
        while(blockSize%magnification!=0){
            blockSize++;
        }
        log.msg("new blockSize is "+blockSize);

        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);

        //ImageStack imsSRBlocks = new ImageStack(blockSize,blockSize,w_SR*h_SR);
        //ImageStack imsRefBlocks = new ImageStack(blockSize/magnification, blockSize/magnification, w_SR*h_SR);

        float[] alphas = new float[w_SR*h_SR];
        float[] betas = new float[w_SR*h_SR];
        float[] sigmas = new float[w_SR*h_SR];

        log.status("doing stuff");
        for (int y = 0; y < h_SR; y++) {
            for (int x = 0; x < w_SR; x++) {
                int ind = y*w_SR+x;
                log.progress(ind/(w_SR*h_SR));
                NTE.execute(new SplitIntoBlocks(fpSR, fpRef, x, y, blockSize, maxSigmaBoundary, alphas, betas, sigmas));
            }
        }

        NTE.finish();

        //new ImagePlus("SR blocks", imsSRBlocks).show();
        //new ImagePlus("Ref blocks", imsRefBlocks).show();

        new ImagePlus("alpha", new FloatProcessor(w_SR, h_SR, alphas)).show();
        new ImagePlus("beta", new FloatProcessor(w_SR, h_SR, betas)).show();
        new ImagePlus("sigma", new FloatProcessor(w_SR, h_SR, sigmas)).show();

    }

    //helper function lifted from imageshelper
    protected static boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
    }

    class SplitIntoBlocks extends Thread{

        private FloatProcessor fpSR, fpRef;
        private int x, y, blockSize;
        private float maxSigmaBoundary;
        private float[] alphas, betas, sigmas;

        public SplitIntoBlocks(FloatProcessor fpSR, FloatProcessor fpRef, int x, int y,
                               int blockSize, float maxSigmaBoundary,
                               float[] alphas, float[] betas, float[] sigmas){
            this.fpSR = fpSR;
            this.fpRef = fpRef;
            this.x = x;
            this.y = y;
            this.blockSize = blockSize;
            this.maxSigmaBoundary = maxSigmaBoundary;
            this.alphas = alphas;
            this.betas = betas;
            this.sigmas = sigmas;
        }

        public void run(){
            int offset = blockSize/2;
            int blockSizeRef = blockSize/magnification;

            float[] pixelsSRBlock = new float[blockSize*blockSize];
            float[] pixelsRefBlock = new float[blockSizeRef*blockSizeRef];

            int p = y*w_SR + x;

            for(int j=0; j<blockSize; j++){
                int y_ = (j+y)-offset;
                if(y_<0 || y_>(h_SR-1)) continue;

                for(int i=0; i<blockSize; i++){
                    int x_ = (i+x)-offset;
                    if(x_<0 || x_>(w_SR-1)) continue;

                    int k = j*blockSize + i;

                    pixelsSRBlock[k] = fpSR.getf(x_,y_);

                    int jRef = j/magnification;
                    int iRef = i/magnification;
                    int yRef_ = y_/magnification;
                    int xRef_ = x_/magnification;

                    int kRef = jRef*blockSizeRef + iRef;

                    pixelsRefBlock[kRef] = fpRef.getf(xRef_, yRef_);
                }
            }

            // work out borders
            int offsetRef = offset/magnification;
            int xBorderStartRef = Math.max((offsetRef-x/magnification), 0);
            int xBorderEndRef = Math.max((x/magnification-offsetRef+blockSizeRef-1)-(w_Ref-1),0);
            int yBorderStartRef = Math.max((offsetRef-y/magnification), 0);
            int yBorderEndRef = Math.max((y/magnification-offsetRef+blockSizeRef-1)-(h_Ref-1),0);

            FloatProcessor fpSRBlock = new FloatProcessor(blockSize, blockSize, pixelsSRBlock);
            FloatProcessor fpRefBlock = new FloatProcessor(blockSizeRef, blockSizeRef, pixelsRefBlock);

            if(xBorderEndRef>0 || xBorderStartRef>0 || yBorderStartRef>0 || yBorderEndRef>0){
                Rectangle rectangleRef = new Rectangle(xBorderStartRef, yBorderStartRef,
                        blockSizeRef-xBorderEndRef-xBorderStartRef, blockSizeRef-yBorderEndRef-yBorderStartRef);
                fpRefBlock.setRoi(rectangleRef);
                fpRefBlock = fpRefBlock.crop().convertToFloatProcessor();

                Rectangle rectangleSR = new Rectangle(xBorderEndRef*magnification, yBorderStartRef*magnification,
                        (blockSizeRef-xBorderEndRef-xBorderStartRef)*magnification,
                        (blockSizeRef-yBorderEndRef-yBorderStartRef)*magnification);
                fpSRBlock.setRoi(rectangleSR);
                fpSRBlock = fpSRBlock.crop().convertToFloatProcessor();
            }

            float[] ones = new float[fpSRBlock.getPixelCount()];
            for(int i=0; i<ones.length; i++) ones[i] = 1;

            // UNIVARIATE OPTIMIZER - LINEAR MATCHING
            /// setup optimizer
            sigmaOptimiseFunction f =  new RollingBlockTest_.sigmaOptimiseFunction(fpSRBlock, fpRefBlock, ones);
            UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);

            /// run optimizer
            UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                    new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, blockSize/2)); //limit to block width
            float sigma_linear = (float) result.getPoint();

            // GET ALPHA AND BETA
            FloatProcessor blurredFp = (FloatProcessor) fpSRBlock.duplicate();
            FloatProcessor blurredOnes = new FloatProcessor(fpSRBlock.getWidth(), fpSRBlock.getHeight(), ones);
            blurredFp.blurGaussian(sigma_linear);
            blurredOnes.blurGaussian(sigma_linear);

            blurredFp = (FloatProcessor) blurredFp.resize(fpRefBlock.getWidth(), fpRefBlock.getHeight());
            blurredOnes = (FloatProcessor) blurredOnes.resize(fpRefBlock.getWidth(), fpRefBlock.getHeight());

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), (float[]) fpRefBlock.getPixels(),
                    (float[]) blurredOnes.getPixels());
            float alpha = aB[0];
            float beta = aB[1];

            // check if sigma hit the boundary

            float alphaBoundary = alpha, betaBoundary = beta, sigmaBoundary = sigma_linear;
            if(sigma_linear>maxSigmaBoundary){
                sigmaBoundary = maxSigmaBoundary;
                // calculate alpha and beta with maxSigmaBoundary as blur
                blurredFp = (FloatProcessor) fpSRBlock.duplicate();
                blurredOnes = new FloatProcessor(fpSRBlock.getWidth(), fpSRBlock.getHeight(), ones);
                blurredFp.blurGaussian(maxSigmaBoundary);
                blurredOnes.blurGaussian(maxSigmaBoundary);

                blurredFp = (FloatProcessor) blurredFp.resize(fpRefBlock.getWidth(), fpRefBlock.getHeight());
                blurredOnes = (FloatProcessor) blurredOnes.resize(fpRefBlock.getWidth(), fpRefBlock.getHeight());

                aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRefBlock, (float[]) blurredOnes.getPixels());
                alphaBoundary = aB[0];
                betaBoundary = aB[1];
            }

            int p_ = y*w_SR + x;

            alphas[p_] = alphaBoundary;
            betas[p_] = betaBoundary;
            sigmas[p_] = sigmaBoundary;

        }



    }

    private class sigmaOptimiseFunction implements UnivariateFunction {

        FloatProcessor fpSRBlock, fpRefBlock;
        float[] ones;
        int blockWidthSR, blockHeightSR, blockWidthRef, blockHeightRef;
        ArrayList<Float> sigmaList = new ArrayList<Float>();
        ArrayList<Float> errorList = new ArrayList<Float>();

        public sigmaOptimiseFunction(FloatProcessor fpSRBlock, FloatProcessor fpRefBlock, float[] ones){
            this.fpSRBlock = fpSRBlock;
            this.fpRefBlock = fpRefBlock;
            this.ones = ones;
            this.blockWidthSR = fpSRBlock.getWidth();
            this.blockHeightSR = fpSRBlock.getHeight();
            this.blockWidthRef = fpRefBlock.getWidth();
            this.blockHeightRef = fpRefBlock.getHeight();
        }

        public double value(double sigma) {

            FloatProcessor blurredFp = (FloatProcessor) fpSRBlock.duplicate();
            FloatProcessor blurredOnes = new FloatProcessor(blockWidthSR, blockHeightSR, ones);
            blurredFp.blurGaussian(sigma);
            blurredOnes.blurGaussian(sigma);

            blurredFp = (FloatProcessor) blurredFp.resize(blockWidthRef, blockHeightRef);
            blurredOnes = (FloatProcessor) blurredOnes.resize(blockWidthRef, blockHeightRef);

            float[] pixelsRefBlock = (float[]) fpRefBlock.getPixels();

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRefBlock, (float[]) blurredOnes.getPixels());

            FloatProcessor finalFpSR = (FloatProcessor) fpSRBlock.duplicate();
            finalFpSR.multiply(aB[0]);
            finalFpSR.add(aB[1]);
            finalFpSR.blurGaussian(sigma);
            FloatProcessor finalFpSRResized = (FloatProcessor) finalFpSR.resize(blockWidthRef, blockHeightRef);

            double error = calculateRMSE(pixelsRefBlock, (float[]) finalFpSRResized.getPixels());
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

}
