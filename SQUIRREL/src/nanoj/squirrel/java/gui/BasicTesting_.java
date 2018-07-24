package nanoj.squirrel.java.gui;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.gui._BaseDialog_;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import static java.lang.Math.abs;
import static nanoj.core.java.array.ArrayCasting.toArray;

public class BasicTesting_ extends _BaseDialog_ {

    ImagePlus impRef, impSR;
    int w_SR, h_SR, w_Ref, h_Ref, nPixelsRef, nPixelsSR, nSlicesSR, magnification;

    double trueAlpha, trueBeta, trueSigma;
    double maxAlpha=10000, maxBeta=5000, maxSigma=15;
    double minAlpha=100, minBeta =50, minSigma=1;

    Random random = new Random();

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        useSettingsObserver = true;

        return true;
    }

    @Override
    public void setupDialog() { }

    @Override
    public boolean loadSettings() {
        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        FloatProcessor fpSR = (FloatProcessor) imp.getProcessor();

        w_SR = fpSR.getWidth();
        h_SR = fpSR.getHeight();
        nPixelsSR = w_SR * h_SR;

        magnification = 5;

        trueAlpha = random.nextDouble()*(maxAlpha-minAlpha) + minAlpha;
        trueBeta = random.nextDouble()*(maxBeta-minBeta) + minBeta;
        //trueSigma = random.nextDouble()*(maxSigma-minSigma) + minSigma;
        trueSigma = 5;

        log.msg("-----Simulated GT parameters-----");
        log.msg("a = "+trueAlpha+", B = "+trueBeta+", sig ="+trueSigma);
        log.msg("---------------------------------");

        FloatProcessor fpRef = (FloatProcessor) fpSR.duplicate();
        fpRef.multiply(trueAlpha);
        fpRef.add(trueBeta);
        fpRef.blurGaussian(trueSigma);
        fpRef = (FloatProcessor) fpRef.resize(w_SR/magnification, h_SR/magnification);

        impRef = new ImagePlus("Widefield", fpRef);
        impRef.show();

        w_Ref = fpRef.getWidth();
        h_Ref = fpRef.getHeight();


        long startTime = System.currentTimeMillis();

        // PREPARATION
        /// images
        float[] pixelsRef = (float[]) fpRef.getPixels();

        /// ones matrix
        float[] ones = new float[nPixelsSR];
        for(int i=0; i<nPixelsSR; i++){ones[i] = 1;}

        // BRUTE FORCE
        float sigmaStep = 0.1f;
        float sigmaRange = (float) (maxSigma-minSigma);
        int nSteps = (int) (sigmaRange/sigmaStep);

        float[] sigmas = new float[nSteps];
        float[] errors = new float[nSteps];
        float[] alphas = new float[nSteps];
        float[] betas = new float[nSteps];

        ImageStack imsConvolved = new ImageStack(w_SR, h_SR, nSteps);

        for(int n=0; n<nSteps; n++){
            float sigma = (float) (minSigma + n*sigmaStep);

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

            imsConvolved.setProcessor(finalFpSR, n+1);

            float error = calculateRMSE(pixelsRef, (float[]) finalFpSRResized.getPixels());

            sigmas[n] = sigma;
            errors[n] = error;
            alphas[n] = aB[0];
            betas[n] = aB[1];
        }

        new Plot("Error vs sigma", "Sigma", "RMSE", sigmas, errors).show();
        new Plot("Alpha vs sigma", "Sigma", "Alpha", sigmas, alphas).show();
        new Plot("Beta vs sigma", "Sigma", "Beta", sigmas, betas).show();


        // UNIVARIATE OPTIMISER
        //setup optimizer
        sigmaOptimiseFunction f = new BasicTesting_.sigmaOptimiseFunction(fpSR, pixelsRef, ones);
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
        UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, 20));
        log.msg("Best sigma is: " + result.getPoint());

        float[] errorList = toArray(f.getErrorList(), 1.0f);
        float[] sigmaList = toArray(f.getSigmaList(), 1.0f);

        new Plot("Optimiser: error vs sigma", "Sigma", "RMSE", sigmaList, errorList).show();

        // GET ALPHA AND BETA
        FloatProcessor blurredFp = (FloatProcessor) fpSR.duplicate();
        FloatProcessor blurredOnes = new FloatProcessor(w_SR, h_SR, ones);
        blurredFp.blurGaussian(result.getPoint());
        blurredOnes.blurGaussian(result.getPoint());

        blurredFp = (FloatProcessor) blurredFp.resize(w_Ref, h_Ref);
        blurredOnes = (FloatProcessor) blurredOnes.resize(w_Ref, h_Ref);

        float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRef, (float[]) blurredOnes.getPixels());
        float alpha = aB[0];
        float beta = aB[1];

        log.msg("Alpha is: " + alpha + ", beta is: " + beta);

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

}
