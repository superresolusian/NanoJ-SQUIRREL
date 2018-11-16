package nanoj.squirrel.java.gui.testing;

import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import nanoj.core.java.gui._BaseDialog_;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class BasicTesting_ extends _BaseDialog_ {

    ImagePlus impRef, impSR;
    int w_SR, h_SR, w_Ref, h_Ref, nPixelsSR;

    double trueAlpha, trueBeta, trueSigma;
    double maxAlpha, maxBeta, maxSigma;
    double minAlpha, minBeta, minSigma;
    int magnification, nRepeats;

    Random random = new Random();

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        useSettingsObserver = true;

        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Test accuracy of α, β and σ");

        gd.addNumericField("Minimum α", getPrefs("minAlpha", 100), 1);
        gd.addNumericField("Maximum α", getPrefs("maxAlpha", 10000), 0);

        gd.addNumericField("Minimum β", getPrefs("minBeta", 50), 1);
        gd.addNumericField("Maximum β", getPrefs("maxBeta", 5000), 0);

        gd.addNumericField("Minimum σ", getPrefs("minSigma", 1), 1);
        gd.addNumericField("Maximum σ", getPrefs("maxSigma", 15), 1);

        gd.addNumericField("Magnification", getPrefs("magnification", 5), 0);

        gd.addNumericField("Number of repeats", getPrefs("nRepeats", 100), 0);

    }

    @Override
    public boolean loadSettings() {
        minAlpha = gd.getNextNumber();
        maxAlpha = gd.getNextNumber();
        minBeta = gd.getNextNumber();
        maxBeta = gd.getNextNumber();
        minSigma = gd.getNextNumber();
        maxSigma = gd.getNextNumber();

        magnification = (int) gd.getNextNumber();
        nRepeats = (int) gd.getNextNumber();

        setPrefs("minAlpha", minAlpha);
        setPrefs("maxAlpha", maxAlpha);
        setPrefs("minBeta", minBeta);
        setPrefs("maxBeta", maxBeta);
        setPrefs("minSigma", minSigma);
        setPrefs("maxSigma", maxSigma);
        setPrefs("magnification", magnification);
        setPrefs("nRepeats", nRepeats);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        FloatProcessor fpSR = (FloatProcessor) imp.getProcessor();

        w_SR = fpSR.getWidth();
        h_SR = fpSR.getHeight();
        nPixelsSR = w_SR * h_SR;

        //set up arrays to store real and estimated parameters
        double[] trueAlphas = new double[nRepeats];
        double[] trueBetas = new double[nRepeats];
        double[] trueSigmas = new double[nRepeats];

        double[] estimatedAlphas = new double[nRepeats];
        double[] estimatedBetas = new double[nRepeats];
        double[] estimatedSigmas = new double[nRepeats];


        for(int n=0; n<nRepeats; n++){
            log.status("Trial "+(n+1)+"/"+nRepeats);
            log.progress(n, nRepeats);
            //get parameters for this run
            trueAlpha = random.nextDouble()*(maxAlpha-minAlpha) + minAlpha;
            trueBeta = random.nextDouble()*(maxBeta-minBeta) + minBeta;
            trueSigma = random.nextDouble()*(maxSigma-minSigma) + minSigma;

            trueAlphas[n] = trueAlpha;
            trueBetas[n] = trueBeta;
            trueSigmas[n] = trueSigma;

            //generate reference image
            FloatProcessor fpRef = (FloatProcessor) fpSR.duplicate();
            fpRef.multiply(trueAlpha);
            fpRef.add(trueBeta);
            fpRef.blurGaussian(trueSigma);
            fpRef = (FloatProcessor) fpRef.resize(w_SR/magnification, h_SR/magnification);

            w_Ref = fpRef.getWidth();
            h_Ref = fpRef.getHeight();


            //prepare arrays for optimisation

            /// images
            float[] pixelsRef = (float[]) fpRef.getPixels();

            /// ones matrix
            float[] ones = new float[nPixelsSR];
            for(int i=0; i<nPixelsSR; i++){ones[i] = 1;}

            // UNIVARIATE OPTIMISER
            //setup optimizer
            sigmaOptimiseFunction f = new BasicTesting_.sigmaOptimiseFunction(fpSR, pixelsRef, ones);
            UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
            UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                    new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, 20));
            double estimatedSigma = result.getPoint();

            // GET ALPHA AND BETA
            FloatProcessor blurredFp = (FloatProcessor) fpSR.duplicate();
            FloatProcessor blurredOnes = new FloatProcessor(w_SR, h_SR, ones);
            blurredFp.blurGaussian(result.getPoint());
            blurredOnes.blurGaussian(result.getPoint());

            blurredFp = (FloatProcessor) blurredFp.resize(w_Ref, h_Ref);
            blurredOnes = (FloatProcessor) blurredOnes.resize(w_Ref, h_Ref);

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRef, (float[]) blurredOnes.getPixels());
            double estimatedAlpha = aB[0];
            double estimatedBeta = aB[1];

            estimatedAlphas[n] = estimatedAlpha;
            estimatedBetas[n] = estimatedBeta;
            estimatedSigmas[n] = estimatedSigma;
        }

        //calculate percentage errors on estimated parameters
        double[] averagePercentageErrors = new double[3], stdevPercentageErrors = new double[3];

        for(int n=0; n<nRepeats; n++){
            double thisErrorAlpha = (trueAlphas[n]-estimatedAlphas[n])/trueAlphas[n];
            double thisErrorBeta = (trueBetas[n]-estimatedBetas[n])/trueBetas[n];
            double thisErrorSigma = (trueSigmas[n]-estimatedSigmas[n])/trueSigmas[n];
            if(thisErrorAlpha>0.1 || thisErrorBeta>0.1 || thisErrorSigma>0.1){
                log.msg("POOR SOLUTION AT DATAPOINT "+(n+1)+": alpha= "+trueAlphas[n]+", beta= "+trueBetas[n]+", sigma = "+trueSigmas[n]);
            }
            averagePercentageErrors[0] += thisErrorAlpha/nRepeats;
            averagePercentageErrors[1] += thisErrorBeta/nRepeats;
            averagePercentageErrors[2] += thisErrorSigma/nRepeats;
        }

        for(int n=0; n<nRepeats; n++){
            double thisDifferenceAlpha = Math.pow((trueAlphas[n]-estimatedAlphas[n])/trueAlphas[n] - averagePercentageErrors[0],2);
            double thisDifferenceBeta = Math.pow((trueBetas[n]-estimatedBetas[n])/trueBetas[n] - averagePercentageErrors[1],2);
            double thisDifferenceSigma = Math.pow((trueSigmas[n]-estimatedSigmas[n])/trueSigmas[n] - averagePercentageErrors[2],2);
            stdevPercentageErrors[0] += thisDifferenceAlpha/nRepeats;
            stdevPercentageErrors[1] += thisDifferenceBeta/nRepeats;
            stdevPercentageErrors[2] += thisDifferenceSigma/nRepeats;
        }
        stdevPercentageErrors[0] = Math.sqrt(stdevPercentageErrors[0]);
        stdevPercentageErrors[1] = Math.sqrt(stdevPercentageErrors[1]);
        stdevPercentageErrors[2] = Math.sqrt(stdevPercentageErrors[2]);

        //plot
        Plot alphaPlot = new Plot("Estimated vs true α", "True α", "Estimated α");
        alphaPlot.addPoints(trueAlphas, estimatedAlphas, Plot.CIRCLE);
        alphaPlot.show();
        Plot betaPlot = new Plot("Estimated vs true β", "True β", "Estimated β");
        betaPlot.addPoints(trueBetas, estimatedBetas, Plot.CIRCLE);
        betaPlot.show();
        Plot sigmaPlot = new Plot("Estimated vs true σ", "True σ", "Estimated σ");
        sigmaPlot.addPoints(trueSigmas, estimatedSigmas, Plot.CIRCLE);
        sigmaPlot.show();

        ResultsTable rt = new ResultsTable();
        String[] labels = new String[]{"α", "β", "σ"};
        for(int i=0; i<3; i++){
            rt.incrementCounter();
            rt.addValue("Parameter", labels[i]);
            rt.addValue("Average % error", averagePercentageErrors[i]*100);
            rt.addValue("Stdev % error", stdevPercentageErrors[i]*100);
        }
        rt.show("Percentage errors (n="+nRepeats+")");

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
