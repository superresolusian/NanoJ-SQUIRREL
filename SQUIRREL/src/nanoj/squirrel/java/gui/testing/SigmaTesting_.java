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

public class SigmaTesting_ extends _BaseDialog_ {

    ImagePlus impRef, impSR;
    int w_SR, h_SR, w_Ref, h_Ref, nPixelsSR;

    double trueSigma;
    double alpha, beta;
    double maxSigma, minSigma;
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
        gd = new NonBlockingGenericDialog("Unit test for σ");

        gd.addNumericField("Fixed α", getPrefs("alpha", 100), 1);

        gd.addNumericField("Fixed β", getPrefs("beta", 50), 1);

        gd.addNumericField("Minimum σ", getPrefs("minSigma", 1), 0);
        gd.addNumericField("Maximum σ", getPrefs("maxSigma", 15), 0);

        gd.addNumericField("Magnification", getPrefs("magnification", 5), 0);

        gd.addNumericField("Number of repeats", getPrefs("nRepeats", 100), 0);

    }

    @Override
    public boolean loadSettings() {
        alpha = gd.getNextNumber();
        beta = gd.getNextNumber();
        minSigma = gd.getNextNumber();
        maxSigma = gd.getNextNumber();

        magnification = (int) gd.getNextNumber();
        nRepeats = (int) gd.getNextNumber();

        setPrefs("alpha", alpha);
        setPrefs("beta", beta);
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
        double[] trueSigmas = new double[nRepeats];
        double[] estimatedSigmas = new double[nRepeats];


        for(int n=0; n<nRepeats; n++){
            log.status("Trial "+(n+1)+"/"+nRepeats);
            log.progress(n, nRepeats);
            //get parameters for this run

            trueSigma = random.nextDouble()*(maxSigma-minSigma) + minSigma;

            trueSigmas[n] = trueSigma;

            //generate reference image
            FloatProcessor fpRef = (FloatProcessor) fpSR.duplicate();
            fpRef.multiply(alpha);
            fpRef.add(beta);
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
            sigmaOptimiseFunction f = new SigmaTesting_.sigmaOptimiseFunction(fpSR, pixelsRef, ones);
            UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
            UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                    new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, 20));
            double estimatedSigma = result.getPoint();

            estimatedSigmas[n] = estimatedSigma;
        }

        //calculate percentage errors on estimated parameters
        double[] averagePercentageErrors = new double[1], stdevPercentageErrors = new double[1];

        for(int n=0; n<nRepeats; n++){
            double thisErrorSigma = (trueSigmas[n]-estimatedSigmas[n])/trueSigmas[n];
            if(thisErrorSigma>0.1){
                log.msg("POOR SOLUTION AT DATAPOINT "+(n+1)+": sigma = "+trueSigmas[n]);
            }
            averagePercentageErrors[0] += thisErrorSigma/nRepeats;
        }

        for(int n=0; n<nRepeats; n++){
            double thisDifferenceSigma = Math.pow((trueSigmas[n]-estimatedSigmas[n])/trueSigmas[n] - averagePercentageErrors[0],2);
            stdevPercentageErrors[0] += thisDifferenceSigma/nRepeats;
        }
        stdevPercentageErrors[0] = Math.sqrt(stdevPercentageErrors[0]);

        //plot
        Plot sigmaPlot = new Plot("Sigma only - Estimated vs true σ", "True σ", "Estimated σ");
        sigmaPlot.addPoints(trueSigmas, estimatedSigmas, Plot.CIRCLE);
        sigmaPlot.show();

        ResultsTable rt = new ResultsTable();
        String[] labels = new String[]{"σ"};
        for(int i=0; i<1; i++){
            rt.incrementCounter();
            rt.addValue("Parameter", labels[i]);
            rt.addValue("Average % error", averagePercentageErrors[i]*100);
            rt.addValue("Stdev % error", stdevPercentageErrors[i]*100);
        }
        rt.show("Percentage errors (n="+nRepeats+")");

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

            float[] aB = new float[] {(float)alpha, (float)beta};

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
    }

}
