package nanoj.squirrel.java.gui.testing;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import nanoj.kernels.Kernel_SquirrelSwarmOptimizer;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.awt.*;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import static java.lang.Math.*;
import static nanoj.kernels.Kernel_BasePSO.*;

public class ImageSizeTest_ extends _BaseSQUIRRELDialog_ {

    String[] structureOptions = new String[]{"0-D (dots)", "1-D (lines)", "2-D (patches)"};
    double maxAlpha, maxBeta, maxSigma;
    double minAlpha, minBeta, minSigma;
    int nRepeats;
    boolean addNoise;

    int magnification=5;

    Random random = new Random();

    private static Kernel_SquirrelSwarmOptimizer kPSOErrorMap = new Kernel_SquirrelSwarmOptimizer();

    @Override
    public boolean beforeSetupDialog(String s) {
        autoOpenImp = false;
        useSettingsObserver = true;
        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Simulation situation");

        gd.addNumericField("Minimum α", getPrefs("minAlpha", 100), 1);
        gd.addNumericField("Maximum α", getPrefs("maxAlpha", 10000), 0);

        gd.addNumericField("Minimum β", getPrefs("minBeta", 50), 1);
        gd.addNumericField("Maximum β", getPrefs("maxBeta", 5000), 0);

        gd.addNumericField("Minimum σ", getPrefs("minSigma", 1), 1);
        gd.addNumericField("Maximum σ", getPrefs("maxSigma", 15), 1);

        gd.addCheckbox("Add noise?", getPrefs("addNoise", false));

        gd.addNumericField("Number of trials", getPrefs("nRepeats", 100), 0);

    }

    @Override
    public boolean loadSettings() {

        minAlpha = gd.getNextNumber();
        maxAlpha = gd.getNextNumber();
        minBeta = gd.getNextNumber();
        maxBeta = gd.getNextNumber();
        minSigma = gd.getNextNumber();
        maxSigma = gd.getNextNumber();

        addNoise = gd.getNextBoolean();
        nRepeats = (int) gd.getNextNumber();

        setPrefs("minAlpha", minAlpha);
        setPrefs("maxAlpha", maxAlpha);
        setPrefs("minBeta", minBeta);
        setPrefs("maxBeta", maxBeta);
        setPrefs("minSigma", minSigma);
        setPrefs("maxSigma", maxSigma);
        setPrefs("addNoise", addNoise);
        setPrefs("nRepeats", nRepeats);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        int[] imSizeRef = new int[]{64, 128, 256, 512, 1048};
        // populate results table
        ResultsTable rt = new ResultsTable();

        // let's go
        for(int n=0; n<5; n++){

            int w_Ref = imSizeRef[n];
            int h_Ref = imSizeRef[n];
            int w_SR = w_Ref*magnification;
            int h_SR = h_Ref*magnification;

            String[] structure = new String[nRepeats];
            double[] coverage = new double[nRepeats];

            double[] UNIQORNTimes = new double[nRepeats];
            double[] SQUIRRELTimes = new double[nRepeats];

            double[] noiseVars = new double[nRepeats];

            for(int j=0; j<nRepeats; j++) {
                //decide on structure
                int randStruct = random.nextInt(3);
                structure[j] = structureOptions[randStruct];
                // decide on coverage
                double thisCoverage = random.nextDouble() * 0.24 + 0.01;
                coverage[j] = thisCoverage;
                double target = thisCoverage * w_SR * h_SR;

                // generate image
                FloatProcessor fpSR;

                if (randStruct == 0) {
                    fpSR = getDotImage(target, w_SR, h_SR);
                } else if (randStruct == 1) {
                    fpSR = getLineImage(target, w_SR, h_SR);
                } else {
                    fpSR = getPatchImage(target, w_SR, h_SR);
                }

                fpSR.multiply(1000);

                // get parameters for this run
                double thisAlpha = random.nextDouble() * (maxAlpha - minAlpha) + minAlpha;
                double thisBeta = random.nextDouble() * (maxBeta - minBeta) + minBeta;
                double thisSigma = random.nextDouble() * (maxSigma - minSigma) + minSigma;

                // generate reference image
                FloatProcessor fpRef = (FloatProcessor) fpSR.duplicate();
                fpRef.multiply(thisAlpha);
                fpRef.add(thisBeta);
                fpRef.blurGaussian(thisSigma);
                fpRef = (FloatProcessor) fpRef.resize(w_Ref, h_Ref);

                if (addNoise) {
                    float[] pixelsRef = (float[]) fpRef.getPixels();
                    float mean = getMean(pixelsRef);
                    double noiseFraction = random.nextDouble() * 0.25;
                    double noiseVar = mean * noiseFraction;
                    noiseVars[j] = noiseVar;
                    for (int p = 0; p < pixelsRef.length; p++) {
                        pixelsRef[p] = (float) max(0, pixelsRef[p] + random.nextGaussian() * noiseVar);
                    }
                    fpRef.setPixels(pixelsRef);
                }

                // PSO
                long SQUIRRELStart = System.currentTimeMillis();

                kPSOErrorMap.maxMagnification = magnification;
                double sigmaGuess = 0;
                double[] results;
                results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
                        new double[]{0.001, -1000, 0.5}, // low boundary
                        new double[]{100, 1000, 10}, // high boundary
                        new double[]{10, 0, sigmaGuess == 0 ? 5 : sigmaGuess}, // best guess
                        new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, sigmaGuess == 0 ? OBEY_BOUNDARY : CONSTANT}, // impose boundaries
                        1e-7);

                long SQUIRRELStop = System.currentTimeMillis();

                SQUIRRELTimes[j] = SQUIRRELStop - SQUIRRELStart;

                // run UNIQORN
                long UNIQORNStart = System.currentTimeMillis();

                float[] pixelsRef = (float[]) fpRef.duplicate().convertToFloatProcessor().getPixels();
                int nPixelsSR = fpSR.getPixelCount();

                float[] ones = new float[nPixelsSR];
                for (int p = 0; p < nPixelsSR; p++) {
                    ones[p] = 1;
                }

                // UNIVARIATE OPTIMIZER - LINEAR MATCHING
                /// setup optimizer
                UnivariateFunction f = new ImageSizeTest_.sigmaOptimiseFunction(fpSR.duplicate().convertToFloatProcessor(), pixelsRef, ones,
                        w_SR, h_SR, w_Ref, h_Ref);
                UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
                /// run optimizer
                UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                        new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, 5 * magnification)); //NYQUIST NOT ASSUMED; limit to 5*magnification
                float sigma_linear = (float) result.getPoint();


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

                long UNIQORNStop = System.currentTimeMillis();

                UNIQORNTimes[j] = UNIQORNStop - UNIQORNStart;
            }

            for (int i = 0; i < nRepeats; i++) {
                rt.incrementCounter();
                rt.addValue("SR image size", w_SR);
                rt.addValue("Trial number", i);
                rt.addValue("Structure type", structure[i]);
                rt.addValue("% coverage", coverage[i] * 100);
                if(addNoise){
                    rt.addValue("Noise variance", noiseVars[i]);
                }
                rt.addValue("SQUIRREL runtime", SQUIRRELTimes[i]);
                rt.addValue("UNIQORN! runtime", UNIQORNTimes[i]);
            }
        }

        rt.show("Results");
    }

    FloatProcessor getDotImage(double targetCoverage, int w_SR, int h_SR){
        float[] pixels = new float[w_SR*h_SR];
        FloatProcessor fp = new FloatProcessor(w_SR, h_SR, pixels);
        int structurePixels = 0;

        while(structurePixels<targetCoverage){
            int x = (int) Math.floor(random.nextDouble()*w_SR);
            int y = (int) Math.floor(random.nextDouble()*h_SR);
            fp.setf(x, y, (float) max(0.01, random.nextGaussian()+1.0));
            structurePixels = getSum(fp);
        }

        return fp;
    }

    FloatProcessor getLineImage(double targetCoverage, int w_SR, int h_SR){
        float[] pixels = new float[w_SR*h_SR];
        FloatProcessor fp = new FloatProcessor(w_SR, h_SR, pixels);
        int structurePixels = 0;

        while(structurePixels<targetCoverage) {
            int x1 = (int) Math.floor(random.nextDouble() * w_SR);
            int x2 = (int) Math.floor(random.nextDouble() * w_SR);
            int y1 = (int) Math.floor(random.nextDouble() * h_SR);
            int y2 = (int) Math.floor(random.nextDouble() * h_SR);
            fp.setColor(max(0.01,random.nextGaussian() + 1));
            fp.drawLine(x1, y1, x2, y2);

            structurePixels = getSum(fp);
        }
        return fp;
    }

    FloatProcessor getPatchImage(double targetCoverage, int w_SR, int h_SR){
        float[] pixels = new float[w_SR*h_SR];
        FloatProcessor fp = new FloatProcessor(w_SR, h_SR, pixels);
        int structurePixels = 0;

        while(structurePixels<targetCoverage){
            int x0 = (int) Math.floor(random.nextDouble()*w_SR);
            int y0 = (int) Math.floor(random.nextDouble()*h_SR);

            int nVertices = (int) Math.floor(random.nextDouble()*9)+3;
            double maxRadius = random.nextDouble()*w_SR*0.05;

            double[] angles = new double[nVertices];
            for(int i=0; i<nVertices; i++){
                angles[i] = random.nextDouble()*2*PI;
            }
            angles = doSort(angles);
            int[] xVals = new int[nVertices];
            int[] yVals = new int[nVertices];

            for(int i=0; i<nVertices; i++){
                double thisRadius = maxRadius*random.nextDouble();
                double thisAngle = angles[i];

                int y = (int) (thisRadius*cos(thisAngle)) + y0;
                int x = (int) (thisRadius*sin(thisAngle)) + x0;

                xVals[i] = x;
                yVals[i] = y;
            }

            Polygon p = new Polygon(xVals, yVals, nVertices);
            fp.setColor(max(0.01,random.nextGaussian()+1.0));
            fp.fillPolygon(p);

            structurePixels = getSum(fp);
        }

        return fp;
    }

    int getSum(FloatProcessor fp){
        float[] pixels = (float[]) fp.getPixels();
        float v=0;
        for(int i=0; i<pixels.length; i++){
            if(pixels[i]>0) v++;;
        }
        return (int) v;
    }

    double[] doSort(double[] array){
        int i = 1;
        double x;
        int j;

        while(i < array.length){
            x = array[i];
            j = i-1;
            while(j>=0 && array[j] > x){
                array[j+1] = array[j];
                j--;
            }
            array[j+1] = x;
            i++;
        }
        return array;
    }

    private class sigmaOptimiseFunction implements UnivariateFunction {

        FloatProcessor fpSR;
        float[] pixelsRef, ones;
        int w_SR, h_SR, w_Ref, h_Ref;
        ArrayList<Float> sigmaList = new ArrayList<Float>();
        ArrayList<Float> errorList = new ArrayList<Float>();

        public sigmaOptimiseFunction(FloatProcessor fpSR, float[] pixelsRef, float[] ones, int w_SR, int h_SR, int w_Ref, int h_Ref){
            this.fpSR = fpSR;
            this.pixelsRef = pixelsRef;
            this.ones = ones;
            this.w_SR = w_SR;
            this.h_SR = h_SR;
            this.w_Ref = w_Ref;
            this.h_Ref = h_Ref;
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
            return error;
        }

    }

    public float[] calculateAlphaBeta(float[] xA, float[] y, float[] oneA){
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

    private float getMean(float[] array){
        int N = array.length;
        float mean = 0;
        for(int i=0; i<N; i++){
            mean += array[i]/N;
        }
        return mean;
    }

}
