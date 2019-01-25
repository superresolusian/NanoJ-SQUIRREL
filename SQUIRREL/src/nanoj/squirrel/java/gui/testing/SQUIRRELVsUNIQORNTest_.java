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

public class SQUIRRELVsUNIQORNTest_ extends _BaseSQUIRRELDialog_ {

    String[] structureOptions = new String[]{"0-D (dots)", "1-D (lines)", "2-D (patches)"};
    double maxAlpha, maxBeta, maxSigma;
    double minAlpha, minBeta, minSigma;
    int nRepeats;
    boolean addNoise;

    int w_SR=500, h_SR=500, w_Ref=100, h_Ref=100, magnification=5;

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

        // set up output images and arrays
        ImageStack imsSR = new ImageStack(w_SR, h_SR, nRepeats);
        ImageStack imsRef = new ImageStack(w_Ref, h_Ref, nRepeats);

        String[] structure = new String[nRepeats];
        double[] coverage = new double[nRepeats];

        double[] trueAlpha = new double[nRepeats];
        double[] trueBeta = new double[nRepeats];
        double[] trueSigma = new double[nRepeats];

        double[] UNIQORNAlpha = new double[nRepeats];
        double[] UNIQORNBeta = new double[nRepeats];
        double[] UNIQORNSigma = new double[nRepeats];

        double[] SQUIRRELAlpha = new double[nRepeats];
        double[] SQUIRRELBeta = new double[nRepeats];
        double[] SQUIRRELSigma = new double[nRepeats];

        double[] UNIQORNTimes = new double[nRepeats];
        double[] SQUIRRELTimes = new double[nRepeats];

        double[] noiseVars = new double[nRepeats];

        // let's go
        for(int i=0; i<nRepeats; i++){

            log.status("Trial "+(i+1)+" - generating structure");

            //decide on structure
            int randStruct = random.nextInt(3);
            structure[i] = structureOptions[randStruct];
            // decide on coverage
            double thisCoverage = random.nextDouble()*0.24 + 0.01;
            coverage[i] = thisCoverage;
            double target = thisCoverage*w_SR*h_SR;

            // generate image
            FloatProcessor fpSR;

            if(randStruct==0){
                fpSR = getDotImage(target);
            }
            else if(randStruct==1){
                fpSR = getLineImage(target);
            }
            else{
                fpSR = getPatchImage(target);
            }

            fpSR.multiply(1000);
            imsSR.setProcessor(fpSR, i+1);

            // get parameters for this run
            double thisAlpha = random.nextDouble()*(maxAlpha-minAlpha) + minAlpha;
            double thisBeta = random.nextDouble()*(maxBeta-minBeta) + minBeta;
            double thisSigma = random.nextDouble()*(maxSigma-minSigma) + minSigma;

            trueAlpha[i] = thisAlpha;
            trueBeta[i] = thisBeta;
            trueSigma[i] = thisSigma;

            // generate reference image
            FloatProcessor fpRef = (FloatProcessor) fpSR.duplicate();
            fpRef.multiply(thisAlpha);
            fpRef.add(thisBeta);
            fpRef.blurGaussian(thisSigma);
            fpRef = (FloatProcessor) fpRef.resize(w_Ref, h_Ref);

            if(addNoise) {
                float[] pixelsRef = (float[]) fpRef.getPixels();
                float mean = getMean(pixelsRef);
                double noiseFraction = random.nextDouble() * 0.25;
                double noiseVar = mean * noiseFraction;
                noiseVars[i] = noiseVar;
                for (int j = 0; j < pixelsRef.length; j++) {
                    pixelsRef[j] = (float) max(0, pixelsRef[j] + random.nextGaussian() * noiseVar);
                }
                fpRef.setPixels(pixelsRef);
            }
            imsRef.setProcessor(fpRef, i+1);

            log.status("Trial "+(i+1)+" - running SQUIRREL");

            // PSO
            long SQUIRRELStart = System.currentTimeMillis();

            kPSOErrorMap.maxMagnification = magnification;
            double sigmaGuess=0;
            double[] results;
            results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
                    new double[]{0.001, -1000, 0.5}, // low boundary
                    new double[]{100, 1000, 10}, // high boundary
                    new double[]{10, 0, sigmaGuess==0?5:sigmaGuess}, // best guess
                    new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, sigmaGuess==0?OBEY_BOUNDARY:CONSTANT}, // impose boundaries
                    1e-7);

            long SQUIRRELStop = System.currentTimeMillis();

            SQUIRRELTimes[i] = SQUIRRELStop-SQUIRRELStart;
            SQUIRRELAlpha[i] = results[0];
            SQUIRRELBeta[i] = results[1];
            SQUIRRELSigma[i] = results[2]*magnification;

            log.status("Trial "+(i+1)+" - running UNIQORN!");

            // run UNIQORN
            long UNIQORNStart = System.currentTimeMillis();

            float[] pixelsRef = (float[]) fpRef.duplicate().convertToFloatProcessor().getPixels();
            int nPixelsSR = fpSR.getPixelCount();

            float[] ones = new float[nPixelsSR];
            for(int n=0; n<nPixelsSR; n++){ones[n] = 1;}

            // UNIVARIATE OPTIMIZER - LINEAR MATCHING
            /// setup optimizer
            UnivariateFunction f =  new SQUIRRELVsUNIQORNTest_.sigmaOptimiseFunction(fpSR.duplicate().convertToFloatProcessor(), pixelsRef, ones);
            UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
            /// run optimizer
            UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                    new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, 5*magnification)); //NYQUIST NOT ASSUMED; limit to 5*magnification
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

            UNIQORNTimes[i] = UNIQORNStop-UNIQORNStart;
            UNIQORNAlpha[i] = alpha;
            UNIQORNBeta[i] = beta;
            UNIQORNSigma[i] = sigma_linear;

        }

        // display images
        new ImagePlus("GT images", imsSR).show();
        new ImagePlus("Reference images", imsRef).show();

        // populate results table
        ResultsTable rt = new ResultsTable();

        for (int i = 0; i < nRepeats; i++) {
            rt.incrementCounter();
            rt.addValue("Trial number", i);
            rt.addValue("Structure type", structure[i]);
            rt.addValue("% coverage", coverage[i] * 100);
            if(addNoise){
                rt.addValue("Noise variance", noiseVars[i]);
            }
            rt.addValue("True alpha", trueAlpha[i]);
            rt.addValue("True beta", trueBeta[i]);
            rt.addValue("True sigma", trueSigma[i]);
            rt.addValue("SQUIRREL alpha", SQUIRRELAlpha[i]);
            rt.addValue("SQUIRREL beta", SQUIRRELBeta[i]);
            rt.addValue("SQUIRREL sigma", SQUIRRELSigma[i]);
            rt.addValue("UNIQORN! alpha", UNIQORNAlpha[i]);
            rt.addValue("UNIQORN! beta", UNIQORNBeta[i]);
            rt.addValue("UNIQORN! sigma", UNIQORNSigma[i]);
            rt.addValue("SQUIRREL runtime", SQUIRRELTimes[i]);
            rt.addValue("UNIQORN! runtime", UNIQORNTimes[i]);
        }

        rt.show("Results");
    }

    FloatProcessor getDotImage(double targetCoverage){
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

    FloatProcessor getLineImage(double targetCoverage){
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

    FloatProcessor getPatchImage(double targetCoverage){
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
