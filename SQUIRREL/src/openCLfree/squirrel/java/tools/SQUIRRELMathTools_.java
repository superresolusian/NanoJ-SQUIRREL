package openCLfree.squirrel.java.tools;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 * Created by sculley on 11/07/2019.
 */
public class SQUIRRELMathTools_ {

    /**
     * returns the value of the mean of a float array
     * @param numbers
     * @return mean(numbers)
     */
    public static float getAverageValue(float[] numbers){
        double v = 0;
        for(int i=0; i<numbers.length; i++){
            v += numbers[i] / numbers.length;
        }
        return (float) v;
    }

    /**
     * returns the value of the mean of a double array
     * @param numbers
     * @return mean(numbers)
     */
    public static double getAverageValue(double[] numbers){
        double v = 0;
        for(int i=0; i<numbers.length; i++){
            v += numbers[i] / numbers.length;
        }
        return v;
    }

    public static double calculatePPMCC(float[] array1, float[] array2, boolean doMeanSubtraction) {
        if (doMeanSubtraction) {
            array1 = array1.clone();
            array2 = array2.clone();
            float mean;
            mean = getAverageValue(array1);
            for (int n=0; n<array1.length; n++) array1[n] -= mean;
            mean = getAverageValue(array2);
            for (int n=0; n<array2.length; n++) array2[n] -= mean;
        }

        double covariance = 0;
        double squareSum1  = 0;
        double squareSum2  = 0;
        for (int n=0; n<array1.length; n++) {
            float v0 = array1[n];
            float v1 = array2[n];
            covariance += v0*v1;
            squareSum1 += v0*v0;
            squareSum2 += v1*v1;
        }
        double similarity = 0;
        if (squareSum1 !=0 && squareSum2 != 0) similarity = covariance / sqrt(squareSum1 * squareSum2);
        return similarity;
    }

    /**
     * @param array1
     * @param array2
     * @return
     */
    public static double calculateMSE(float[] array1, float[] array2) {
        double MSE = 0;
        int counter = 1;

        for (int n=0; n<array1.length; n++) {
            float v0 = array1[n];
            float v1 = array2[n];
            if (Float.isNaN(v0) || Float.isNaN(v1)) continue;

            MSE += (pow(v0-v1,2)-MSE)/counter;
            counter++;
        }
        return MSE;
    }
}
