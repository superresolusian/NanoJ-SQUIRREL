package openCLfree.squirrel.java.tools;

import ij.IJ;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;

import static java.lang.Math.*;

/**
 * Created by sculley on 11/07/2019.
 */
public class SQUIRREL_CrossCorrelationMap {
    
    public static boolean showProgress = false;

    public static ImageProcessor calculateCrossCorrelationMap(ImageProcessor ip1, ImageProcessor ip2, boolean normalized) {
        int w1 = ip1.getWidth();
        int h1 = ip1.getHeight();
        int w2 = ip2.getWidth();
        int h2 = ip2.getHeight();
        if (w1!=w2 && h1!=h2) {
            IJ.error("Both comparison images don't have same size");
            return null;
        }
        if (!SQUIRREL_FHT.isEvenSquare(ip1)) {
            int size = SQUIRREL_FHT.getClosestEvenSquareSize(ip1);
            ip1.setRoi((w1-size)/2, (h1-size)/2, size, size);
            ip2.setRoi((w1-size)/2, (h1-size)/2, size, size);
            ip1 = ip1.crop();
            ip2 = ip2.crop();
        }

        return _calculateCrossCorrelationMap(ip1.convertToFloatProcessor(), ip2.convertToFloatProcessor(), normalized);
    }

    /**
     * Assumes ip1 and ip2 are already even square
     * @param ip1
     * @param ip2
     */
    private static FloatProcessor _calculateCrossCorrelationMap(FloatProcessor ip1, FloatProcessor ip2, boolean normalized) {
        FloatProcessor h1 = SQUIRREL_FHT.forwardFHT(ip1);
        FloatProcessor h2 = SQUIRREL_FHT.forwardFHT(ip2);
        FloatProcessor ipCCM = SQUIRREL_FHT.conjugateMultiply(h1, h2);
        ipCCM = SQUIRREL_FHT.inverseFHT(ipCCM, false);
        SQUIRREL_FHT.swapQuadrants(ipCCM);
        SQUIRREL_FHT.flip(ipCCM);

        if (normalized) ipCCM = _normalizeCrossCorrelationMap(ip1, ip2, ipCCM);
        ipCCM.setRoi(new Rectangle(0, 0, ipCCM.getWidth()-1, ipCCM.getHeight()-1));
        ipCCM = (FloatProcessor) ipCCM.crop();

        return ipCCM;
    }

    private static FloatProcessor _normalizeCrossCorrelationMap(FloatProcessor ip1, FloatProcessor ip2, FloatProcessor ipCCM) {
        float[] ccmPixels = (float[]) ipCCM.getPixels();

        int w = ipCCM.getWidth();
        int h = ipCCM.getHeight();
        float vMax = -Float.MAX_VALUE;
        float vMin = Float.MAX_VALUE;
        int pMax = 0;
        int pMin = 0;

        for (int n=0; n<ccmPixels.length; n++) {
            float v = ccmPixels[n];
            if (v > vMax) {
                vMax = v;
                pMax = n;
            }
            if (v < vMin) {
                vMin = v;
                pMin = n;
            }
        }

        int shiftXMax = (pMax % w) - w/2;
        int shiftYMax = (pMax / w) - h/2;
        int shiftXMin = (pMin % w) - w/2;
        int shiftYMin = (pMin / w) - h/2;

        float maxPPMCC = calculatePPMCC(ip1, ip2, shiftXMax, shiftYMax);
        float minPPMCC = calculatePPMCC(ip1, ip2, shiftXMin, shiftYMin);

        // calculate max and min Pearson product-moment correlation coefficient
        float deltaV = vMax - vMin;
        float deltaP = maxPPMCC - minPPMCC;
        for (int n=0; n<ccmPixels.length; n++) {
            float v = (ccmPixels[n] - vMin) / deltaV;
            v = (v * (maxPPMCC - minPPMCC)) + minPPMCC;
            ccmPixels[n] = v;
        }

        return new FloatProcessor(w, h, ccmPixels);
    }

    public static float calculatePPMCC(FloatProcessor ip1, FloatProcessor ip2, int shiftX, int shiftY) {
        int w = ip1.getWidth();
        int h = ip1.getHeight();

        int newW = w - abs(shiftX);
        int newH = h - abs(shiftY);

        // shift ips and crop as needed
        int x0 = max(0, -shiftX);
        int y0 = max(0, -shiftY);
        int x1 = x0 + shiftX;
        int y1 = y0 + shiftY;
        ip1.setRoi(x0, y0, newW, newH);
        ip2.setRoi(x1, y1, newW, newH);
        ip1 = (FloatProcessor) ip1.crop();
        ip2 = (FloatProcessor) ip2.crop();
        float[] pixels1 = (float[]) ip1.getPixels();
        float[] pixels2 = (float[]) ip2.getPixels();

        // calculate means
        double mean1 = 0;
        double mean2 = 0;
        for (int n=0; n<pixels1.length; n++) {
            mean1 += (pixels1[n] - mean1) / (n+1);
            mean2 += (pixels2[n] - mean2) / (n+1);
        }

        // calculate correlation
        double covariance = 0;
        double squareSum1 = 0;
        double squareSum2 = 0;
        for (int n=0; n<pixels1.length; n++) {
            double v1 = pixels1[n] - mean1;
            double v2 = pixels2[n] - mean2;

            covariance += v1*v2;
            squareSum1 += v1*v1;
            squareSum2 += v2*v2;
        }

        if (squareSum1 == 0 || squareSum2 == 0) return 0;
        double PPMCC = covariance / sqrt(squareSum1 * squareSum2);
        return (float) PPMCC;
    }

    public static FloatProcessor cropCCM(FloatProcessor ipCCM, int radius) {
        if (radius != 0 && radius*2+1 < ipCCM.getWidth() && radius*2+1 < ipCCM.getHeight()) {
            int xStart = ipCCM.getWidth()/2 - radius;
            int yStart = ipCCM.getHeight()/2 - radius;
            ipCCM.setRoi(xStart, yStart, radius*2+1, radius*2+1);
            return (FloatProcessor) ipCCM.crop();
        }
        return ipCCM;
    }
}
