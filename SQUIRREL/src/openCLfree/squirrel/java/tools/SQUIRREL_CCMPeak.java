package openCLfree.squirrel.java.tools;

import ij.ImageStack;
import ij.measure.Minimizer;
import ij.measure.UserFunction;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/**
 * Created by sculley on 11/07/2019.
 */
public class SQUIRREL_CCMPeak {

    public static final int CENTER_OF_MASS = 0;
    public static final int MAXIMUM = 1;
    public static final int MAX_FITTING = 2;

    public SQUIRREL_CCMPeak() {
    }

    public static float[] getShiftFromCrossCorrelationPeak(FloatProcessor CCMap, int method) {
        int width = CCMap.getWidth();
        int height = CCMap.getHeight();
        float radiusX = (float)width / 2.0F;
        float radiusY = (float)height / 2.0F;
        float[] drift;
        if(method == 0) {
            drift = getCenterOfMass(CCMap);
        } else if(method == 2) {
            drift = getMaxFindByOptimization(CCMap);
        } else {
            drift = getMax(CCMap);
        }

        float shiftX = radiusX - drift[0] - 0.5F;
        float shiftY = radiusY - drift[1] - 0.5F;
        float shiftXY = (float)Math.sqrt(Math.pow((double)shiftX, 2.0D) + Math.pow((double)shiftY, 2.0D));
        float similarity = drift[2];
        return new float[]{shiftXY, shiftX, shiftY, similarity};
    }

    public static float[][] getShiftFromCrossCorrelationPeak(ImageStack imsCCMap, int method) {
        int nSlices = imsCCMap.getSize();
        float[] shiftXY = new float[nSlices];
        float[] shiftX = new float[nSlices];
        float[] shiftY = new float[nSlices];
        float[] similarity = new float[nSlices];

        for(int s = 0; s < nSlices; ++s) {
            FloatProcessor fp = imsCCMap.getProcessor(s + 1).convertToFloatProcessor();
            float[] shift = getShiftFromCrossCorrelationPeak(fp, method);
            shiftXY[s] = shift[0];
            shiftX[s] = shift[1];
            shiftY[s] = shift[2];
            similarity[s] = shift[3];
        }

        return new float[][]{shiftXY, shiftX, shiftY, similarity};
    }

    public static float[][] getShiftAndTiltFromRotationAndCorrelationPeak(ImageStack[] imsRCCMap, float angleStep) {
        int width = imsRCCMap[0].getWidth();
        int height = imsRCCMap[0].getHeight();
        int nSlices = imsRCCMap.length;
        double radiusX = (double)(width / 2);
        double radiusY = (double)(height / 2);
        float[] shiftXY = new float[nSlices];
        float[] shiftX = new float[nSlices];
        float[] shiftY = new float[nSlices];
        float[] theta = new float[nSlices];
        float[] similarity = new float[nSlices];

        for(int s = 0; s < nSlices; ++s) {
            ImageStack imsRCCMapSlice = imsRCCMap[s];
            float[] shiftAndTilt = getCenterOfMass(imsRCCMapSlice);
            theta[s] = shiftAndTilt[2] * angleStep;
            double x = (double)shiftAndTilt[0] - radiusX;
            double y = (double)shiftAndTilt[1] - radiusY;
            double ca = Math.cos((double)theta[s]);
            double sa = Math.sin((double)theta[s]);
            double xs = x * ca - y * sa;
            double ys = x * sa + y * ca;
            shiftX[s] = -((float)xs);
            shiftY[s] = -((float)ys);
            shiftXY[s] = (float)Math.sqrt(Math.pow((double)shiftX[s], 2.0D) + Math.pow((double)shiftY[s], 2.0D));
            similarity[s] = shiftAndTilt[3];
        }

        return new float[][]{shiftXY, shiftX, shiftY, theta, similarity};
    }

    public static float[] getMaxFindByOptimization(FloatProcessor fp) {
        MaxFindByOptimization MFO = new MaxFindByOptimization(fp);
        return MFO.calculate();
    }

    public static float[] getCenterOfMass(FloatProcessor fp) {
        float[] maxValues = getMax(fp);
        int xMax = (int)maxValues[0];
        int yMax = (int)maxValues[1];
        int xStart = Math.max(xMax - 2, 0);
        int yStart = Math.max(yMax - 2, 0);
        int xEnd = Math.min(xMax + 3, fp.getWidth());
        int yEnd = Math.min(yMax + 3, fp.getHeight());
        float xCM = 0.0F;
        float yCM = 0.0F;
        float sSum = 0.0F;
        float v = 0.0F;
        float vMax = 1.4E-45F;
        float vMin = 3.4028235E38F;

        int j;
        for(int vMaxMinusMin = yStart; vMaxMinusMin < yEnd; ++vMaxMinusMin) {
            for(j = xStart; j < xEnd; ++j) {
                v = fp.getf(j, vMaxMinusMin);
                vMax = Math.max(vMax, v);
                vMin = Math.min(vMin, v);
            }
        }

        float var17 = vMax - vMin;

        for(j = yStart; j < yEnd; ++j) {
            for(int i = xStart; i < xEnd; ++i) {
                v = (fp.getf(i, j) - vMin) / var17;
                xCM += (float)i * v;
                yCM += (float)j * v;
                sSum += v;
            }
        }

        xCM /= sSum;
        yCM /= sSum;
        return new float[]{xCM, yCM, vMax};
    }

    public static float[] getCenterOfMass(ImageStack ims) {
        if(ims.getProcessor(1).getBitDepth() != 32) {
            ims = ims.convertToFloat();
        }

        int width = ims.getWidth();
        int height = ims.getHeight();
        int depth = ims.getSize();
        float vMax = 0.0F;
        int xMax = 0;
        int yMax = 0;
        int zMax = 0;

        float v;
        int xStart;
        int zStart;
        for(xStart = 0; xStart < depth; ++xStart) {
            float[] yStart = (float[])((float[])ims.getPixels(xStart + 1));

            for(zStart = 0; zStart < yStart.length; ++zStart) {
                v = yStart[zStart];
                if(v > vMax) {
                    vMax = v;
                    xMax = zStart % width;
                    yMax = zStart / width;
                    zMax = xStart;
                }
            }
        }

        xStart = Math.max(xMax - 2, 0);
        int var27 = Math.max(yMax - 2, 0);
        zStart = zMax - 2;
        int xEnd = Math.min(xMax + 3, width);
        int yEnd = Math.min(yMax + 3, height);
        int zEnd = zMax + 3;
        double xCM = 0.0D;
        double yCM = 0.0D;
        double sSum = 0.0D;
        float minValue = vMax * 0.5F;
        FloatProcessor fp = (FloatProcessor)ims.getProcessor(zMax + 1);

        for(int zCM = var27; zCM < yEnd; ++zCM) {
            for(int i = xStart; i < xEnd; ++i) {
                v = fp.getf(i, zCM) - minValue;
                if(v >= 0.0F) {
                    xCM += (double)((float)i * v);
                    yCM += (double)((float)zCM * v);
                    sSum += (double)v;
                }
            }
        }

        xCM /= sSum;
        yCM /= sSum;
        double var28 = 0.0D;
        sSum = 0.0D;

        for(int z = zStart; z < zEnd; ++z) {
            int z_;
            if(z < 0) {
                z_ = ims.getSize() + z;
            } else if(z >= ims.getSize()) {
                z_ = z - ims.getSize();
            } else {
                z_ = z;
            }

            fp = (FloatProcessor)ims.getProcessor(z_ + 1);
            v = fp.getf(xMax, yMax) - minValue;
            if(v >= 0.0F) {
                var28 += (double)((float)z * v);
                sSum += (double)v;
            }
        }

        var28 /= sSum;
        return new float[]{(float)xCM, (float)yCM, (float)var28, vMax};
    }

    static public float[] getMax(ImageProcessor ip) {
        float vMax = -Float.MAX_VALUE;
        int pMax = 0;

        for (int p=0; p<ip.getPixelCount(); p++) {
            float v = ip.getf(p);
            if (v > vMax) {
                vMax = v;
                pMax = p;
            }
        }

        int w = ip.getWidth();
        int xMax = pMax % w;
        int yMax = pMax / w;

        return new float[] {xMax, yMax, vMax};
    }
}

class MaxFindByOptimization implements UserFunction {

    private final FloatProcessor fp;
    private int w, h;
    public float v, x, y;

    public MaxFindByOptimization(FloatProcessor fp) {
        this.fp = fp;
        this.w = fp.getWidth();
        this.h = fp.getHeight();
        fp.setInterpolationMethod(ImageProcessor.BICUBIC);
    }

    public float[] calculate() {
        float[] max = getMax(fp);
        float xMax = max[0];
        float yMax = max[1];

        double[] initialParameters = new double[] {xMax, yMax}; // sigma
        double[] initialParametersVariation = new double[] {1, 1};

        Minimizer min = new Minimizer();
        min.setFunction(this, 2);
        min.setMaxIterations(1000);
        min.setStatusAndEsc("Estimating max subpixel position: Iteration ", true);
        min.minimize(initialParameters, initialParametersVariation);
        x = (float) min.getParams()[0];
        y = (float) min.getParams()[1];
        v = (float) fp.getInterpolatedPixel((double) x, (double) y);

        return new float[]{x, y, v};
    }

    @Override
    public double userFunction(double[] params, double p) {
        double x = params[0];
        double y = params[1];
        if (x <= 1 || x >= w-2) return Double.NaN;
        if (y <= 1 || y >= h-2) return Double.NaN;
        return -fp.getInterpolatedPixel(x, y);
    }

    static public float[] getMax(ImageProcessor ip) {
        float vMax = -Float.MAX_VALUE;
        int pMax = 0;

        for (int p=0; p<ip.getPixelCount(); p++) {
            float v = ip.getf(p);
            if (v > vMax) {
                vMax = v;
                pMax = p;
            }
        }

        int w = ip.getWidth();
        int xMax = pMax % w;
        int yMax = pMax / w;

        return new float[] {xMax, yMax, vMax};
    }
}