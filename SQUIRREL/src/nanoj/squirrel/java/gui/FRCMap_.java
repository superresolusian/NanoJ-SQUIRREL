package nanoj.squirrel.java.gui;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import nanoj.core.java.image.analysis.FRC;
import nanoj.core.java.threading.NanoJThreadExecutor;
import nanoj.core.java.tools.NJ_LUT;
import nanoj.kernels.Kernel_VoronoiImage;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import java.io.IOException;
import java.util.ArrayList;

import static java.lang.Math.*;

/**
 * Created by paxcalpt on 14/03/2017.
 */
public class FRCMap_ extends _BaseSQUIRRELDialog_ {
    private int blocksPerXAxis, blocksPerYAxis, nSlices, w, h;
    private ImageStack ims;
    private final static Kernel_VoronoiImage VI = new Kernel_VoronoiImage();
    private double pixelSize;

    String[] thresholdMethods = new String[]{"1/7", "Half bit", "One bit", "Two bit", "One sigma", "Two sigma",
                                            "Three sigma", "Four sigma"};
    String thresholdMethod;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        useSettingsObserver = true;
        return true;
    }

    @Override
    public void setupDialog() {
        w = imp.getWidth();
        h = imp.getHeight();

        ims = imp.getImageStack();
        nSlices = ims.getSize();

        if (nSlices < 2) {
            log.error("Image-stack for FRC analysis needs at least 2 frames");
            return;
        }

        gd = new NonBlockingGenericDialog("Calculate FRC Map...");
        gd.addNumericField("Blocks per axis (default: 10)", getPrefs("blocksPerAxis", 10), 0);
        gd.addNumericField("Pixel size (nm)", getPrefs("pixelSize", 25), 2);
        gd.addChoice("Threshold method", thresholdMethods, getPrefs("thresholdMethod", thresholdMethods[0]));

        gd.addCheckbox("Show preview", false);

    }

    @Override
    public boolean loadSettings() {
        int blocksPerAxis = (int) max(gd.getNextNumber(), 1);
        pixelSize = gd.getNextNumber();

        blocksPerXAxis = blocksPerAxis;
        blocksPerYAxis = (int) (blocksPerXAxis * (((double) h)/w));
        thresholdMethod = gd.getNextChoice();

        showPreview = gd.getNextBoolean();

        setPrefs("blocksPerAxis", blocksPerAxis);
        setPrefs("pixelSize", pixelSize);
        setPrefs("thresholdMethod", thresholdMethod);
        return true;
    }

    @Override
    public void doPreview() {
        log.status("calculating preview...");

        FloatProcessor ip1 = ims.getProcessor(1).convertToFloatProcessor();
        FloatProcessor ip2 = ims.getProcessor(2).convertToFloatProcessor();
        FRCData results = calculateFRCMap(ip1, ip2);

        if (impPreview != null) impPreview.setProcessor(results.FRCMap);
        else impPreview = new ImagePlus("Preview...", results.FRCMap);
        impPreview.show();

        log.msg("FRC-Resolution: Mean="+round(results.meanFRC)+"nm Std-Dev="+round(results.stdDevFRC)+"nm Min="+round(results.minFRC)+"nm N-Blocks="+results.nValues);
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        ImageStack imsMap = new ImageStack(w, h);
        ImageStack imsMask = new ImageStack(w, h);
        ResultsTable rt = new ResultsTable();

        for (int s=1; s<nSlices; s++) {
            log.status("Calculating FRC: frame "+s+"-"+(s+1));
            log.progress(s, nSlices);
            FloatProcessor ip1 = ims.getProcessor(s).convertToFloatProcessor();
            FloatProcessor ip2 = ims.getProcessor(s+1).convertToFloatProcessor();
            FRCData results = calculateFRCMap(ip1, ip2);
            imsMap.addSlice(results.FRCMap);
            imsMask.addSlice(results.FRCMask);

            rt.incrementCounter();
            rt.addValue("N-Blocks", results.nValues);
            rt.addValue("Mean (nm)", results.meanFRC);
            rt.addValue("Std-Dev (nm)", results.stdDevFRC);
            rt.addValue("Min FRC (nm)", results.minFRC);
            rt.addValue("Max FRC (nm)", results.maxFRC);
        }

        ImagePlus impMap = new ImagePlus("FRC Map", imsMap);
        NJ_LUT.applyLUT_SQUIRREL_FRC(impMap);
        impMap.show();

        //new ImagePlus("FRC Mask", imsMask).show();
        rt.show("FRC-Resolution");

        log.msg("Threshold method is "+thresholdMethod);
    }

    public FRCData calculateFRCMap(FloatProcessor ip1, FloatProcessor ip2) {

        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
        FloatProcessor ipMask = new FloatProcessor(w, h);

        ArrayList<double[]> vectors = new ArrayList<double[]>();
        for (int nYB = 0; nYB < blocksPerYAxis; nYB++) {
            for (int nXB = 0; nXB < blocksPerXAxis; nXB++) {
                NTE.execute(new CalculateFRCBlock(nXB, nYB, ip1, ip2, ipMask, vectors));
            }
        }
        NTE.finish();

        int nValues = vectors.size();
        float[] xPositions = new float[nValues];
        float[] yPositions = new float[nValues];
        float[] values = new float[nValues];
        int counter = 0;
        double mean = 0, max = -Double.MAX_VALUE, min = Double.MAX_VALUE;
        for (double[] vector: vectors) {
            xPositions[counter] = (float) vector[0];
            yPositions[counter] = (float) vector[1];
            values[counter] = (float) vector[2];
            mean += vector[2];
            max = Math.max(vector[2], max);
            min = Math.min(vector[2], min);
            counter++;
        }
        mean /= nValues;
        // now get stdDev
        double stdDev = 0;
        for(int n=0;  n<nValues; n++){
            stdDev += Math.pow(values[n]-mean,2);
        }
        stdDev = (float) sqrt(stdDev/nValues);

        FRCData data = new FRCData();
        data.nValues = nValues;
        data.meanFRC = mean;
        data.stdDevFRC = stdDev;
        data.maxFRC = max;
        data.minFRC = min;
        data.FRCMask = ipMask;
        data.FRCMap = VI.calculateImage(w, h, xPositions, yPositions, values);
        return data;
    }

    class CalculateFRCBlock extends Thread {
        private final int nXB;
        private final int nYB;
        private final FRC myFRC = new FRC();
        private final FloatProcessor ip1;
        private final FloatProcessor ip2;
        private final FloatProcessor ipMask;
        private final ArrayList<double[]> vectors;

        public CalculateFRCBlock(int nXB, int nYB, FloatProcessor ip1, FloatProcessor ip2, FloatProcessor ipMask, ArrayList<double[]> vectors) {
            this.nXB = nXB;
            this.nYB = nYB;
            this.ip1 = ip1;
            this.ip2 = ip2;
            this.ipMask = ipMask;
            this.vectors = vectors;
        }

        private FloatProcessor getROI(FloatProcessor ip, int x, int y, int w, int h) {
            FloatProcessor ipCrop = new FloatProcessor(w,h);
            for (int j=0; j<h; j++) {
                for (int i=0; i<w; i++) {
                    ipCrop.setf(i, j, ip.getf(x+i, y+j));
                }
            }
            return ipCrop;
        }

        public void run() {
            int blockWidth = ip1.getWidth() / blocksPerXAxis;
            int blockHeight = ip1.getHeight() / blocksPerYAxis;
            int xStart = nXB * blockWidth;
            int yStart = nYB * blockHeight;
            blockWidth = Math.min(blockWidth, ip1.getWidth()-xStart);
            blockHeight = Math.min(blockHeight, ip1.getHeight()-yStart);

            FloatProcessor ipROI1 = getROI(ip1, xStart, yStart, blockWidth, blockHeight);
            FloatProcessor ipROI2 = getROI(ip2, xStart, yStart, blockWidth, blockHeight);
            double resolution = myFRC.calculateFireNumber(ipROI1, ipROI2, setMethod(thresholdMethod));

            if (Double.isNaN(resolution)) return;

            for (int j=0; j<blockHeight; j++) {
                for (int i=0; i<blockWidth; i++) {
                    ipMask.setf(xStart+i, yStart+j, 1);
                }
            }
            vectors.add(new double[]{xStart+blockWidth/2., yStart+blockHeight/2., resolution*pixelSize});
        }
    }

    public FRC.ThresholdMethod setMethod(String method){
        if(method==thresholdMethods[0]) return FRC.ThresholdMethod.FIXED_1_OVER_7;
        if(method==thresholdMethods[1]) return FRC.ThresholdMethod.HALF_BIT;
        if(method==thresholdMethods[2]) return FRC.ThresholdMethod.ONE_BIT;
        if(method==thresholdMethods[3]) return FRC.ThresholdMethod.TWO_BIT;
        if(method==thresholdMethods[4]) return FRC.ThresholdMethod.ONE_SIGMA;
        if(method==thresholdMethods[5]) return FRC.ThresholdMethod.TWO_SIGMA;
        if(method==thresholdMethods[6]) return FRC.ThresholdMethod.THREE_SIGMA;
        if(method==thresholdMethods[7]) return FRC.ThresholdMethod.FOUR_SIGMA;
        return FRC.ThresholdMethod.FIXED_1_OVER_7;
    }

}

class FRCData {
    public FloatProcessor FRCMap;
    public FloatProcessor FRCMask;
    public double meanFRC;
    public double stdDevFRC;
    public double minFRC;
    public double maxFRC;
    public int nValues;
}
