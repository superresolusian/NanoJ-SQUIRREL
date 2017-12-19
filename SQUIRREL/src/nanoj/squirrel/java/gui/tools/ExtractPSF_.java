package nanoj.squirrel.java.gui.tools;

import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.array.ArrayCasting;
import nanoj.core.java.array.ArrayMath;
import nanoj.core.java.featureExtraction.ExtractRois;
import nanoj.core.java.featureExtraction.Peaks;
import nanoj.core.java.threading.NanoJThreadExecutor;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;
import nanoj.squirrel.java.minimizers.GaussianFitMinimizer;

import java.awt.*;
import java.io.IOException;

import static nanoj.core.java.featureExtraction.Peaks.getROIs;
import static nanoj.core.java.featureExtraction.Peaks.populateRoiManagerWithPeaks;

/**
 * Created with IntelliJ IDEA.
 * User: Ricardo Henriques <paxcalpt@gmail.com>
 * Date: 02/06/17
 * Time: 14:06
 */
public class ExtractPSF_ extends _BaseSQUIRRELDialog_ {

    private RoiManager rm = null;
    private int radius, nROIs;
    private boolean showPeaks, crop;
    private float[][] peaks;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        useSettingsObserver = true;
        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Extract PSF...");
        gd.addMessage(
                "Note: this plugin is intended to run over a image containing discrete fluorophore detections.\n");

        gd.addNumericField("Max FWHM (pixels)", getPrefs("radius", 5), 0);
        gd.addNumericField("Number of ROIs to use", getPrefs("nROIs", 100), 0);

        gd.addCheckbox("Crop PSF size to meaningful data", getPrefs("crop", true));
        gd.addCheckbox("Show detected peaks", getPrefs("showPeaks", false));
    }

    @Override
    public boolean loadSettings() {
        // Grab data from dialog
        radius = (int) gd.getNextNumber();
        nROIs = (int) gd.getNextNumber();

        crop = gd.getNextBoolean();
        showPeaks = gd.getNextBoolean();

        setPrefs("radius", radius);
        setPrefs("nROIs", nROIs);
        setPrefs("crop", crop);

        setPrefs("showPeaks", showPeaks);

        prefs.savePreferences();

        showPreview = true;

        return true;
    }

    public void doPreview() {

        log.status("calculating preview...");

        FloatProcessor ip = imp.getProcessor().convertToFloatProcessor();
        peaks = Peaks.getPeaks(ip, nROIs, radius, 0.25f);
        rm = ExtractRois.getRoiManager();
        populateRoiManagerWithPeaks(peaks, radius, rm);
        rm.runCommand("Associate", "false");
    }

    @Override
    public void execute() throws InterruptedException, IOException {
        if (rm == null) doPreview();

        log.status("extracting ROIs...");

        ImageProcessor ip = imp.getProcessor();

        Roi[] rois = getROIs(peaks, radius+2);

        // Grabbing ROIs
        ImageStack imsRoisBigger = new ImageStack(2*(radius+2)+1,2*(radius+2)+1);

        for (Roi r: rois) {
            if (r == null) continue;
            Rectangle rect = r.getBounds();
            if (rect.x+rect.width >= ip.getWidth() || rect.y+rect.height >= ip.getHeight()) {
                continue;
            }
            ip.setRoi(r);
            imsRoisBigger.addSlice(ip.crop());
        }

        int wBigger = imsRoisBigger.getWidth();
        int hBigger = imsRoisBigger.getHeight();
        int nRois = imsRoisBigger.getSize();

        double[] sigmas = new double[nRois];
        double[] amplitudes = new double[nRois];
        ImageStack imsRoisBiggerAligned = new ImageStack(wBigger, hBigger, nRois);

        log.status("estimating PSFs...");
        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
        for (int n=0; n<nRois; n++) {
            ThreadedFitterAndRealigner t = new ThreadedFitterAndRealigner(imsRoisBigger, imsRoisBiggerAligned, n, sigmas, amplitudes);
            NTE.execute(t);
        }
        NTE.finish();

        int r2p1 = 2*radius+1;
        ImageStack imsRoisAligned = imsRoisBiggerAligned.crop(2,2,0,r2p1,r2p1,nRois);

        Plot plot = new Plot("Amplitude vs Sigma", "Sigma", "Amplitude");
        plot.add("cross", sigmas, amplitudes);
        plot.show();

        log.status("Calculating PSF estimate...");

        int nPixels = r2p1*r2p1;
        float[] pixelPSF = new float[nPixels];
        double amplitudesSum = ArrayMath.sum(amplitudes);

        for (int s=1; s<=nRois; s++) {
            float[] pixelsROI = (float[]) imsRoisAligned.getPixels(s);
            for (int p=0; p<nPixels;p++) pixelPSF[p] += pixelsROI[p] * amplitudes[s-1] / amplitudesSum;
            log.progress(s, nROIs);
        }
        ImageProcessor ipPSF = new FloatProcessor(r2p1, r2p1, pixelPSF);


        // iterate new model redoing weights based on correlation
        for (int i=0; i<3; i++) {
            log.status("Calculating PSF estimate... iteration "+(i+1)+"/3");

            double[] newPixelPSF = new double[nPixels];
            double weightSum = 0;

            for (int s=1; s<=nRois; s++) {
                float[] pixelsROI = (float[]) imsRoisAligned.getPixels(s);
                double weight = Math.pow(ArrayMath.calculatePPMCC(pixelsROI, pixelPSF, true), 4);
                imsRoisAligned.setSliceLabel("Weight = "+String.format("%.3g%n", weight), s);
                for (int p=0; p<nPixels;p++) newPixelPSF[p] += pixelsROI[p] * weight;
                weightSum += weight;
                log.progress(s, nROIs);
            }
            for (int p=0; p<nPixels;p++) newPixelPSF[p] /= weightSum;
            pixelPSF = ArrayCasting.doubleToFloat(newPixelPSF);
        }

        log.progress(1);
        log.status("Done...");

        if (showPeaks) {
            new ImagePlus("Detected Aligned Peaks", imsRoisAligned).show();
        }

        if (crop) {
            ImageProcessor ipPSFCropped = null;
            Rectangle roi = new Rectangle(0, 0, ipPSF.getWidth(), ipPSF.getHeight());

            float intensity = ArrayMath.getAbsAverageValue((float[]) ipPSF.getPixels()) * ipPSF.getPixelCount();
            float croppedIntensity = intensity;

            while (croppedIntensity > (intensity * 0.99f)) {
                roi.x += 1;
                roi.y += 1;
                roi.width -= 2;
                roi.height -= 2;
                ipPSF.setRoi(roi);
                ipPSFCropped = ipPSF.crop();
                croppedIntensity = ArrayMath.getAbsAverageValue((float[]) ipPSFCropped.getPixels()) * ipPSFCropped.getPixelCount();
            }

            roi.x -= 1;
            roi.y -= 1;
            roi.width += 2;
            roi.height += 2;
            ipPSF.setRoi(roi);
            ipPSF = ipPSF.crop();
        }
        ImagePlus impPSF = new ImagePlus("PSF - " + imp.getTitle(), ipPSF);
        impPSF.show();
    }

    class ThreadedFitterAndRealigner extends Thread {
        private final ImageStack imsOut;
        private final FloatProcessor ip;
        private final int n;
        private final double[] sigmas, amplitudes;
        private final int w, h;
        GaussianFitMinimizer fitter;

        ThreadedFitterAndRealigner(ImageStack imsIn, ImageStack imsOut, int n, double[] sigmas, double[] amplitudes) {
            this.imsOut = imsOut;
            this.n = n;
            this.sigmas = sigmas;
            this.amplitudes = amplitudes;
            this.ip = imsIn.getProcessor(n+1).convertToFloatProcessor();
            this.fitter = new GaussianFitMinimizer(ip, 1.5, ip.getWidth()/2d, ip.getHeight()/2d);
            this.w = ip.getWidth();
            this.h = ip.getHeight();
        }

        @Override
        public void run() {
            fitter.calculate();
            double xc = fitter.xc;
            double yc = fitter.yc;

            double deltaX = w/2d - xc - 0.5;
            double deltaY = h/2d - yc - 0.5;

            FloatProcessor ipAlignedNormalised = new FloatProcessor(w, h);

            ip.setInterpolationMethod(ip.BICUBIC);

            float vMax = -Float.MAX_VALUE;
            float vMin = Float.MAX_VALUE;

            for (int j=0; j<h; j++) {
                for (int i=0; i<w; i++) {
                    float v = (float) ip.getInterpolatedPixel(i-deltaX, j-deltaY); // not sure if the right direction?
                    vMax = Math.max(vMax, v);
                    vMin = Math.min(vMin, v);
                    ipAlignedNormalised.setf(i, j, v);
                }
            }

            // normalise
            float[] pixels = (float[]) ipAlignedNormalised.getPixels();
            for (int n=0; n<pixels.length; n++) pixels[n] = (pixels[n]-vMin)/(vMax-vMin);

            this.imsOut.setProcessor(ipAlignedNormalised, this.n+1);
            this.sigmas[this.n] = fitter.sigma;
            this.amplitudes[this.n] = fitter.amplitude;
        }
    }
}