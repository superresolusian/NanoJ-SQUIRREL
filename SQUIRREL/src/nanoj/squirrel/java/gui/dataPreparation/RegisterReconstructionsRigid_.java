package nanoj.squirrel.java.gui.dataPreparation;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import static nanoj.core.java.image.drift.EstimateShiftAndTilt.MAX_FITTING;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.getShiftFromCrossCorrelationPeak;
import static nanoj.core.java.image.transform.CrossCorrelationMap.calculateCrossCorrelationMap;
import static nanoj.core.java.image.transform.CrossCorrelationMap.cropCCM;
import static nanoj.squirrel.java.gui.ImagesHelper.magnify;


/**
 * Created by sculley on 08/09/2016.
 */
public class RegisterReconstructionsRigid_ extends _BaseSQUIRRELDialog_ {

    private int maxExpectedMisalignment;
    private ImagePlus impRef;
    private ImagePlus impSR;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = false;
        useSettingsObserver = false;
        if (WindowManager.getImageTitles().length < 2) {
            log.error("Need at least 2 images open...");
            return false;
        }
        return true;
    }

    @Override
    public void setupDialog() {
        String[] imageTitles = WindowManager.getImageTitles();

        gd = new NonBlockingGenericDialog("Super-Resolution Image Registration - Rigid");
        gd.addChoice("Reference image", imageTitles, imageTitles[0]);
        gd.addChoice("Super-Resolution reconstructions", imageTitles, imageTitles[1]);
        gd.addNumericField("Maximum expected misalignment (0-auto)", getPrefs("maxExpectedMisalignment", 0), 0);
    }

    @Override
    public boolean loadSettings() {
        String titleRefImage = gd.getNextChoice();
        String titleSRImage = gd.getNextChoice();
        maxExpectedMisalignment = (int) gd.getNextNumber();

        setPrefs("maxExpectedMisalignment", maxExpectedMisalignment);

        impRef = WindowManager.getImage(titleRefImage);
        impSR = WindowManager.getImage(titleSRImage);

        savePrefs();
        return true;
    }

    @Override
    public void execute() {

        // Get images
        ImageStack imsRef = impRef.getImageStack();
        ImageStack imsSR = impSR.getImageStack();

        int nSlicesSR = imsSR.getSize();
        int nSlicesRef = imsRef.getSize();
        int width = imsSR.getWidth();
        int height = imsSR.getHeight();

        // Magnify ref image
        ImageStack imsRefResized = magnify(imsRef, width, height);

        // Cross-correlation realignment
        FloatProcessor fpRef = imsRefResized.getProcessor(1).convertToFloatProcessor();
        ImageStack imsSRTranslated = new ImageStack(width, height, nSlicesSR);
        ResultsTable rt = new ResultsTable();

        for(int s=1; s<=nSlicesSR; s++) {
            log.progress(s, nSlicesSR);

            if (imsRef.getSize() > 1) fpRef = imsRefResized.getProcessor(s).convertToFloatProcessor();

            FloatProcessor fpSR = imsSR.getProcessor(s).convertToFloatProcessor();

            log.status("calculating cross-correlation...");
            FloatProcessor ipCCM = (FloatProcessor) calculateCrossCorrelationMap(fpRef, fpSR, true);

            if (maxExpectedMisalignment != 0) ipCCM = cropCCM(ipCCM, maxExpectedMisalignment);

            log.status("calculating cross-correlation peaks...");

            float[] drift = getShiftFromCrossCorrelationPeak(ipCCM, MAX_FITTING);
            float driftX = drift[1];
            float driftY = drift[2];
            rt.incrementCounter();
            rt.addValue("Frame", s);
            rt.addValue("x-shift (pixels)", driftX);
            rt.addValue("y-shift (pixels)", driftY);

            log.status("translating SR image...");
            fpSR.setInterpolationMethod(fpSR.BICUBIC);
            fpSR.translate(driftX, driftY);
            imsSRTranslated.setProcessor(fpSR, s);
        }

        ImagePlus impSRTranslated = impSR.duplicate(); // to copy the settings of the original image
        impSRTranslated.setTitle(impSR.getTitle().replace(".tif", "")+" - Aligned");
        impSRTranslated.setStack(imsSRTranslated);
        impSRTranslated.show();
        IJ.run(impSRTranslated, "Enhance Contrast...", "saturated=0.3");
        rt.show("Translation");
    }
}
