package nanoj.squirrel.java.gui.dataPreparation;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import java.io.IOException;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static nanoj.core.java.image.registration.CrossCorrelationElastic.applyElasticTransform;
import static nanoj.core.java.image.registration.CrossCorrelationElastic.calculateTranslationMask;
import static nanoj.squirrel.java.gui.ImagesHelper.magnify;

/**
 * Created with IntelliJ IDEA.
 * User: Ricardo Henriques <paxcalpt@gmail.com>
 * Date: 16/10/15
 * Time: 16:51
 */
public class RegisterReconstructionsElastic_ extends _BaseSQUIRRELDialog_ {

    private int blocksPerAxis, w, h, maxExpectedMisalignment;
    private double minSimilarity;

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

        gd = new NonBlockingGenericDialog("Super-Resolution Image Registration - Elastic");
        gd.addChoice("Reference image", imageTitles, imageTitles[0]);
        gd.addChoice("Super-Resolution reconstructions", imageTitles, imageTitles[1]);
        gd.addNumericField("Maximum expected misalignment (0-auto)", getPrefs("maxExpectedMisalignment", 0), 0);
        gd.addNumericField("Blocks per axis (default: 5)", getPrefs("blocksPerAxis", 5), 0);
        gd.addNumericField("Min similarity (default: 0.5, range 0-1)", getPrefs("minSimilarity", 0.5), 2);
    }

    @Override
    public boolean loadSettings() {
        String titleRefImage = gd.getNextChoice();
        String titleSRImage = gd.getNextChoice();
        maxExpectedMisalignment = (int) max(gd.getNextNumber(), 0);
        blocksPerAxis = (int) max(gd.getNextNumber(), 1);
        minSimilarity = min(max(gd.getNextNumber(), 0), 1);


        impRef = WindowManager.getImage(titleRefImage);
        impSR = WindowManager.getImage(titleSRImage);

        setPrefs("maxExpectedMisalignment", maxExpectedMisalignment);
        setPrefs("blocksPerAxis", blocksPerAxis);
        setPrefs("minSimilarity", minSimilarity);
        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {
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
        ImageStack imsTranslationMask = null;
        ImageStack imsBlockCrossCorrelation = null;

        for(int s=1; s<=nSlicesSR; s++) {
            log.progress(s, nSlicesSR);

            if (imsRef.getSize() > 1) fpRef = imsRefResized.getProcessor(s).convertToFloatProcessor();

            FloatProcessor fpSR = imsSR.getProcessor(s).convertToFloatProcessor();

            log.status("Calculating registration mask...");
            FloatProcessor[] fps = calculateTranslationMask(fpSR, fpRef, blocksPerAxis, maxExpectedMisalignment, minSimilarity);
            FloatProcessor fpTranslationMask = fps[0];
            FloatProcessor fpBCC = fps[1];
            FloatProcessor fpSRTranslated = applyElasticTransform(fpSR, fpTranslationMask);

            if (imsTranslationMask == null) imsTranslationMask = new ImageStack(fpTranslationMask.getWidth(), fpTranslationMask.getHeight(), nSlicesSR);
            if (imsBlockCrossCorrelation == null) imsBlockCrossCorrelation = new ImageStack(fpBCC.getWidth(), fpBCC.getHeight(), nSlicesSR);

            imsSRTranslated.setProcessor(fpSRTranslated, s);
            imsTranslationMask.setProcessor(fpTranslationMask, s);
            imsBlockCrossCorrelation.setProcessor(fpBCC, s);
        }

        ImagePlus impSRTranslated = impSR.duplicate(); // to copy the settings of the original image
        impSRTranslated.setTitle(impSR.getTitle().replace(".tif", "")+" - Aligned");
        impSRTranslated.setStack(imsSRTranslated);
        impSRTranslated.show();
        IJ.run(impSRTranslated, "Enhance Contrast...", "saturated=0.3");

        new ImagePlus("Translation Mask", imsTranslationMask).show();
        new ImagePlus("Block-Wise Cross-Correlation", imsBlockCrossCorrelation).show();
    }
}
