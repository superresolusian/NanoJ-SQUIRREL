package nanoj.squirrel.java.gui;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import java.awt.*;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static nanoj.squirrel.java.gui.ErrorMap_.contains;

/**
 * Created by sculley on 26/07/2016.
 */
public class FuseSRImages_ extends _BaseSQUIRRELDialog_ {

    private String titleSRImage = "", titleErrorImage = "";
    private ImagePlus impSR, impEM;
    private ImageStack imsSR, imsEM;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = false;
        useSettingsObserver = false;
        int nImages = WindowManager.getImageCount();
        if(nImages<2){
            log.error("At least 2 images are required to run this code!");
            return false;
        }
        return true;
    }

    @Override
    public void setupDialog() {
        String[] imageTitles = WindowManager.getImageTitles();

        gd = new NonBlockingGenericDialog("SQUIRREL image fusion");
        gd.addHelp("https://bitbucket.org/rhenriqueslab/nanoj-squirrel");

        if (titleSRImage.equals("")) {
            titleSRImage = getPrefs("titleSRImage", "");
            titleErrorImage = getPrefs("titleErrorImage", "");
        }

        if (!titleSRImage.equals("") && contains(titleSRImage, imageTitles))
            gd.addChoice("Super-Resolution Reconstructions", imageTitles, titleSRImage);
        else
            gd.addChoice("Super-Resolution Reconstructions", imageTitles, imageTitles[0]);

        gd.addMessage("The input super-resolution reconstructions must be aligned and intensity-matched", new Font("Arial", Font.ITALIC, 12));

        if (!titleErrorImage.equals("") && contains(titleErrorImage, imageTitles))
            gd.addChoice("Error Maps", imageTitles, titleErrorImage);
        else
            gd.addChoice("Error Maps", imageTitles, imageTitles[0]);

    }

    @Override
    public boolean loadSettings() {
        titleSRImage = gd.getNextChoice();
        titleErrorImage = gd.getNextChoice();

        impSR = WindowManager.getImage(titleSRImage);
        impEM = WindowManager.getImage(titleErrorImage);

        if (impSR.getNSlices() < 2) {
            log.warning("Need at least 2 Super-Resolution frames to fuse...");
            return false;
        }

        if (impSR.getNSlices() != impEM.getNSlices()) {
            log.warning("Number of Super-Resolution frames needs to match number of Error-Maps...");
            return false;
        }

        imsSR = impSR.getImageStack();
        imsEM = impEM.getImageStack();

        setPrefs("titleSRImage", titleSRImage);
        setPrefs("titleErrorImage", titleErrorImage);
        savePrefs();

        return true;
    }

    @Override
    public void execute() {

        // Get images
        int width = imsSR.getWidth();
        int height = imsSR.getHeight();
        int widthM1 = width-1;
        int heightM1 = height-1;
        int nSlices = imsSR.getSize();
        int nPixels = width * height;

        double[][] pixelsR = new double[nSlices][nPixels];
        double[] pixelsRMax = new double[nPixels];
        double[] pixelsWSum = new double[nPixels];
        double[] pixelsFusion = new double[nPixels];

        // calculate equation S7 & S8
        log.status("Smoothing error maps...");
        for (int s=0; s<nSlices; s++) {
            FloatProcessor ipEM = imsEM.getProcessor(s+1).convertToFloatProcessor();

            for (int y=0; y<height; y++) {
                for (int x=0; x<width; x++) {

                    double v = 0;
                    // 3 x 3 smoothing
                    for (int dy=-1; dy<=1; dy++) {
                        for (int dx=-1; dx<=1; dx++) {
                            v += ipEM.getf(
                                    min(max(x+dx, 0), widthM1),
                                    min(max(y+dy, 0), heightM1));
                        }
                    }
                    v /= 9;

                    int p = y * width + x;

                    pixelsR[s][p] = v;
                    pixelsRMax[p] = max(pixelsRMax[p], v);
                }
            }
        }

        // calculate equation S9 & S10
        log.status("Calculating weight matrix...");
        for (int s=0; s<nSlices; s++) {
            FloatProcessor ipSR = imsSR.getProcessor(s+1).convertToFloatProcessor();

            for (int p = 0; p < nPixels; p++) {
                double w = pixelsRMax[p] / max(pixelsR[s][p], 1);
                pixelsWSum[p] += w;
                pixelsFusion[p] += ipSR.getf(p) * w;
            }
        }
        log.status("Fusing data...");
        for (int p = 0; p < nPixels; p++) {
            pixelsFusion[p] /= pixelsWSum[p];
        }

        FloatProcessor ipFusion = new FloatProcessor(width, height, pixelsFusion);

        new ImagePlus(impSR.getTitle()+" - Fused", ipFusion).show();
    }
}