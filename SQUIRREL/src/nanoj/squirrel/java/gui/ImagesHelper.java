package nanoj.squirrel.java.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import nanoj.core.java.tools.Log;
import nanoj.core.java.tools.Prefs;

import static nanoj.core.java.array.ArrayMath.normalizeIntegratedIntensity;
import static nanoj.core.java.imagej.FilesAndFoldersTools.getOpenPath;

/**
 * Created by paxcalpt on 02/06/2017.
 */
public class ImagesHelper {

    private String prefsHeader = this.getClass().getName();
    private Prefs prefs = new Prefs();
    private Log log = new Log();

    private String[] imageTitles;
    public String titleSRImage = "";
    public String titleRefImage = "";
    public String titleRSFImage = "";
    public ImagePlus impRef, impSR, impRSF;
    public int nSlicesRef, nSlicesSR, nSlicesRSF;
    public ImageStack imsRef;
    public ImageStack imsSR;
    public ImageStack imsRSF;
    public boolean check_sameNumberFrames_RefAndSR = true;
    public boolean isRSFMagnified = false;

    public boolean preset_SR_CSR_RSF() {
        int nImages = WindowManager.getImageCount();

        if (nImages < 2) {
            String path;
            ImagePlus imp;

            log.status("Requesting reference, SR and RSF images...");

            // get Reference image
            path = getOpenPath("Select Reference image", null);
            if (path == null) return false;
            imp = IJ.openImage(path);
            imp.show();
            titleRefImage = imp.getTitle();

            // get SR image
            path = getOpenPath("Select Super-Resolved image", null);
            if (path == null) return false;
            imp = IJ.openImage(path);
            imp.show();
            titleSRImage = imp.getTitle();

            // get RSF image
            path = getOpenPath("Select RSF image", null);
            if (path == null) return false;
            imp = IJ.openImage(path);
            imp.show();
            titleRSFImage = imp.getTitle();

            log.status("");
        }

        return true;
    }

    public boolean set_SR_CSR_RSF_Dialog(NonBlockingGenericDialog gd) {

        if (titleRefImage.equals("")) {
            titleRefImage = prefs.get(prefsHeader + ".titleRefImage", "");
            titleSRImage = prefs.get(prefsHeader + ".titleSRImage", "");
            titleRSFImage = prefs.get(prefsHeader + ".titleRSFImage", "");
        }

        imageTitles = WindowManager.getImageTitles();

        if (!titleRefImage.equals("") && contains(titleRefImage, imageTitles))
            gd.addChoice("Reference image", imageTitles, titleRefImage);
        else
            gd.addChoice("Reference image", imageTitles, imageTitles[0]);

        if (!titleSRImage.equals("") && contains(titleSRImage, imageTitles))
            gd.addChoice("Super-Resolution reconstructions", imageTitles, titleSRImage);
        else
            gd.addChoice("Super-Resolution reconstructions", imageTitles, imageTitles[0]);

        if (!titleRSFImage.equals("") && contains(titleRSFImage, imageTitles))
            gd.addChoice("RSF image", imageTitles, titleRSFImage);
        else
            gd.addChoice("RSF image", imageTitles, imageTitles[0]);

        gd.addChoice("RSF size is scaled to", new String[]{"Reference", "Super-Resolution Reconstruction"}, prefs.get(prefsHeader + ".RSFScaledTo", "Reference"));

        return true;
    }

    public boolean get_SR_CSR_RSF_Dialog(NonBlockingGenericDialog gd) {
        titleRefImage = gd.getNextChoice();
        titleSRImage = gd.getNextChoice();
        titleRSFImage = gd.getNextChoice();
        String RSFScaledTo = gd.getNextChoice();

        if (RSFScaledTo.equals("Reference")) isRSFMagnified = false;
        else isRSFMagnified = true;

        prefs.set(prefsHeader + ".titleRefImage", titleRefImage);
        prefs.set(prefsHeader + ".titleSRImage", titleSRImage);
        prefs.set(prefsHeader + ".titleRSFImage", titleRSFImage);
        prefs.set(prefsHeader + ".RSFScaledTo", RSFScaledTo);
        prefs.savePreferences();

        impRef = WindowManager.getImage(titleRefImage);
        impSR = WindowManager.getImage(titleSRImage);
        impRSF = WindowManager.getImage(titleRSFImage);

        imsRef = impRef.getImageStack().convertToFloat();
        imsSR = impSR.getImageStack().convertToFloat();
        imsRSF = impRSF.getImageStack().convertToFloat();

        nSlicesSR = impSR.getStackSize();
        nSlicesRef = impRef.getStackSize();
        nSlicesRSF = impRSF.getStackSize();

        // check number of frames in Ref
        if (check_sameNumberFrames_RefAndSR && nSlicesRef != 1 && nSlicesSR != nSlicesRef) {
            log.error("Number of reference frames need to be 1 or the same as number of SR frames");
            return false;
        }

        // check number of frames in RSF
        if (imsRSF.getSize() != 1 && nSlicesRSF != nSlicesSR){
            IJ.error("Number of slices in RSF image(s) needs to be 1 or the same number of SR reconstructions.");
            return false;
        }

        // expand RSF and normalize
        if(imsRSF.getSize() == 1) {
            imsRSF = new ImageStack(impRSF.getWidth(), impRSF.getHeight(), nSlicesSR);
            for(int r=1; r<=nSlicesSR; r++){
                imsRSF.setProcessor(impRSF.getProcessor().convertToFloatProcessor(), r);
            }
        }
        for(int r=1; r<=nSlicesSR; r++){
            normalizeIntegratedIntensity((float[]) imsRSF.getProcessor(r).getPixels(), 1);
        }

        return true;
    }

    public ImageStack magnifyRef() {
        int w = imsSR.getWidth();
        int h = imsSR.getHeight();

        ImageStack imsRefResized = new ImageStack(w, h, nSlicesRef);
        for (int s = 1; s <= nSlicesRef; s++) {
            FloatProcessor fpRefResized = (FloatProcessor) imsRef.getProcessor(s);
            fpRefResized.setInterpolationMethod(fpRefResized.BICUBIC);
            fpRefResized = (FloatProcessor) fpRefResized.resize(w, h);
            imsRefResized.setProcessor(fpRefResized, s);
        }
        return imsRefResized;
    }

    public ImageStack magnifyRSF() {
        int w = imsRSF.getWidth();
        int h = imsRSF.getHeight();
        int magnification = imsSR.getWidth()/imsRef.getWidth();
        int newW = w*magnification;
        int newH = w*magnification;
        if (newW % 2 == 0) newW++;
        if (newH % 2 == 0) newH++;

        ImageStack imsRSFResized = new ImageStack(newW, newH, nSlicesRSF);
        for (int s = 1; s <= nSlicesRef; s++) {
            FloatProcessor fpRSFResized = (FloatProcessor) imsRSF.getProcessor(s);
            fpRSFResized.setInterpolationMethod(fpRSFResized.BICUBIC);
            fpRSFResized = (FloatProcessor) fpRSFResized.resize(newW, newH);
            normalizeIntegratedIntensity((float[]) fpRSFResized.getPixels(), 1);
            imsRSFResized.setProcessor(fpRSFResized, s);
        }
        return imsRSFResized;
    }

    private boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
    }


    public static ImageStack magnify(ImageStack ims, int w, int h) {

        ImageStack imsResized = new ImageStack(w, h, ims.getSize());
        for (int s = 1; s <= ims.getSize(); s++) {
            FloatProcessor fpRefResized = ims.getProcessor(s).convertToFloatProcessor();
            fpRefResized.setInterpolationMethod(fpRefResized.BICUBIC);
            fpRefResized = (FloatProcessor) fpRefResized.resize(w, h);
            imsResized.setProcessor(fpRefResized, s);
        }
        return imsResized;
    }
}
