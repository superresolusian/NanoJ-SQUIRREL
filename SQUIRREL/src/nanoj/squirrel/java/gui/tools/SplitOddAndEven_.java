package nanoj.squirrel.java.gui.tools;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import java.io.IOException;

/**
 * Created by Henriques-lab on 23/03/2017.
 */
public class SplitOddAndEven_ extends _BaseSQUIRRELDialog_ {
    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        useSettingsObserver = false;
        return true;
    }

    @Override
    public void setupDialog() {

    }

    @Override
    public boolean loadSettings() {
        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {
        ImageStack ims = imp.getImageStack();
        if (ims.getSize() < 2) {
            log.error("Image must be a stack!!");
            return;
        }

        int w = ims.getWidth();
        int h = ims.getHeight();

        ImageStack imsOdd = new ImageStack(w, h);
        ImageStack imsEven = new ImageStack(w, h);

        for (int n=1; n<=ims.getSize(); n++) {
            ImageProcessor ip = ims.getProcessor(n);
            if (n%2!=0) imsOdd.addSlice(ip);
            else imsEven.addSlice(ip);
        }

        new ImagePlus(imp.getTitle()+" - Odd Frames", imsOdd).show();
        new ImagePlus(imp.getTitle()+" - Even Frames", imsEven).show();
    }
}
