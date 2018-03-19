package nanoj.squirrel.java.gui;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.process.FloatProcessor;
import nanoj.squirrel.java.ISuckAtCL;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import java.io.IOException;

public class AddTestCL_ extends _BaseSQUIRRELDialog_ {

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = false;
        useSettingsObserver = true;
        return true;
    }

    @Override
    public void setupDialog(){}

    @Override
    public boolean loadSettings(){
        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {
        imp = WindowManager.getCurrentImage();
        ImageStack ims = imp.getImageStack();

        FloatProcessor fp1, fp2;

        fp1 = ims.getProcessor(1).convertToFloatProcessor();
        fp2 = ims.getProcessor(2).convertToFloatProcessor();

        ISuckAtCL iSuckAtCL = new ISuckAtCL(imp.getWidth(), imp.getHeight());

        FloatProcessor fpSummed = iSuckAtCL.justAdd(fp1, fp2);

        new ImagePlus("ta-da, I guess", fpSummed).show();

    }

}
