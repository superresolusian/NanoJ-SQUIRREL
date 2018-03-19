package nanoj.squirrel.java.gui.tools;
import ij.gui.NonBlockingGenericDialog;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import java.awt.*;

public class SetMaximumStackSize_ extends _BaseSQUIRRELDialog_ {

    int maxSRImages;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = false;
        useSettingsObserver = true;
        return true;
    }


    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Set maximum expected stack size in SQUIRREL...");
        gd.addMessage("SQUIRREL checks stack sizes before running analysis." +
                " To speed up this process, you can set the maximum stack size for super-resolution" +
                " images to be checked in this dialog.", new Font("Arial", Font.ITALIC, 12));
        gd.addNumericField("Maximum number of super-resolution images in stack", getPrefs("maxSRImages", 50), 0);
    }

    @Override
    public boolean loadSettings() {
        maxSRImages = (int) gd.getNextNumber();
        setPrefs("maxSRImages", maxSRImages);
        prefs.savePreferences();
        return true;
    }

    @Override
    public void execute(){
    }


}
