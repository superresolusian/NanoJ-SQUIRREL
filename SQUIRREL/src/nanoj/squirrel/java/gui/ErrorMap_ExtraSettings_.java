package nanoj.squirrel.java.gui;

import ij.gui.NonBlockingGenericDialog;
import nanoj.squirrel.java._BaseSQUIRRELDialog_;

import java.awt.*;
import java.io.IOException;

/**
 * Created by sculley on 14/06/2017.
 */
public class ErrorMap_ExtraSettings_ extends _BaseSQUIRRELDialog_ {

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = false;
        useSettingsObserver = true;
        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Calculate Error Map Advanced Settings...");
        gd.hideCancelButton();

        gd.addMessage("-=-= Frame purge =-=-");
        gd.addCheckbox("Check and purge empty frames?", getPrefs("framePurge", false));

        gd.addMessage("-=-= Image Resizing =-=-\n", headerFont);
        gd.addCheckbox("Crop black borders from super-resolution image (default: active)", getPrefs("borderControl", true));

        gd.addMessage("-=-= Registration Options =-=-\n", headerFont);
        gd.addCheckbox("Enable registration (default: active)", getPrefs("doRegistration", true));
        gd.addNumericField("Maximum expected misalignment (0-auto)", getPrefs("maxExpectedMisalignment", 0), 0);

        gd.addMessage("-=-= Output Image Options =-=-\n", headerFont);
        gd.addCheckbox("Show intensity-normalised and cropped super-resolution image(s) (default:active)", getPrefs("showIntensityNormalised", true));
        gd.addMessage("The above output is required for image fusion", new Font("Arial", Font.ITALIC, 12));
        gd.addCheckbox("Show_RSF-convolved super-resolution image(s) (default: active)", getPrefs("showConvolved", true));
        gd.addCheckbox("Show_RSF image(s) (default: disabled)", getPrefs("showRSF", false));
        gd.addCheckbox("Show_positive and negative contributions to error map (default: disabled)", getPrefs("showPositiveNegative", false));
        gd.addMessage("If the above checkbox is disabled, error map will just contain absolute values of errors.", new Font("Arial", Font.ITALIC, 12));


    }

    public boolean loadSettings(){

        boolean framePurge = gd.getNextBoolean();
        boolean borderControl = gd.getNextBoolean();

        boolean doRegistration = gd.getNextBoolean();
        int maxExpectedMisalignment = (int) gd.getNextNumber();

        boolean showIntensityNormalised = gd.getNextBoolean();
        boolean showConvolved = gd.getNextBoolean();
        boolean showRSF = gd.getNextBoolean();
        boolean showPositiveNegative = gd.getNextBoolean();

        setPrefs("framePurge", framePurge);
        setPrefs("borderControl", borderControl);

        setPrefs("doRegistration", doRegistration);
        setPrefs("maxExpectedMisalignment", maxExpectedMisalignment);

        setPrefs("showIntensityNormalised", showIntensityNormalised);
        setPrefs("showConvolved", showConvolved);
        setPrefs("showRSF", showRSF);
        setPrefs("showPositiveNegative", showPositiveNegative);

        prefs.savePreferences();

        return true;
    }

    public boolean loadSettingsFromPrefs(){
        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

    }

}
