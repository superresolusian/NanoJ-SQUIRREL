package openCLfree.squirrel.java.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.FolderOpener;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.File;
import java.io.IOException;


/**
 * Created by sculley on 11/07/2019.
 */
public abstract class _BaseSQUIRRELDialogNoCL_ extends Thread implements PlugIn {

    public NonBlockingGenericDialog gd = null;
    protected Font headerFont = new Font("Arial", Font.BOLD, 16);

    public ImagePlus imp = null;
    protected ImagePlus impPreview = null;
    protected boolean showPreview = false;

    private boolean firstRun = true;

    protected boolean useSettingsObserver = false;
    protected boolean autoOpenImp = true;

    public String prefsHeader = null, impPath = null;
    protected Prefs prefs = new Prefs();
    protected RefreshViewObserver RVO = new RefreshViewObserver();
    protected String arg;

    protected long typingWaitDelay = 500;

    public void run() {
        run("");
    }

    @Override
    public void run(String arg) {
        this.arg = arg;

        if (!beforeSetupDialog(arg)) return;
        if (autoOpenImp && imp == null) {
            setupImp();
            if (imp == null) return;
        }

        setupDialog();

        if (gd != null) {
            // Add listener to dialog
            MyDialogListener dl = new MyDialogListener();
            gd.addDialogListener(dl);

            if (useSettingsObserver) RVO.start();
            // Show dialog
            gd.showDialog();
            if (useSettingsObserver) RVO.quit();

            if (gd.wasCanceled()) {
                return;
            }
        } else {
            if (!loadSettings()) return;
        }

        afterDoPreview();
        if (!loadSettingsFromPrefs()) return; // also called within the RVO

        try {
            if (impPreview != null) {
                impPreview.close();
                impPreview = null;
            }
            execute();
        } catch (InterruptedException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        prefs.savePreferences();
    }

    abstract public boolean beforeSetupDialog(String arg);

    abstract public void setupDialog();

    abstract public boolean loadSettings();

    public boolean loadSettingsFromPrefs() {
        return true;
    }

    public void doPreview() {
    }

    public void afterDoPreview() {
        if (impPreview != null)
            impPreview.hide();
    }

    abstract public void execute() throws InterruptedException, IOException;

    protected void setupImp() {
        //get open image processor
        imp = WindowManager.getCurrentImage();

        if (imp == null) {
            openImp();
        }
        if (imp == null) return;

        if (imp.isComposite())
            IJ.log("WARNING: detected composite image. NanoJ is optimized for pure image stacks and may yield errors.");
        if (imp.isHyperStack())
            IJ.log("WARNING: detected hyperstack. NanoJ is optimized for pure image stacks and may yield errors.");
    }


    protected void openImp() {

        impPath = IJ.getFilePath("Choose data to load...");
        if (impPath == null || impPath.equals("")) return;

        else if (impPath.endsWith(".tif")) {
            if (impPath.contains("00")) {
                //open image sequence in directory selected and show
                imp = FolderOpener.open(new File(impPath).getParent());
            } else {
                imp = IJ.openImage(impPath);
            }
        } else {
            imp = IJ.openImage(impPath);
        }
        if (imp == null) return;
        imp.show();
    }

    class MyDialogListener implements DialogListener {
        @Override
        public boolean dialogItemChanged(GenericDialog gd, AWTEvent awtEvent) {
            if (!loadSettings()) return false;
            RVO.needsUpdate = true;
            return true;
        }
    }

    private class RefreshViewObserver extends Thread {

        boolean _stop = false;
        boolean needsUpdate = false;

        int nslice = -1;

        public void start() {
            super.start();
        }

        public void quit() {
            _stop = true;
        }

        public void run() {
            while (!_stop) {

                boolean changedSlice = false;
                // ~~ See if something changed ~~

                if (imp != null && imp.getSlice() != nslice) { // if the user changed the current frame in the image
                    needsUpdate = true;
                    changedSlice = true;
                    nslice = imp.getSlice();
                }
                // ~~ Act on change ~~
                if (needsUpdate) { // the user is typing, give the boss time to finish
                    if (firstRun) {
                        try {
                            sleep(100);
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                        firstRun = false;
                    } else if (changedSlice) {
                        try {
                            sleep(3000);
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    } else {
                        try {
                            sleep(typingWaitDelay);
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    }

                    if (!loadSettingsFromPrefs()) gd.getButtons()[0].setEnabled(false);
                    else gd.getButtons()[0].setEnabled(true);
                    updateView();
                    needsUpdate = false;
                }
                // give it 30ms until reprocessing
                try {
                    sleep(30);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

        public void updateView() {
            if (showPreview) {
                String title = gd.getTitle();
                gd.setTitle("(calculating)... " + title);
                doPreview();
                gd.setTitle(title);
            }
        }
    }

    public String getClassName() {
        return this.getClass().getName();
    }

    public int getPrefs(String key, int defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (int) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public float getPrefs(String key, float defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (float) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public double getPrefs(String key, double defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return (double) prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public boolean getPrefs(String key, boolean defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public String getPrefs(String key, String defaultValue) {
        if (prefsHeader == null) prefsHeader = getClassName();
        return prefs.get(prefsHeader + "." + key, defaultValue);
    }

    public void setPrefs(String key, int value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, float value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, double value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, boolean value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void setPrefs(String key, String value) {
        if (prefsHeader == null) prefsHeader = getClassName();
        prefs.set(prefsHeader + "." + key, value);
    }

    public void savePrefs() {
        prefs.savePreferences();
    }
}
