package nanoj.squirrel.java.gui.testing;

import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import nanoj.core.java.gui._BaseDialog_;
import nanoj.core.java.threading.NanoJThreadExecutor;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by sculley on 18/01/2019.
 */
public class RollingBlockTest_ extends _BaseDialog_ {

    ArrayList<String> titles = new ArrayList<String>();
    String[] imageTitles;

    String titleRefImage = "", titleSRImage = "";

    ImagePlus impRef, impSR;
    ImageStack imsRef, imsSR;
    int nSlicesRef, nSlicesSR;
    int w_SR, h_SR, w_Ref, h_Ref;
    int magnification;

    @Override
    public boolean beforeSetupDialog(String s) {
        autoOpenImp = false;
        useSettingsObserver = true;
        int nImages = WindowManager.getImageCount();
        if(nImages<2){
            log.error("At least 2 images are required to run this code!");
            return false;
        }

        imageTitles = WindowManager.getImageTitles();
        for(int n=0; n<nImages; n++){
            ImagePlus thisImp = WindowManager.getImage(imageTitles[n]);
            titles.add(thisImp.getTitle());
        }

        imageTitles = new String[titles.size()];
        for(int n=0; n<titles.size(); n++){
            imageTitles[n] = titles.get(n);
        }

        return true;
    }

    @Override
    public void setupDialog() {

        gd = new NonBlockingGenericDialog("Test blocking");

        if (titleRefImage.equals("")) {
            titleRefImage = getPrefs("titleRefImage", "");
            titleSRImage = getPrefs("titleSRImage", "");
        }

        if (!titleRefImage.equals("") && contains(titleRefImage, imageTitles))
            gd.addChoice("Reference Image", imageTitles, titleRefImage);
        else
            gd.addChoice("Reference Image", imageTitles, imageTitles[0]);

        if (!titleSRImage.equals("") && contains(titleSRImage, imageTitles))
            gd.addChoice("Super-Resolution Reconstruction(s)", imageTitles, titleSRImage);
        else
            gd.addChoice("Super-Resolution Reconstruction(s)", imageTitles, imageTitles[0]);

    }

    @Override
    public boolean loadSettings() {

        titleRefImage = gd.getNextChoice();
        titleSRImage = gd.getNextChoice();

        setPrefs("titleRefImage", titleRefImage);
        setPrefs("titleSRImage", titleSRImage);

        impRef = WindowManager.getImage(titleRefImage);
        impSR = WindowManager.getImage(titleSRImage);

        imsSR = impSR.getImageStack();
        imsRef = impRef.getImageStack();

        w_SR = imsSR.getWidth();
        h_SR = imsSR.getHeight();

        w_Ref = imsRef.getWidth();
        h_Ref = imsRef.getHeight();

        magnification = w_SR / w_Ref;

        nSlicesSR = impSR.getStackSize();
        nSlicesRef = impRef.getStackSize();

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        FloatProcessor fpSR = imsSR.getProcessor(1).convertToFloatProcessor();
        FloatProcessor fpRef = imsRef.getProcessor(1).convertToFloatProcessor();

        float maxSigmaBoundary = 7.5f;

        int blockSize = (int) (5*maxSigmaBoundary);
        log.msg("blockSize is "+blockSize);
        while(blockSize%magnification!=0){
            blockSize++;
        }
        log.msg("new blockSize is "+blockSize);

        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);

        ImageStack imsSRBlocks = new ImageStack(blockSize,blockSize,w_SR*h_SR);
        ImageStack imsRefBlocks = new ImageStack(blockSize/magnification, blockSize/magnification, w_SR*h_SR);

        for (int y = 0; y < h_SR; y++) {
            for (int x = 0; x < w_SR; x++) {
                NTE.execute(new SplitIntoBlocks(fpSR, fpRef, x, y, blockSize, imsSRBlocks, imsRefBlocks));
            }
        }

        NTE.finish();

        new ImagePlus("SR blocks", imsSRBlocks).show();
        new ImagePlus("Ref blocks", imsRefBlocks).show();

    }

    //helper function lifted from imageshelper
    protected static boolean contains(String key, String[] strings) {
        for (String s: strings) {
            if (s.equals(key)) return true;
        }
        return false;
    }

    class SplitIntoBlocks extends Thread{

        private FloatProcessor fpSR, fpRef;
        private int x, y, blockSize;
        private ImageStack imsSRBlocks, imsRefBlocks;

        public SplitIntoBlocks(FloatProcessor fpSR, FloatProcessor fpRef, int x, int y,
                               int blockSize, ImageStack imsSRBlocks, ImageStack imsRefBlocks){
            this.fpSR = fpSR;
            this.fpRef = fpRef;
            this.x = x;
            this.y = y;
            this.blockSize = blockSize;
            this.imsSRBlocks = imsSRBlocks;
            this.imsRefBlocks = imsRefBlocks;
        }

        public void run(){
            int offset = blockSize/2;
            int blockSizeRef = blockSize/magnification;

            float[] pixelsSRBlock = new float[blockSize*blockSize];
            float[] pixelsRefBlock = new float[blockSizeRef*blockSizeRef];

            int p = y*w_SR + x;

            for(int j=0; j<blockSize; j++){
                int y_ = (j+y)-offset;
                if(y_<0 || y_>(h_SR-1)) continue;

                for(int i=0; i<blockSize; i++){
                    int x_ = (i+x)-offset;
                    if(x_<0 || x_>(w_SR-1)) continue;

                    int k = j*blockSize + i;

                    pixelsSRBlock[k] = fpSR.getf(x_,y_);

                    int jRef = j/magnification;
                    int iRef = i/magnification;
                    int yRef_ = y_/magnification;
                    int xRef_ = x_/magnification;

                    int kRef = jRef*blockSizeRef + iRef;

                    pixelsRefBlock[kRef] = fpRef.getf(xRef_, yRef_);
                }
            }

            imsSRBlocks.setProcessor(new FloatProcessor(blockSize, blockSize, pixelsSRBlock), p+1);
            imsRefBlocks.setProcessor(new FloatProcessor(blockSizeRef, blockSizeRef, pixelsRefBlock), p+1);

        }



    }
}
