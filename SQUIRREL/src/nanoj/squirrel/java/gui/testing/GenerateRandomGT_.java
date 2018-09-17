package nanoj.squirrel.java.gui.testing;

import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import ij.process.FloatProcessor;
import nanoj.core.java.gui._BaseDialog_;

import java.awt.*;
import java.io.IOException;
import java.util.Random;

import static java.lang.Math.*;

/**
 * Created by sculley on 14/09/2018.
 */
public class GenerateRandomGT_ extends _BaseDialog_ {

    String[] structureOptions = new String[]{"0-D (dots)", "1-D (lines)", "2-D (patches)"};
    int width;
    String structure;
    double coverage;

    Random random = new Random();

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = false;
        useSettingsObserver = true;
        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Generate random ground truth test structure");

        gd.addNumericField("Image width (pixels)", getPrefs("width", 500), 0);
        gd.addRadioButtonGroup("Type of structure", structureOptions, 1, 3, getPrefs("structure", structureOptions[0]));
        gd.addNumericField("Approximate fraction of image occupied by structure", getPrefs("coverage", 0.05), 2);
    }

    @Override
    public boolean loadSettings() {
        width = (int) gd.getNextNumber();
        structure = gd.getNextRadioButton();
        coverage = gd.getNextNumber();

        setPrefs("width", width);
        setPrefs("structure", structure);
        setPrefs("coverage", coverage);

        return true;
    }

    @Override
    public void execute() throws InterruptedException, IOException {

        float[] pixels = new float[width*width];
        FloatProcessor fp = new FloatProcessor(width, width, pixels);
        fp.setColor(1);
        double target = coverage*width*width;

        int structurePixels = 0;

        //0-D case - random dots
        if(structure.equals(structureOptions[0])){

            while(structurePixels<target){
                int x = (int) Math.floor(random.nextDouble()*width);
                int y = (int) Math.floor(random.nextDouble()*width);
                fp.setf(x, y, 1.0f);
                structurePixels = getSum(fp);
            }
        }
        // 1-D case - random lines
        else if(structure.equals(structureOptions[1])){

            while(structurePixels<target) {
                int x1 = (int) Math.floor(random.nextDouble() * width);
                int x2 = (int) Math.floor(random.nextDouble() * width);
                int y1 = (int) Math.floor(random.nextDouble() * width);
                int y2 = (int) Math.floor(random.nextDouble() * width);

                fp.drawLine(x1, y1, x2, y2);

                structurePixels = getSum(fp);
            }
        }
        //2-D case - random patches
        else if(structure.equals(structureOptions[2])){

            while(structurePixels<target){
                int x0 = (int) Math.floor(random.nextDouble()*width);
                int y0 = (int) Math.floor(random.nextDouble()*width);

                int nVertices = (int) Math.floor(random.nextDouble()*9)+3;
                double maxRadius = random.nextDouble()*width*0.05;

                double[] angles = new double[nVertices];
                for(int i=0; i<nVertices; i++){
                    angles[i] = random.nextDouble()*2*PI;
                }
                angles = doSort(angles);
                int[] xVals = new int[nVertices];
                int[] yVals = new int[nVertices];

                for(int i=0; i<nVertices; i++){
                    double thisRadius = maxRadius*random.nextDouble();
                    double thisAngle = angles[i];

                    int y = (int) (thisRadius*cos(thisAngle)) + y0;
                    int x = (int) (thisRadius*sin(thisAngle)) + x0;

                    xVals[i] = x;
                    yVals[i] = y;
                }

                Polygon p = new Polygon(xVals, yVals, nVertices);
                fp.fillPolygon(p);

                structurePixels = getSum(fp);
            }

        }

        double actualCoverage = (double) structurePixels/(width*width);

        new ImagePlus("GT - true coverage="+actualCoverage, fp).show();

    }

    int getSum(FloatProcessor fp){
        float[] pixels = (float[]) fp.getPixels();
        float v=0;
        for(int i=0; i<pixels.length; i++){
            v += pixels[i];
        }
        return (int) v;
    }

    double[] doSort(double[] array){
        int i = 1;
        double x;
        int j;

        while(i < array.length){
            x = array[i];
            j = i-1;
            while(j>=0 && array[j] > x){
                array[j+1] = array[j];
                j--;
            }
            array[j+1] = x;
            i++;
        }
        return array;
    }
}
