package nanoj.squirrel.java;

import com.jogamp.opencl.*;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.io.IOException;
import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;


public class ISuckAtCL {

    public boolean DEBUG = false;

    static CLContext context;
    static CLProgram program;
    static CLKernel kernelJustAdd;
    static CLCommandQueue queue;

    private final int width, height;

    private CLBuffer<FloatBuffer> clBufferPx1, clBufferPx2, clBufferSum;

    public ISuckAtCL(int width, int height){

        this.width = width;
        this.height = height;

        context = CLContext.create();
        System.out.println("created " + context);

        CLDevice device = context.getMaxFlopsDevice();
        System.out.println("using " + device);

        queue = device.createCommandQueue();

        try{
            program = context.createProgram(ISuckAtCL.class.getResourceAsStream("/JustAdd.cl")).build();
        } catch (IOException e) {
            e.printStackTrace();
        }

        kernelJustAdd = program.createCLKernel("justAdd");

        clBufferPx1 = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferPx2 = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferSum = context.createFloatBuffer(width*height, READ_ONLY);
    }

    public synchronized FloatProcessor justAdd(ImageProcessor ip1, ImageProcessor ip2){
        FloatProcessor fpPx1 = ip1.convertToFloatProcessor();
        FloatProcessor fpPx2 = ip2.convertToFloatProcessor();
        FloatProcessor fpSum = new FloatProcessor(width,height);

        fillBuffer(clBufferPx1, fpPx1);
        fillBuffer(clBufferPx2, fpPx2);

        int n = 0;
        kernelJustAdd.setArg(n++, clBufferPx1 ); // make sure type is the same !!
        kernelJustAdd.setArg(n++, clBufferPx2 ); // make sure type is the same !!
        kernelJustAdd.setArg(n++, clBufferSum);
        kernelJustAdd.setArg(n++, width);
        kernelJustAdd.setArg(n++, height);

        queue.putWriteBuffer(clBufferPx1, false );
        queue.putWriteBuffer(clBufferPx2, false );

        queue.put2DRangeKernel(kernelJustAdd, 0, 0, width, height, 0, 0);
        queue.finish();

        queue.putReadBuffer(clBufferSum, true);
        grabBuffer(clBufferSum, fpSum);

        return fpSum;
    }

    public static void fillBuffer(CLBuffer<FloatBuffer> clBuffer, ImageProcessor ip) {
        FloatBuffer buffer = clBuffer.getBuffer();
        for(int n=0; n<ip.getPixelCount(); n++) buffer.put(n, ip.getf(n));
    }

    public static void grabBuffer(CLBuffer<FloatBuffer> clBuffer, FloatProcessor fp) {
        FloatBuffer buffer = clBuffer.getBuffer();
        for(int n=0; n<fp.getPixelCount(); n++) fp.setf(n, buffer.get(n));
    }


}
