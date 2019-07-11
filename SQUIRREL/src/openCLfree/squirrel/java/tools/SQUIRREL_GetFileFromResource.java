package openCLfree.squirrel.java.tools;

import java.io.*;
import java.util.Random;

/**
 * Created by sculley on 11/07/2019.
 */
public class SQUIRREL_GetFileFromResource {

    private SQUIRREL_GetFileFromResource() {
    }

    public static File getLocalFileFromResource(String path) throws IOException {
        if (!path.startsWith("/")) {
            throw new IllegalArgumentException("The path to be absolute (start with '/').");
        }

        // Obtain filename from path
        String[] parts = path.split("/");
        String filename = (parts.length > 1) ? parts[parts.length - 1] : null;

        File tDir = null;
        while (true) {
            int randomNum = (new Random()).nextInt(1000000000);
            tDir = new File(System.getProperty("java.io.tmpdir"), "NanoJ" + randomNum);
            if (tDir.exists() && !deleteDirectory(tDir)) {
                continue;
            }
            break;
        }
        tDir.mkdir();
        tDir.deleteOnExit();

        File temp = new File(tDir, filename);
        if (!temp.exists()) temp.createNewFile();
        else return temp;

        if (!temp.exists()) {
            throw new FileNotFoundException("File " + temp.getAbsolutePath() + " does not exist.");
        }

        // Prepare buffer for data copying
        byte[] buffer = new byte[1024];
        int readBytes;

        // Open and check input stream
        InputStream is = SQUIRREL_GetFileFromResource.class.getResourceAsStream(path);
        if (is == null) {
            throw new FileNotFoundException("File " + path + " was not found inside JAR.");
        }

        // Open output stream and copy data between source file in JAR and the temporary file
        OutputStream os = new FileOutputStream(temp);
        try {
            while ((readBytes = is.read(buffer)) != -1) {
                os.write(buffer, 0, readBytes);
            }
        } finally {
            // If read/write fails, close streams safely before throwing an exception
            os.close();
            is.close();
        }

        temp.deleteOnExit();
        return temp;
    }

    public static boolean deleteDirectory(File dir) {
        if (dir.isDirectory()) {
            File[] children = dir.listFiles();
            for (int i = 0; i < children.length; i++) {
                boolean success = deleteDirectory(children[i]);
                if (!success) {
                    return false;
                }
            }
        }
        // either file or an empty directory
        System.out.println("removing file or directory : " + dir.getName());
        return dir.delete();
    }
}
