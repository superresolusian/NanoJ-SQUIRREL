package nanoj.squirrel.java;

import ij.Prefs;
import nanoj.core.java._Version_;

/**
 * Created with IntelliJ IDEA.
 * User: Ricardo Henriques <paxcalpt@gmail.com>
 * Date: 12/12/14
 * Time: 16:33
 */

public class Version extends _Version_ {

    protected final static String tagInPrefs = "NJ.SQUIRREL.version";
    protected final static int major = 1;
    protected final static int minor = 1;
    protected final static int status = 1; // 0 - alpha, 1 - beta, 2 - release candidate, 3 - stable
    protected final static int release = 3;
    protected final static String codename = "";
    public final static String header = "NanoJ-SQUIRREL: ";
    public final static String WHATS_NEW =
            nanoj.core.java.Version.WHATS_NEW +  "\n \n" +
                    "What's new in NanoJ-SQUIRREL " + headlessGetVersion() + ":\n" +
                    "- Stability improvement for larger datasets...\n";

    public static String getVersion() {
        return header+headlessGetVersion();
    }

    public static String headlessGetVersion() {
        String text = major+"."+minor+getStatus()+release;
        if (codename!=null && !codename.equals(""))
            text += " \""+codename+"\"";
        return text;
    }

    public static String headlessGetVersionSmall() {
        return  major+"."+minor;
    }

    public static boolean isNewVersion() {
        if (!headlessGetVersion().equals(Prefs.get(tagInPrefs, ""))) {
            Prefs.set(tagInPrefs, headlessGetVersion());
            Prefs.savePreferences();
            return true;
        }
        return false;
    }

    public static String getStatus(){
        switch (status) {
            case 0: return "Alpha";
            case 1: return "Beta";
            case 2: return "RC";
            case 3: return "Stable";
        }
        return "";
    }
}

