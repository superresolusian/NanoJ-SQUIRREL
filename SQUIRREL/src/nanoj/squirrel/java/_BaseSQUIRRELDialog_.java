package nanoj.squirrel.java;

import com.boxysystems.jgoogleanalytics.FocusPoint;
import com.boxysystems.jgoogleanalytics.JGoogleAnalyticsTracker;
import nanoj.core.java.Version;
import nanoj.core.java.gui._BaseDialog_;

import static nanoj.core.java.Version.headlessGetVersionSmall;

/**
 * Created by Henriques-lab on 24/06/2017.
 */
public abstract class _BaseSQUIRRELDialog_ extends _BaseDialog_{

    protected JGoogleAnalyticsTracker tracker = new JGoogleAnalyticsTracker("NanoJ-SQUIRREL", headlessGetVersionSmall(), "UA-61590656-3");
    protected FocusPoint parentFocusPoint;

    protected void track(String value) {
        if (parentFocusPoint==null) parentFocusPoint = new FocusPoint(getClassName());
        FocusPoint focus = new FocusPoint(value);
        focus.setParentTrackPoint(parentFocusPoint);
        tracker.trackAsynchronously(focus);
    }

    protected boolean isNewVersion() {
        return Version.isNewVersion() || nanoj.squirrel.java.Version.isNewVersion();
    }

    protected String getWhatsNew() {
        return nanoj.squirrel.java.Version.WHATS_NEW;
    }
}
