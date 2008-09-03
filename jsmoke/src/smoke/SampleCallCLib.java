package smoke;
import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
/**
 *
 * @author huub
 */
public class SampleCallCLib {    
    private static MyNativeLibrary INSTANCE;
    
    static {
        INSTANCE = (MyNativeLibrary)Native.loadLibrary("c",MyNativeLibrary.class);
    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.err.println(INSTANCE.atol("123"));
    }
    
}

interface MyNativeLibrary extends Library {
    public int atol(String s);
}

