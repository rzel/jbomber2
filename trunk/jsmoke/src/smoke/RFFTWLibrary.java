package smoke;

import com.sun.jna.Library;
import com.sun.jna.Pointer;
import com.sun.jna.Structure;

public interface RFFTWLibrary extends Library {
    int FFTW_REAL_TO_COMPLEX=-1; // FFTW_FORWARD    
    int FFTW_COMPLEX_TO_REAL= 1; // FFTW_BACKWARD
    int FFTW_IN_PLACE=8;
    
    Pointer rfftw2d_create_plan(int nx, int ny, int dir, int flags);
    void rfftwnd_one_real_to_complex_in_place(Pointer p, double[] u);
    void rfftwnd_one_complex_to_real_in_place(Pointer p, double[] u);
}