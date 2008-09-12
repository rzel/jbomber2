import com.sun.jna.Native;
import com.sun.jna.Pointer;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.util.Vector;
import javax.swing.BoxLayout;
import javax.media.opengl.GL;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLCapabilities;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.GLCanvas;
import javax.media.opengl.glu.GLU;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import java.awt.event.*;
import java.awt.*;
import javax.swing.table.*;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.colorchooser.AbstractColorChooserPanel;

/** Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
 *        the velocity field at the mouse location. Press the indicated keys to change options
 */
class Smoke {
    static RFFTWLibrary RFFTW = (RFFTWLibrary)Native.loadLibrary("rfftw",RFFTWLibrary.class);

//--- SIMULATION PARAMETERS ------------------------------------------------------------------------
    // final int DIM = 16;		//size of simulation grid
    final int DIM = 50;		//size of simulation grid
    double dt = 0.4;		//simulation time step
    double visc = 0.001;	//fluid viscosity
    double/*fftw_real*/ vx [], vy  [];        //(vx,vy)   = velocity field at the current moment
    double/*fftw_real*/ vx0[], vy0 [];        //(vx0,vy0) = velocity field at the previous moment
    double/*fftw_real*/ fx [], fy  [];	//(fx,fy)   = user-controlled simulation forces, steered with the mouse
    double/*fftw_real*/ rho[], rho0[];	//smoke density at the current (rho) and previous (rho0) moment
    /*rfftwnd_plan*/ Pointer plan_rc, plan_cr;  //simulation domain discretization


//--- VISUALIZATION PARAMETERS ---------------------------------------------------------------------
    int   winWidth, winHeight;     //size of the graphics window, in pixels
    int   color_dir = 0;           //use direction color-coding or not
    float vec_scale = 1000;        //scaling of hedgehogs
    boolean   draw_smoke = true;  //draw the smoke or not
    boolean   draw_vecs = false;    //draw the vector field or not
    final int COLOR_BLACKWHITE=0;  //different types of color mapping: black-and-white, rainbow, banded
    final int COLOR_RAINBOW=1;
    final int COLOR_BANDS=2;
		final int COLOR_CUSTOM=3;      // wtf enum?
    int   scalar_col = COLOR_CUSTOM;           //method for scalar coloring
    boolean frozen = false;         //toggles on/off the animation
    final int DATASET_RHO = 0;
    final int DATASET_F = 1;
    final int DATASET_V = 2;
    int dataset = 0;                // Rho |V| of|f|
    final int SCALE_NONE = 0;
    final int SCALE_CLAMP = 1;
    final int SCALE_SCALE = 2;
    int scaleMode = SCALE_SCALE;
    int ncolors = 255;
    double minClamp = 0.0d;
    double maxClamp = 1.0d;

//------ SIMULATION CODE STARTS HERE -----------------------------------------------------------------

    /**init_simulation: Initialize simulation data structures as a function of the grid size 'n'.
     *                 Although the simulation takes place on a 2D grid, we allocate all data structures as 1D arrays,
     *                 for compatibility with the FFTW numerical library.
     */
    void init_simulation(int n) {
        int dim     = n * 2*(n/2+1);              //Allocate data structures
        vx       = new double/*fftw_real*/[dim];
        vy       = new double/*fftw_real*/[dim];
        vx0      = new double/*fftw_real*/[dim];
        vy0      = new double/*fftw_real*/[dim];
        for (int i = 0; i < dim; i++)                      //Initialize data structures to 0
        { vx[i] = vy[i] = vx0[i] = vy0[i] =0.0; } // float or double??

        dim     = n * n; // * sizeof(double/*fftw_real*/);
        fx      = new double/*fftw_real*/[dim];
        fy      = new double/*fftw_real*/[dim];
        rho     = new double/*fftw_real*/[dim];
        rho0    = new double/*fftw_real*/[dim];
        plan_rc = RFFTW.rfftw2d_create_plan(n, n, RFFTW.FFTW_REAL_TO_COMPLEX, RFFTW.FFTW_IN_PLACE);
        plan_cr = RFFTW.rfftw2d_create_plan(n, n, RFFTW.FFTW_COMPLEX_TO_REAL, RFFTW.FFTW_IN_PLACE);

        for (int i = 0; i < dim; i++)               //Initialize data structures to 0
        { fx[i] = fy[i] = rho[i] = rho0[i] = 0.0; } // float or double??
    }


//FFT: Execute the Fast Fourier Transform on the dataset 'vx'.
//     'direction' indicates if we do the direct (1) or inverse (-1) Fourier Transform
    void FFT(int direction,double/*fftw_real*/[] vx) {
        if(direction==1) RFFTW.rfftwnd_one_real_to_complex_in_place(plan_rc,(double[])vx);
        else             RFFTW.rfftwnd_one_complex_to_real_in_place(plan_cr,(double[])vx);
    }

    int clamp(double x) {
        return ((x)>=0.0?((int)(x)):(-((int)(1-(x))))); }

//solve: Solve (compute) one step of the fluid flow simulation
    void solve(int n, double/*fftw_real*/[] vx, double/*fftw_real*/[] vy, double/*fftw_real*/[] vx0, double/*fftw_real*/[] vy0, double/*fftw_real*/ visc, double dt) {
        double/*fftw_real*/ x, y, x0, y0, f, r, U[]=new double/*fftw_real*/[2], V[]=new double/*fftw_real*/[2], s, t;
        int i, j, i0, j0, i1, j1;

        for (i=0;i<n*n;i++) {
            vx[i] += dt*vx0[i]; vx0[i] = vx[i]; vy[i] += dt*vy0[i]; vy0[i] = vy[i]; }

        for ( x=0.5/n,i=0 ; i<n ; i++,x+=1.0/n ) {
            for ( y=0.5/n,j=0 ; j<n ; j++,y+=1.0/n ) {
                x0 = n*(x-dt*vx0[i+n*j])-0.5f;
                y0 = n*(y-dt*vy0[i+n*j])-0.5f;
                i0 = clamp(x0); s = x0-i0;
                i0 = (n+(i0%n))%n;
                i1 = (i0+1)%n;
                j0 = clamp(y0); t = y0-j0;
                j0 = (n+(j0%n))%n;
                j1 = (j0+1)%n;
                vx[i+n*j] = (1-s)*((1-t)*vx0[i0+n*j0]+t*vx0[i0+n*j1])+s*((1-t)*vx0[i1+n*j0]+t*vx0[i1+n*j1]);
                vy[i+n*j] = (1-s)*((1-t)*vy0[i0+n*j0]+t*vy0[i0+n*j1])+s*((1-t)*vy0[i1+n*j0]+t*vy0[i1+n*j1]);
            }
        }

        for(i=0; i<n; i++) {
            for(j=0; j<n; j++) {
                vx0[i+(n+2)*j] = vx[i+n*j]; vy0[i+(n+2)*j] = vy[i+n*j]; }
        }

        FFT(1,vx0);
        FFT(1,vy0);

        for (i=0;i<=n;i+=2) {
            x = 0.5*i;
            for (j=0;j<n;j++) {
                y = j<=n/2 ? (double)j : (double)(j-n);
                r = x*x+y*y;
                if ( r==0.0f ) continue;
                f = Math.exp(-r*dt*visc);
                U[0] = vx0[i  +(n+2)*j]; V[0] = vy0[i  +(n+2)*j];
                U[1] = vx0[i+1+(n+2)*j]; V[1] = vy0[i+1+(n+2)*j];

                vx0[i  +(n+2)*j] = f*((1-x*x/r)*U[0]     -x*y/r *V[0]);
                vx0[i+1+(n+2)*j] = f*((1-x*x/r)*U[1]     -x*y/r *V[1]);
                vy0[i+  (n+2)*j] = f*(  -y*x/r *U[0] + (1-y*y/r)*V[0]);
                vy0[i+1+(n+2)*j] = f*(  -y*x/r *U[1] + (1-y*y/r)*V[1]);
            }
        }

        FFT(-1,vx0);
        FFT(-1,vy0);

        f = 1.0/(n*n);
        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) {
                vx[i+n*j] = f*vx0[i+(n+2)*j]; vy[i+n*j] = f*vy0[i+(n+2)*j];
            }
        }
    }


// diffuse_matter: This function diffuses matter that has been placed in the velocity field. It's almost identical to the
// velocity diffusion step in the function above. The input matter densities are in rho0 and the result is written into rho.
    void diffuse_matter(int n, double/*fftw_real*/[] vx, double/*fftw_real*/[] vy, double/*fftw_real*/[] rho, double/*fftw_real*/[] rho0, double dt) {
        double/*fftw_real*/ x, y, x0, y0, s, t;
        int i, j, i0, j0, i1, j1;

        for ( x=0.5/n,i=0 ; i<n ; i++,x+=1.0/n ) {
            for ( y=0.5/n,j=0 ; j<n ; j++,y+=1.0/n ) {
                x0 = n*(x-dt*vx[i+n*j])-0.5;
                y0 = n*(y-dt*vy[i+n*j])-0.5;
                i0 = clamp(x0);
                s = x0-i0;
                i0 = (n+(i0%n))%n;
                i1 = (i0+1)%n;
                j0 = clamp(y0);
                t = y0-j0;
                j0 = (n+(j0%n))%n;
                j1 = (j0+1)%n;
                rho[i+n*j] = (1-s)*((1-t)*rho0[i0+n*j0]+t*rho0[i0+n*j1])+s*((1-t)*rho0[i1+n*j0]+t*rho0[i1+n*j1]);
            }
        }
    }

//set_forces: copy user-controlled forces to the force vectors that are sent to the solver.
//            Also dampen forces and matter density to get a stable simulation.
    void set_forces() {
        int i;
        for (i = 0; i < DIM * DIM; i++) {
            rho0[i]  = 0.995f * rho[i];
            fx[i] *= 0.85f;
            fy[i] *= 0.85f;
            vx0[i]    = fx[i];
            vy0[i]    = fy[i];
        }
    }


//do_one_simulation_step: Do one complete cycle of the simulation:
//      - set_forces:
//      - solve:            read forces from the user
//      - diffuse_matter:   compute a new set of velocities
//      - gluPostRedisplay: draw a new visualization frame
    void do_one_simulation_step() {
        if (!frozen) {
            set_forces();
            solve(DIM, vx, vy, vx0, vy0, visc, dt);

            diffuse_matter(DIM, vx, vy, rho, rho0, dt);
            panel.repaint(50);
        }
    }


//------ VISUALIZATION CODE STARTS HERE -----------------------------------------------------------------


    //rainbow: Implements a color palette, mapping the scalar 'value' to a rainbow color RGB
    void rainbow(float value,float[] color) {
        final float dx=0.8f;
        if (value<0) value=0; if (value>1) value=1;
        value = (6-2*dx)*value+dx;
        color[0] = Math.max(0.0f,(3-Math.abs(value-4)-Math.abs(value-5))/2);
        color[1] = Math.max(0.0f,(4-Math.abs(value-2)-Math.abs(value-4))/2);
        color[2] = Math.max(0.0f,(3-Math.abs(value-1)-Math.abs(value-2))/2);
    }

		private static float maxvy=1.0f;
		private static float maxvy_lastframe=1.0f;
		private static float minvy=0.0f;
		private static float minvy_lastframe=0.0f;

		void custom_gradient(float f, float[] rgb) {
//			maxvy = f > maxvy ? f : maxvy;
//			minvy = f < minvy ? f : minvy;

	//		f=(f-minvy_lastframe)/(maxvy_lastframe-minvy_lastframe); // Autoscale!
			f = f<0.0f ? 0.0f : f>1.0f ? 1.0f : f; // Clamp!

			int n = colortable.getRowCount();
			float r = (n-1) * f;

			Color a = (Color)colortable.getValueAt((int)r, 0);
			Color b = (Color)colortable.getValueAt(Math.min((int)r+1, n-1), 0);

			float h = r - (int)r;
			float l = 1 - h;
			rgb[0] = (a.getRed()   * l + b.getRed()   * h) / 255.0f;
			rgb[1] = (a.getGreen() * l + b.getGreen() * h) / 255.0f;
			rgb[2] = (a.getBlue()  * l + b.getBlue()  * h) / 255.0f;
		}

    //set_colormap: Sets three different types of colormaps
    void set_colormap(GL gl, float vy) {
        float[] rgb = new float[3];

        // Clamp
        if ((scaleMode & SCALE_CLAMP) == SCALE_CLAMP) {
            vy = (float)(vy<minClamp ? minClamp : vy>maxClamp ? maxClamp : vy); // Clamp!
        }

        maxvy = vy > maxvy ? vy : maxvy;
        minvy = vy < minvy ? vy : minvy;

        if ((scaleMode & SCALE_SCALE) == SCALE_SCALE) {
					vy = (vy - minvy_lastframe) / (maxvy_lastframe - minvy_lastframe);
        }

        // Band color count
        vy *= ncolors; vy = (int)(vy); vy/= ncolors;

        if (scalar_col==COLOR_BLACKWHITE)
            rgb[0]=rgb[1]=rgb[2] = vy;
        else if (scalar_col==COLOR_RAINBOW)
            rainbow(vy,rgb);
        else if (scalar_col==COLOR_BANDS) {
            final int NLEVELS = 7;
            vy *= NLEVELS; vy = (int)(vy); vy/= NLEVELS;
            rainbow(vy,rgb);
        }
				else if(scalar_col==COLOR_CUSTOM) {
					custom_gradient(vy, rgb);
				}

        gl.glColor3f(rgb[0], rgb[1], rgb[2]);
    }


    //direction_to_color: Set the current color by mapping a direction vector (x,y), using
    //                    the color mapping method 'method'. If method==1, map the vector direction
    //                    using a rainbow colormap. If method==0, simply use the white color
    void direction_to_color(GL gl, float x, float y, int method) {
        float r,g,b;
        if (method==1) {
            float f = (float)(Math.atan2(y,x) / Math.PI + 1);
            r = f;
            if(r > 1) r = 2 - r;
            g = f + .66667f;
            if(g > 2) g -= 2;
            if(g > 1) g = 2 - g;
            b = f + 2 * .66667f;
            if(b > 2) b -= 2;
            if(b > 1) b = 2 - b;
        } else { // method==0
            r = g = b = 1; }
        gl.glColor3f(r,g,b);
    }

    private float getDatasetColor(int idx) {
         switch (dataset) {
            case DATASET_RHO:
                return (float)rho[idx];
            case DATASET_F:
                return (float)Math.sqrt(fx[idx] * fx[idx] + fy[idx] * fy[idx]) * DIM;
            case DATASET_V:
                return (float)Math.sqrt(vx[idx] * vx[idx] + vy[idx] * vy[idx]) * DIM;
        }

        return 0.0f;
    }

    //visualize: This is the main visualization function
    void visualize(GL gl) {
        int        i, j, idx; double px,py;
        double/*fftw_real*/  wn = winWidth / (double)(DIM + 1);   // Grid cell width
        double/*fftw_real*/  hn = winHeight / (double)(DIM + 1);  // Grid cell heigh

        if (draw_smoke) {
            gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
            for (j = 0; j < DIM - 1; j++)			//draw smoke
            {
                gl.glBegin(GL.GL_TRIANGLE_STRIP);

                i = 0;
                px = wn + (float)i * wn;
                py = hn + (float)j * hn;
                idx = (j * DIM) + i;
                gl.glColor3d(rho[idx],rho[idx],rho[idx]);
                gl.glVertex2d(px,py);

                for (i = 0; i < DIM - 1; i++) {
                    px = wn + i * wn;
                    py = hn + (j + 1) * hn;
                    idx = ((j + 1) * DIM) + i;

                    set_colormap(gl, getDatasetColor(idx));
                    gl.glVertex2d(px, py);
                    px = wn + (i + 1) * wn;
                    py = hn + j * hn;
                    idx = (j * DIM) + (i + 1);
                    set_colormap(gl, getDatasetColor(idx));
                    gl.glVertex2d(px, py);
                }

                px = wn + (float)(DIM - 1) * wn;
                py = hn + (float)(j + 1) * hn;
                idx = ((j + 1) * DIM) + (DIM - 1);
                set_colormap(gl, getDatasetColor(idx));
                gl.glVertex2d(px, py);
                gl.glEnd();
            }
        }

        if (draw_vecs) {
            gl.glBegin(GL.GL_LINES);				//draw velocities
            for (i = 0; i < DIM; i++)
                for (j = 0; j < DIM; j++) {
                idx = (j * DIM) + i;
                direction_to_color(gl, (float)(double)vx[idx],(float)(double)vy[idx],color_dir);
                gl.glVertex2d(wn + i * wn, hn + j * hn);
                gl.glVertex2d((wn + i * wn) + vec_scale * vx[idx], (hn + j * hn) + vec_scale * vy[idx]);
                }
            gl.glEnd();
        }
        gl.glFlush(); // forces all opengl commands to complete. Blocking!!
				maxvy_lastframe = maxvy;
				minvy_lastframe = minvy;
				maxvy = Float.MIN_VALUE;
				minvy = Float.MAX_VALUE;
    }




    //display: Handle window redrawing events. Simply delegates to visualize().
    void display(GL gl) {
        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
        gl.glMatrixMode(GL.GL_MODELVIEW);
        gl.glLoadIdentity();
        visualize(gl);

//        /** draw white diagonal. */
//        gl.glBegin(GL.GL_LINES);
//        gl.glColor3f(1,1,1);
//        gl.glVertex2f(0,0);
//        gl.glVertex2f(100,100);
//        gl.glEnd();
//        gl.glFlush();
//
        //glutSwapBuffers();
    }

    //reshape: Handle window resizing (reshaping) events
    void reshape(GL gl, int w, int h) {
        GLU glu = new GLU();
        gl.glViewport(0, 0,w, h);
        gl.glMatrixMode(GL.GL_PROJECTION);
        gl.glLoadIdentity();
        glu.gluOrtho2D(0.0, (double)w, 0.0, (double)h);
        winWidth = w; winHeight = h;
    }

//------ INTERACTION CODE STARTS HERE -----------------------------------------------------------------

//keyboard: Handle key presses
    void keyboard(char key, int x, int y) {
        switch (key) {
            case 't': dt -= 0.001; break;
            case 'T': dt += 0.001; break;
            case 'c': color_dir = 1 - color_dir; break;
            case 'S': vec_scale *= 1.2; break;
            case 's': vec_scale *= 0.8; break;
            case 'V': visc *= 5; break;
            case 'v': visc *= 0.2; break;
            case 'x': draw_smoke = !draw_smoke;
            if (!draw_smoke) draw_vecs = true; break;
            case 'y': draw_vecs = !draw_vecs;
            if (!draw_vecs) draw_smoke = true; break;
            case 'm': scalar_col++; if (scalar_col>COLOR_CUSTOM) scalar_col=COLOR_BLACKWHITE; break;
            case 'a': frozen = !frozen; break;
            case 'q': System.exit(0);
        }
    }



// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
//       cursor movement. Also inject some new matter into the field at the mouse location.
    int lmx=0,lmy=0;				//remembers last mouse location
    void drag(int mx, int my) {
        int xi,yi,X,Y; double  dx, dy, len;

        // Compute the array index that corresponds to the cursor location
        xi = (int)clamp((double)(DIM + 1) * ((double)mx / (double)winWidth));
        yi = (int)clamp((double)(DIM + 1) * ((double)(winHeight - my) / (double)winHeight));

        X = xi; Y = yi;

        if (X > (DIM - 1))  X = DIM - 1; if (Y > (DIM - 1))  Y = DIM - 1;
        if (X < 0) X = 0; if (Y < 0) Y = 0;

        // Add force at the cursor location
        my = winHeight - my;
        dx = mx - lmx; dy = my - lmy;
        len = Math.sqrt(dx * dx + dy * dy);
        if (len != 0.0) {  dx *= 0.1 / len; dy *= 0.1 / len; }
        fx[Y * DIM + X] += (float)dx;
        fy[Y * DIM + X] += (float)dy;
        rho[Y * DIM + X] = 10.0d;
        lmx = mx; lmy = my;
    }


    //main: The main program
    public static void main(String[] argv) {
        System.out.println("Fluid Flow Simulation and Visualization");
        System.out.println("=======================================");
        System.out.println("Click and drag the mouse to steer the flow!");
        System.out.println("T/t:   increase/decrease simulation timestep");
        System.out.println("S/s:   increase/decrease hedgehog scaling");
        System.out.println("c:     toggle direction coloring on/off");
        System.out.println("V/v:   increase decrease fluid viscosity");
        System.out.println("x:     toggle drawing matter on/off");
        System.out.println("y:     toggle drawing hedgehogs on/off");
        System.out.println("m:     toggle thru scalar coloring");
        System.out.println("a:     toggle the animation on/off\n");
        System.out.println("q:     quit");

        new Smoke();
    }

    class DatasetSelectorListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
            if (e.getActionCommand().equals("DATASET_RHO"))
            {
                dataset = DATASET_RHO;
            }
            else if (e.getActionCommand().equals("DATASET_F"))
            {
                dataset = DATASET_F;
            }
            else if  (e.getActionCommand().equals("DATASET_V"))
            {
                dataset = DATASET_V;
            }
            else
            {
                System.out.println("Dataset: " + e.getActionCommand());
            }
        }
    }

    class ColorSelectorListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
            if (e.getActionCommand().equals("COLORMAP_RAINBOW"))
            {
                scalar_col = COLOR_RAINBOW;
            }
            else if (e.getActionCommand().equals("COLORMAP_GRAYSCALE"))
            {
                scalar_col = COLOR_BLACKWHITE;
            }
            else if  (e.getActionCommand().equals("COLORMAP_DEFINED"))
            {
                scalar_col = COLOR_BANDS;
            }
            else if  (e.getActionCommand().equals("COLORMAP_CUSTOM"))
            {
                scalar_col = COLOR_CUSTOM;
            }
            else if  (e.getActionCommand().equals("SIMULATION_ON"))
            {
                frozen = false;
                do_one_simulation_step();
            }
            else if  (e.getActionCommand().equals("SIMULATION_OFF"))
            {
                frozen = true;
            }
            else
            {
                System.out.println("Colormap: " + e.getActionCommand());
            }
        }
    }

    class SmokeSelectorListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
            if (e.getActionCommand().equals("SMOKE_TOGGLE"))
            {
                draw_smoke = !draw_smoke;
            }
            else if (e.getActionCommand().equals("VECTOR_TOGGLE"))
            {
                draw_vecs = !draw_vecs;
            }
            else
            {
                System.out.println("Smoke: " + e.getActionCommand());
            }
        }
    }

    private JPanel initDatasetSelectPanel() {
        JRadioButton rhoButton = new JRadioButton("rho");
        rhoButton.setMnemonic(KeyEvent.VK_R);
        rhoButton.setActionCommand("DATASET_RHO");
        rhoButton.setSelected(true);
        dataset = DATASET_RHO;
        rhoButton.addActionListener(new DatasetSelectorListener());

        JRadioButton fButton = new JRadioButton("|f|");
        fButton.setMnemonic(KeyEvent.VK_F);
        fButton.setActionCommand("DATASET_F");
        fButton.addActionListener(new DatasetSelectorListener());

        JRadioButton vButton = new JRadioButton("|v|");
        vButton.setMnemonic(KeyEvent.VK_V);
        vButton.setActionCommand("DATASET_V");
        vButton.addActionListener(new DatasetSelectorListener());

        ButtonGroup datasetSelectGroup = new ButtonGroup();
        datasetSelectGroup.add(rhoButton);
        datasetSelectGroup.add(fButton);
        datasetSelectGroup.add(vButton);

        JPanel datasetSelectPanel = new JPanel();
        datasetSelectPanel.setLayout(new GridLayout(1, 4));
        datasetSelectPanel.add(new JLabel("Dataset:"));
        datasetSelectPanel.add(rhoButton);
        datasetSelectPanel.add(fButton);
        datasetSelectPanel.add(vButton);
        return datasetSelectPanel;
    }

    private JPanel initColorMapSelectPanel() {
        // Initialize colormap selection
        JRadioButton rainbowButton = new JRadioButton("Rainbow");
        rainbowButton.setMnemonic(KeyEvent.VK_R);
        rainbowButton.setActionCommand("COLORMAP_RAINBOW");
        rainbowButton.addActionListener(new ColorSelectorListener());

        JRadioButton grayscaleButton = new JRadioButton("Grayscale");
        grayscaleButton.setMnemonic(KeyEvent.VK_G);
        grayscaleButton.setActionCommand("COLORMAP_GRAYSCALE");
        grayscaleButton.addActionListener(new ColorSelectorListener());

        JRadioButton definedButton   = new JRadioButton("Defined");
        definedButton.setMnemonic(KeyEvent.VK_D);
        definedButton.setActionCommand("COLORMAP_DEFINED");
        definedButton.addActionListener(new ColorSelectorListener());

        JRadioButton customButton   = new JRadioButton("Custom");
        customButton.setMnemonic(KeyEvent.VK_C);
        customButton.setActionCommand("COLORMAP_CUSTOM");
        customButton.addActionListener(new ColorSelectorListener());

        ButtonGroup colorMapSelectGroup = new ButtonGroup();
        colorMapSelectGroup.add(rainbowButton);
        colorMapSelectGroup.add(grayscaleButton);
        colorMapSelectGroup.add(definedButton);
        colorMapSelectGroup.add(customButton);

        JPanel colorMapSelectPanel = new JPanel();
        colorMapSelectPanel.setLayout(new GridLayout(4,1));
        colorMapSelectPanel.setBorder(new TitledBorder("Colormaps"));
        colorMapSelectPanel.add(rainbowButton);
        colorMapSelectPanel.add(grayscaleButton);
        colorMapSelectPanel.add(definedButton);
        colorMapSelectPanel.add(customButton);


        JPanel onoffpanel = new JPanel();
        onoffpanel.setBorder(new TitledBorder("Simulation"));

        JRadioButton simOnButton = new JRadioButton("On");
        simOnButton.setMnemonic(KeyEvent.VK_N);
        simOnButton.setActionCommand("SIMULATION_ON");
        simOnButton.addActionListener(new ColorSelectorListener());

        JRadioButton simOffButton = new JRadioButton("Off");
        simOffButton.setMnemonic(KeyEvent.VK_F);
        simOffButton.setActionCommand("SIMULATION_OFF");
        simOffButton.addActionListener(new ColorSelectorListener());

        ButtonGroup simOnOffGroup = new ButtonGroup();
        simOnOffGroup.add(simOnButton);
        simOnOffGroup.add(simOffButton);

        onoffpanel.add(simOnButton);
        onoffpanel.add(simOffButton);

        customButton.setSelected(scalar_col == COLOR_CUSTOM);
        simOnButton.setSelected(!frozen);
        simOffButton.setSelected(frozen);
        rainbowButton.setSelected(scalar_col == COLOR_RAINBOW);

        JPanel blaat = new JPanel();
        blaat.add(colorMapSelectPanel);
        blaat.add(onoffpanel);
        return blaat;
    }

    JLabel colorCountLabel = new JLabel("Limit colors to 255");
    private JPanel initScalingSelectPanel() {
        JCheckBox clampButton = new JCheckBox("clamping");
        clampButton.setMnemonic(KeyEvent.VK_C);
        clampButton.setActionCommand("SCALE_CLAMP");
        clampButton.addActionListener(new ScaleSelectListener());

        JCheckBox scaleButton = new JCheckBox("scaling");
        scaleButton.setMnemonic(KeyEvent.VK_S);
        scaleButton.setActionCommand("SCALE_SCALE");
        scaleButton.addActionListener(new ScaleSelectListener());
        scaleButton.setSelected(true);

        JSlider colorCountSlider = new JSlider(JSlider.HORIZONTAL, 0, 255, 255);
        colorCountSlider.addChangeListener(new CustomColorPanelHandler());

        JPanel clampSelectPanel = new JPanel();
        clampSelectPanel.setLayout(new GridLayout(3,2));
        clampSelectPanel.add(clampButton);
        clampSelectPanel.add(new JLabel());
        clampSelectPanel.add(new JLabel("Low: ", SwingConstants.RIGHT));
        clampSelectPanel.add(new JSpinner(new MinClampSelectSpinnerModel(0.0d)));
        clampSelectPanel.add(new JLabel("High: ", SwingConstants.RIGHT));
        clampSelectPanel.add(new JSpinner(new MaxClampSelectSpinnerModel(1.0d)));

        JPanel scaleSelectPanel = new JPanel();
        scaleSelectPanel.setLayout(new BoxLayout(scaleSelectPanel, BoxLayout.Y_AXIS));
        scaleButton.setAlignmentX(Component.LEFT_ALIGNMENT);
        clampSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        colorCountSlider.setAlignmentX(Component.LEFT_ALIGNMENT);
        scaleSelectPanel.setBorder(new TitledBorder("Scaling:"));
        scaleSelectPanel.add(scaleButton);
        scaleSelectPanel.add(clampSelectPanel);
        scaleSelectPanel.add(colorCountLabel);
        scaleSelectPanel.add(colorCountSlider);
        scaleSelectPanel.setAlignmentX(Component.RIGHT_ALIGNMENT);
        return scaleSelectPanel;


    }

    class MinClampSelectSpinnerModel implements SpinnerModel {
        private Vector listeners = new Vector();

        private double value = 0.0d;

        public MinClampSelectSpinnerModel(double value) {
            this.value = value;
            minClamp = value;
        }

        public void removeChangeListener(ChangeListener listener) {
            listeners.remove(listener);
        }

        public void addChangeListener(ChangeListener listener) {
            listeners.add(listener);
        }

        public Object getPreviousValue() {
            return new Double(value - 0.1d);
        }

        public Object getNextValue() {
            return new Double(value + 0.1d);
        }

        public void setValue(Object value) {
            this.value = ((Double)value).doubleValue();

            for(int i = 0; i < listeners.size(); ++i) {
                ChangeListener l = (ChangeListener)listeners.get(i);
                l.stateChanged(new ChangeEvent(this));
            }

            minClamp = this.value;
        }

        public Object getValue() {
            return new Double(value);
        }
    }

    class MaxClampSelectSpinnerModel implements SpinnerModel {
        private Vector listeners = new Vector();

        private double value = 0.0d;

        public MaxClampSelectSpinnerModel(double value) {
            this.value = value;
            maxClamp = value;
        }

        public void removeChangeListener(ChangeListener listener) {
            listeners.remove(listener);
        }

        public void addChangeListener(ChangeListener listener) {
            listeners.add(listener);
        }

        public Object getPreviousValue() {
            return new Double(value - 0.1d);
        }

        public Object getNextValue() {
            return new Double(value + 0.1d);
        }

        public void setValue(Object value) {
            this.value = ((Double)value).doubleValue();

            for(int i = 0; i < listeners.size(); ++i) {
                ChangeListener l = (ChangeListener)listeners.get(i);
                l.stateChanged(new ChangeEvent(this));
            }

            maxClamp = this.value;
        }

        public Object getValue() {
            return new Double(value);
        }
    }

    class ScaleSelectListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
            if (e.getActionCommand().equals("SCALE_CLAMP")) {
                scaleMode ^= SCALE_CLAMP;
            }
            else if (e.getActionCommand().equals("SCALE_SCALE")) {
                scaleMode ^= SCALE_SCALE;
            }
        }
    }

    private JPanel initSmokeSelectPanel()
    {
        JCheckBox smokeButton = new JCheckBox("Smoke");
        smokeButton.setMnemonic(KeyEvent.VK_S);
        smokeButton.setActionCommand("SMOKE_TOGGLE");
        smokeButton.addActionListener(new SmokeSelectorListener());

        JCheckBox vectorButton = new JCheckBox("Vectors");
        vectorButton.setMnemonic(KeyEvent.VK_V);
        vectorButton.setActionCommand("VECTOR_TOGGLE");
        vectorButton.setSelected(true);
        vectorButton.addActionListener(new SmokeSelectorListener());

        JPanel smokeSelectPanel = new JPanel();
        smokeSelectPanel.setLayout(new GridLayout(3,1));
        smokeSelectPanel.add(new JLabel("Enable smoke and vectors:"));
        smokeSelectPanel.add(smokeButton);
        smokeSelectPanel.add(vectorButton);

        smokeButton.setSelected(true);
        vectorButton.setSelected(false);

        return smokeSelectPanel;
    }


		protected ColorSelector colorselector;
		ColorTable colortable;
		private JPanel initCustomColorPanel(JFrame frame) {
			colorselector = new ColorSelector(frame);
			colorselector.addChangeListener(new CustomColorPanelHandler());

			DefaultTableModel tablemodelcolors= new DefaultTableModel(1,1);

			tablemodelcolors.addTableModelListener(new ColorTableTableModelListener());
			colortable = new ColorTable(tablemodelcolors) {
				private static final long serialVersionUID = 1L; // prevent warning
				public Class getColumnClass(int column) { //enable JTable to use different renderers, eg Checkbox for Boolean
					return getValueAt(0, column).getClass();
				}
			};
			colortable.addMouseListener(new ColorTableRightClick(new CustomColorPanelHandler()));

			colortable.getColumnModel().getColumn(colortable.convertColumnIndexToView(0)).setHeaderValue("Colors");
			colortable.doLayout();

			TableColumnModel cm = colortable.getColumnModel();
			ColorTableRenderer r = new ColorTableRenderer();
			TableColumn c = cm.getColumn(0);
			c.setCellRenderer(r);

			tablemodelcolors.setRowCount(0);
			Object[] o = new Object[1];
			o[0] = new Color(0,0,255);
			tablemodelcolors.addRow(o);
			o[0] = new Color(0,255,0);
			tablemodelcolors.addRow(o);
			o[0] = new Color(255,255,0);
			tablemodelcolors.addRow(o);
			o[0] = new Color(255,0,0);
			tablemodelcolors.addRow(o);


			JScrollPane scroll  = new JScrollPane(colortable);

			JPanel panel = new JPanel();
			panel.setBorder(new TitledBorder("Current gradient"));
			BoxLayout panellayout = new BoxLayout(panel, BoxLayout.X_AXIS);
			panel.setLayout(panellayout);

			panel.add(scroll);

			return panel;
		}

		class CustomColorPanelHandler implements ActionListener, ChangeListener {

			public void actionPerformed(ActionEvent e) {
				if (e.getActionCommand().equals("COLORTABLE_PICK_COLOR")) {
					colorselector.setVisible(true);
				}
				else if (e.getActionCommand().equals("COLORTABLE_ADD_NEW_COLOR")) {
					int row = colortable.getSelectedRow();
					Object[] o=new Object[1];
					o[0]=new Color((int)(Math.random()*256), (int)(Math.random()*256), (int)(Math.random()*256));
					((DefaultTableModel)(colortable.getModel())).insertRow(row, o);
					colortable.changeSelection(row, 0, false, false);
				}
				else if (e.getActionCommand().equals("COLORTABLE_REMOVE_COLOR")) {
					int row = colortable.getSelectedRow();
					((DefaultTableModel)(colortable.getModel())).removeRow(row);
				}
			}
			public void stateChanged(ChangeEvent e) {
				if (e.getSource().getClass().getName().equals("java.awt.Color")) {
					int row = colortable.getSelectedRow();
					int column = colortable.getSelectedColumn();
					if(row>=0 && column>=0) {
						colortable.setValueAt(e.getSource(), row, column);
					}
				}
                                else if (e.getSource().getClass().getName().equals("javax.swing.JSlider")) {
                                    int value = ((JSlider)e.getSource()).getValue();
                                    colorCountLabel.setText("Limit colors to " + value);
                                    ncolors = value;
                                }
			}
		}

    private JComponent initOptionPanel(JFrame frame)
    {
        JTabbedPane tabPane = new JTabbedPane();

        // Initialize option panel
        JPanel optionPanel = new JPanel();
        optionPanel.setLayout(new BoxLayout(optionPanel, BoxLayout.Y_AXIS));
        optionPanel.add(initDatasetSelectPanel());
        optionPanel.add(initColorMapSelectPanel());
        optionPanel.add(initScalingSelectPanel());
        optionPanel.add(initSmokeSelectPanel());
				optionPanel.add(initCustomColorPanel(frame));

        tabPane.addTab("Colors", optionPanel);

//         frame.add(optionPanel, BorderLayout.EAST);
	return tabPane;
    }

	GLCanvas panel;

    public Smoke() {
        init_simulation(DIM);	//initialize the simulation data structures

        // initialize GUI
        JFrame frame = new JFrame("Real-time smoke simulation and visualization");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        GLCapabilities caps = new GLCapabilities();
        caps.setDoubleBuffered(true);

        // initialize opengl Panel
        panel = new GLCanvas(caps);
        panel.addGLEventListener(new MyGLEventListener());
        panel.addMouseMotionListener(new MouseListener());
        panel.addKeyListener(new MyKeyListener());
        panel.setFocusable(true);
				//panel.setPreferredSize(new Dimension(500,0));

        // add panel to window
        frame.setLayout(new BorderLayout());
				JComponent anotherpanel = initOptionPanel(frame);
				anotherpanel.setPreferredSize(new Dimension(256,0));
        frame.add(anotherpanel, BorderLayout.EAST);
        frame.add(panel,BorderLayout.CENTER);
        frame.setSize(1000,750);
				frame.doLayout();

				Point     p = new Point();
				Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
				frame.setLocation((int) (d.getWidth()/2) - (frame.getWidth()/2),(int) (d.getHeight()/2)-(frame.getHeight()/2));

        // show window
        frame.setVisible(true);
    }

    class MyKeyListener extends KeyAdapter {
        public void keyTyped(KeyEvent e) {
            char key = e.getKeyChar();
            Smoke.this.keyboard(key, 0, 0);
        }
    }

    class MouseListener extends MouseMotionAdapter {
        public void mouseDragged(MouseEvent e) {
            drag(e.getX(), e.getY());
            //Flow.this.panel.repaint();
        }
    }

    class MyGLEventListener implements GLEventListener {
        public void init(GLAutoDrawable drawable) {  }

        public void display(GLAutoDrawable drawable) {
            Smoke.this.do_one_simulation_step();
            GL gl = drawable.getGL();
            Smoke.this.display(gl);
            gl.glFlush();
        }

        public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {
            Smoke.this.reshape(drawable.getGL(),width,height);
        }

        public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged) { }
    }
}
