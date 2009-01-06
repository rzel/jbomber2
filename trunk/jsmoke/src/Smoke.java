import com.sun.opengl.util.texture.*;
import com.sun.jna.Native;
import com.sun.jna.Pointer;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.GradientPaint;
import java.awt.color.ColorSpace;
import java.util.Vector;
import java.util.*;
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
import javax.swing.DefaultBoundedRangeModel;
import java.awt.event.*;
import java.awt.*;
import javax.swing.table.*;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.colorchooser.AbstractColorChooserPanel;
import java.nio.*;
import com.sun.opengl.util.*;
import javax.imageio.*;
import java.awt.image.*;
import java.io.*;
import javax.swing.event.MouseInputListener;
/** Usage: Drag with the mouse to add smoke to the fluid. This will also move a "rotor" that disturbs
 *        the velocity field at the mouse location. Press the indicated keys to change options
 */
class Smoke {
	boolean updateisotexture = true;
	static RFFTWLibrary RFFTW = (RFFTWLibrary)Native.loadLibrary("rfftw", RFFTWLibrary.class);

//--- SIMULATION PARAMETERS ------------------------------------------------------------------------
	// final int DIM = 16;		//size of simulation grid
	int DIM = 76;		//size of simulation grid
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
	boolean draw_smoke     = true;  //draw the smoke or not
	boolean draw_vecs      = false;    //draw the vector field or not
	boolean draw_iso_lines = true;    //draw the iso lines or not
	static final int VECTOR_TYPE_HEDGEHOG = 0;
	static final int VECTOR_TYPE_ARROW    = VECTOR_TYPE_HEDGEHOG + 1;
	int vector_type = VECTOR_TYPE_ARROW;
	boolean frozen = false;         //toggles on/off the animation
	private static int TEXTURE_COLORMAP_SMOKE    = 0;
	private static int TEXTURE_COLORMAP_ISOLINES = 1;
	private static int TEXTURE_ARROW_1           = 2;
	private static int TEXTURE_ARROW_2           = 3;
	private static int TEXTURE_ARROW_3           = 4;
	private static int TEXTURE_COUNT             = 5;
	int[] textures        = new int[TEXTURE_COUNT];
	double[] texture_fill = new double[TEXTURE_COUNT];

	private static double mousex_prev, mousey_prev, mousex, mousey, rotationx, rotationy, scale;

	int[][] msquare_lookup = {
		{-1, -1, -1, -1, -1},   // 0: Case 0
		{ 0,  3, -1, -1, -1},   // 1: Case 8
		{ 0,  1, -1, -1, -1},   // 2: Case 4
		{ 1,  3, -1, -1, -1},   // 3: Case 12
		{ 1,  2, -1, -1, -1},   // 4: Case 2
		{ 0,  3,  1,  2, -1},   // 5: Case 10 <= ambigu
		{ 0,  2, -1, -1, -1},   // 6: Case 6
		{ 2,  3, -1, -1, -1},   // 7: Case 14
		{ 2,  3, -1, -1, -1},   // 8: Case 1
		{ 0,  2, -1, -1, -1},   // 9: Case 9
		{ 0,  1,  2,  3, -1},   // 10: Case 5 <= ambigu
		{ 1,  2, -1, -1, -1},   // 11: Case 13
		{ 1,  3, -1, -1, -1},   // 12: Case 3
		{ 0,  1, -1, -1, -1},   // 13: Case 11
		{ 0,  3, -1, -1, -1},   // 14: Case 7
		{-1, -1, -1, -1, -1},   // 15: Case 15
	};



	private ColormapSelectPanel smokeColormapSelectPanel;
	private VectorOptionSelectPanel vectorOptionSelectPanel;
        private IsoLineSelectPanel isoLineSelectPanel;
        private HeightplotSelectPanel heightplotSelectPanel;
	//------ SIMULATION CODE STARTS HERE -----------------------------------------------------------------

	/**init_simulation: Initialize simulation data structures as a function of the grid size 'n'.
	 *                 Although the simulation takes place on a 2D grid, we allocate all data structures as 1D arrays,
	 *                 for compatibility with the FFTW numerical library.
	 */
	void init_simulation(int n) {
		int dim     = n * 2 * (n / 2 + 1);        //Allocate data structures
		vx       = new double/*fftw_real*/[dim];
		vy       = new double/*fftw_real*/[dim];
		vx0      = new double/*fftw_real*/[dim];
		vy0      = new double/*fftw_real*/[dim];
		for (int i = 0; i < dim; i++)                      //Initialize data structures to 0
			{ vx[i] = vy[i] = vx0[i] = vy0[i] = 0.0; } // float or double??

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
	void FFT(int direction, double/*fftw_real*/[] vx) {
		if (direction == 1) RFFTW.rfftwnd_one_real_to_complex_in_place(plan_rc, (double[])vx);
		else             RFFTW.rfftwnd_one_complex_to_real_in_place(plan_cr, (double[])vx);
	}

	int clamp(double x) {
		return ((x) >= 0.0 ? ((int)(x)) : (-((int)(1 - (x)))));
	}

//solve: Solve (compute) one step of the fluid flow simulation
	void solve(int n, double/*fftw_real*/[] vx, double/*fftw_real*/[] vy, double/*fftw_real*/[] vx0, double/*fftw_real*/[] vy0, double/*fftw_real*/ visc, double dt) {
		double/*fftw_real*/ x, y, x0, y0, f, r, U[] = new double/*fftw_real*/[2], V[] = new double/*fftw_real*/[2], s, t;
		int i, j, i0, j0, i1, j1;

		for (i = 0;i < n*n;i++) {
			vx[i] += dt * vx0[i]; vx0[i] = vx[i]; vy[i] += dt * vy0[i]; vy0[i] = vy[i];
		}

		for (x = 0.5 / n, i = 0 ; i < n ; i++, x += 1.0 / n) {
			for (y = 0.5 / n, j = 0 ; j < n ; j++, y += 1.0 / n) {
				x0 = n * (x - dt * vx0[i+n*j]) - 0.5f;
				y0 = n * (y - dt * vy0[i+n*j]) - 0.5f;
				i0 = clamp(x0); s = x0 - i0;
				i0 = (n + (i0 % n)) % n;
				i1 = (i0 + 1) % n;
				j0 = clamp(y0); t = y0 - j0;
				j0 = (n + (j0 % n)) % n;
				j1 = (j0 + 1) % n;
				vx[i+n*j] = (1 - s) * ((1 - t) * vx0[i0+n*j0] + t * vx0[i0+n*j1]) + s * ((1 - t) * vx0[i1+n*j0] + t * vx0[i1+n*j1]);
				vy[i+n*j] = (1 - s) * ((1 - t) * vy0[i0+n*j0] + t * vy0[i0+n*j1]) + s * ((1 - t) * vy0[i1+n*j0] + t * vy0[i1+n*j1]);
			}
		}

		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				vx0[i+(n+2)*j] = vx[i+n*j]; vy0[i+(n+2)*j] = vy[i+n*j];
			}
		}

		FFT(1, vx0);
		FFT(1, vy0);

		for (i = 0;i <= n;i += 2) {
			x = 0.5 * i;
			for (j = 0;j < n;j++) {
				y = j <= n / 2 ? (double)j : (double)(j - n);
				r = x * x + y * y;
				if (r == 0.0f) continue;
				f = Math.exp(-r * dt * visc);
				U[0] = vx0[i  +(n+2)*j]; V[0] = vy0[i  +(n+2)*j];
				U[1] = vx0[i+1+(n+2)*j]; V[1] = vy0[i+1+(n+2)*j];

				vx0[i  +(n+2)*j] = f * ((1 - x * x / r) * U[0]     - x * y / r * V[0]);
				vx0[i+1+(n+2)*j] = f * ((1 - x * x / r) * U[1]     - x * y / r * V[1]);
				vy0[i+ (n+2)*j] = f * (-y * x / r * U[0] + (1 - y * y / r) * V[0]);
				vy0[i+1+(n+2)*j] = f * (-y * x / r * U[1] + (1 - y * y / r) * V[1]);
			}
		}

		FFT(-1, vx0);
		FFT(-1, vy0);

		f = 1.0 / (n * n);
		for (i = 0;i < n;i++) {
			for (j = 0;j < n;j++) {
				vx[i+n*j] = f * vx0[i+(n+2)*j]; vy[i+n*j] = f * vy0[i+(n+2)*j];
			}
		}
	}


// diffuse_matter: This function diffuses matter that has been placed in the velocity field. It's almost identical to the
// velocity diffusion step in the function above. The input matter densities are in rho0 and the result is written into rho.
	void diffuse_matter(int n, double/*fftw_real*/[] vx, double/*fftw_real*/[] vy, double/*fftw_real*/[] rho, double/*fftw_real*/[] rho0, double dt) {
		double/*fftw_real*/ x, y, x0, y0, s, t;
		int i, j, i0, j0, i1, j1;

		for (x = 0.5 / n, i = 0 ; i < n ; i++, x += 1.0 / n) {
			for (y = 0.5 / n, j = 0 ; j < n ; j++, y += 1.0 / n) {
				x0 = n * (x - dt * vx[i+n*j]) - 0.5;
				y0 = n * (y - dt * vy[i+n*j]) - 0.5;
				i0 = clamp(x0);
				s = x0 - i0;
				i0 = (n + (i0 % n)) % n;
				i1 = (i0 + 1) % n;
				j0 = clamp(y0);
				t = y0 - j0;
				j0 = (n + (j0 % n)) % n;
				j1 = (j0 + 1) % n;
				rho[i+n*j] = (1 - s) * ((1 - t) * rho0[i0+n*j0] + t * rho0[i0+n*j1]) + s * ((1 - t) * rho0[i1+n*j0] + t * rho0[i1+n*j1]);
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


	//direction_to_color: Set the current color by mapping a direction vector (x,y), using
	//                    the color mapping method 'method'. If method==1, map the vector direction
	//                    using a rainbow colormap. If method==0, simply use the white color
	void direction_to_color(GL gl, float x, float y, int method) {
		float r, g, b;
		if (method == 1) {
			float f = (float)(Math.atan2(y, x) / Math.PI + 1);
			r = f;
			if (r > 1) r = 2 - r;
			g = f + .66667f;
			if (g > 2) g -= 2;
			if (g > 1) g = 2 - g;
			b = f + 2 * .66667f;
			if (b > 2) b -= 2;
			if (b > 1) b = 2 - b;
		} else { // method==0
			r = g = b = 1;
		}
		gl.glColor3f(r, g, b);
	}

	private int[] getXYFromIdx(int idx) {
			int y = idx / DIM;
			int x = idx - (DIM * y);
			return new int[] {x, y};
	}

	private int getIdxFromXY(int x, int y) {
			x = (x + DIM) % DIM;
			y = (y + DIM) % DIM;
			return (x + (y * DIM));
	}

	private float getDatasetColor(int idx, ColormapSelectPanel panel) {
		float dataset_value = 0.0f;
		switch (panel.getDataset()) {
			case ColormapSelectPanel.DATASET_RHO:
				dataset_value = (float)rho[idx];
				break;
			case ColormapSelectPanel.DATASET_F:
				dataset_value = (float)Math.sqrt(fx[idx] * fx[idx] + fy[idx] * fy[idx]) * DIM;
				break;
			case ColormapSelectPanel.DATASET_V:
				dataset_value = (float)Math.sqrt(vx[idx] * vx[idx] + vy[idx] * vy[idx]) * DIM;
				break;
			case ColormapSelectPanel.DATASET_F_DIV: {
				int[] pos = getXYFromIdx(idx);
				double div = 0.5 * ((fx[getIdxFromXY(pos[0] - 1, pos[1])] - fx[getIdxFromXY(pos[0] + 1, pos[1])])
				                 + (fy[getIdxFromXY(pos[0], pos[1] - 1)] - fy[getIdxFromXY(pos[0], pos[1] + 1)]));
				dataset_value = (float)div;
				}; break;
			case ColormapSelectPanel.DATASET_V_DIV: {
				int[] pos = getXYFromIdx(idx);
				double div = 0.5 * ((vx[getIdxFromXY(pos[0] - 1, pos[1])] - vx[getIdxFromXY(pos[0] + 1, pos[1])]) +
				                 + (vy[getIdxFromXY(pos[0], pos[1] - 1)] - vy[getIdxFromXY(pos[0], pos[1] + 1)]));
				dataset_value = (float)div;
			} break;
		}

		// Clamp
		if ((panel.getScalemode() & panel.SCALE_CLAMP) == panel.SCALE_CLAMP) {
			double minClamp = panel.getMinClamp();
			double maxClamp = panel.getMaxClamp();
			dataset_value = (float)(dataset_value < minClamp ? minClamp : dataset_value > maxClamp ? maxClamp : dataset_value); // Clamp!
		}

		panel.update_maxdataset_value(dataset_value);
		panel.update_mindataset_value(dataset_value);

		if ((panel.getScalemode() & panel.SCALE_SCALE) == panel.SCALE_SCALE) {
			dataset_value = (dataset_value - panel.get_mindataset_value()) / (panel.get_maxdataset_value() - panel.get_mindataset_value());
		}

		return dataset_value;
	}

	private double[] sampleDataset(int x, int y, ColormapSelectPanel panel, double gridx, double gridy) {

		double[] vfx, vfy;

		if (vectorOptionSelectPanel.getVectorField() == VectorOptionSelectPanel.VECTOR_FIELD_FORCE) {
		    vfx = fx;
		    vfy = fy;
		}
		else {
		    vfx = vx;
		    vfy = vy;
		}

		// x and y here represent the (x,y) id of the grid-rectangle we're sampling
// 		double gridx = vectorOptionSelectPanel.getVectorGridX();
// 		double gridy = vectorOptionSelectPanel.getVectorGridY();
		double cells_per_sample_x = DIM / gridx;
		double cells_per_sample_y = DIM / gridy;
		double dim = (double)DIM;
		double x_sample_centre_pos = (dim/(gridx+1))*(x+1);//((x+1) / (gridx+1)) * DIM;
		double y_sample_centre_pos = (dim/(gridy+1))*(y+1);//((y+1) / (gridy+1)) * DIM;

		//
		// 2           4           6           8         10           12
		//   2.66 3.33   4.66 5.33   6.66 7.33 ..


		double weight_sx  = 0.5 - Math.sqrt(((x_sample_centre_pos - (int)x_sample_centre_pos) - 0.5)
		                                    * ((x_sample_centre_pos - (int)x_sample_centre_pos) - 0.5));
		double weight_sy  = 0.5 - Math.sqrt(((y_sample_centre_pos - (int)y_sample_centre_pos) - 0.5)
		                                    * ((y_sample_centre_pos - (int)y_sample_centre_pos) - 0.5));
		double weight_nnx = 1.0 - weight_sx;
		double weight_nny = 1.0 - weight_sy;

		int nearest_neighbour_x = (int)(/*1.5 **/ ((int)((x_sample_centre_pos - (int)x_sample_centre_pos) + 0.5)) - 1.0);
		int nearest_neighbour_y = (int)(/*1.5 **/ ((int)((y_sample_centre_pos - (int)y_sample_centre_pos) + 0.5)) - 1.0);

		double avgs = 0.0;
		double avgx = 0.0;
		double avgy = 0.0;
		int samples = 0;
		double[] result = new double[5];
		for (double sample_x = -(cells_per_sample_x / 2.0); sample_x < cells_per_sample_x / 2.0; sample_x += 1.0) {
			for (double sample_y = -(cells_per_sample_y / 2.0); sample_y < cells_per_sample_y / 2.0; sample_y += 1.0) {
				double px, py;
				px = x_sample_centre_pos + sample_x;
				py = y_sample_centre_pos + sample_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell1s  = getDatasetColor((int)(px)+(int)(py)*DIM, panel);
				double cell1vx = vfx[(int)(px)+(int)(py)*DIM];
				double cell1vy = vfy[(int)(px)+(int)(py)*DIM];
				px = x_sample_centre_pos + sample_x + nearest_neighbour_x;
				py = y_sample_centre_pos + sample_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell2s  = getDatasetColor((int)(px)+(int)(py)*DIM, panel);
				double cell2vx = vfx[(int)(px)+(int)(py)*DIM];
				double cell2vy = vfy[(int)(px)+(int)(py)*DIM];
				px = x_sample_centre_pos + sample_x;
				py = y_sample_centre_pos + sample_y + nearest_neighbour_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell3s  = getDatasetColor((int)(px)+(int)(py)*DIM, panel);
				double cell3vx = vfx[(int)(px)+(int)(py)*DIM];
				double cell3vy = vfy[(int)(px)+(int)(py)*DIM];
				px = x_sample_centre_pos + sample_x + nearest_neighbour_x;
				py = y_sample_centre_pos + sample_y + nearest_neighbour_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell4s  = getDatasetColor((int)(px)+(int)(py)*DIM, panel);
				double cell4vx = vfx[(int)(px)+(int)(py)*DIM];
				double cell4vy = vfy[(int)(px)+(int)(py)*DIM];
				double cavg = (weight_sx * cell1s + weight_nnx * cell2s) * weight_sy + (weight_sx * cell3s + weight_nnx * cell4s) * weight_nny;
				avgs += cavg;
				avgx += (weight_sx * cell1vx + weight_nnx * cell2vx) * weight_sy + (weight_sx * cell3vx + weight_nnx * cell4vx) * weight_nny;
				avgy += (weight_sy * cell1vy + weight_nny * cell2vy) * weight_sy + (weight_sy * cell3vy + weight_nny * cell4vy) * weight_nny;
				++samples;
				result[4] = Math.max(result[4], cavg);
			}
		}

		double len = Math.sqrt(avgx*avgx + avgy*avgy);
		result[0] = avgx;
		result[1] = avgy;
		result[2] = len;
		result[3] = avgs / ((double)samples);
		return result;
	}

	double[] calculateNormal(double[] a, double[] b) {
		double[] n = {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
		return n;
	}

	void setNormal(GL gl, int idx, ColormapSelectPanel panel, double zscale) {
		int[] xy = getXYFromIdx(idx);
		int x = xy[0];
		int y = xy[1];
		double[] v = {
			getDatasetColor(getIdxFromXY(x  ,y  ), panel) * zscale, // 0
			getDatasetColor(getIdxFromXY(x-1,y  ), panel) * zscale, // 1
			getDatasetColor(getIdxFromXY(x  ,y-1), panel) * zscale, // 2
			getDatasetColor(getIdxFromXY(x+1,y  ), panel) * zscale, // 3
			getDatasetColor(getIdxFromXY(x  ,y+1), panel) * zscale, // 4
		};
		double[][] n = {
			calculateNormal(new double[]{ 0, 1, v[4] - v[0]}, new double[]{-1, 0, v[1] - v[0]}),
			calculateNormal(new double[]{-1, 0, v[1] - v[0]}, new double[]{ 0,-1, v[2] - v[0]}),
			calculateNormal(new double[]{ 0,-1, v[2] - v[0]}, new double[]{ 1, 0, v[3] - v[0]}),
			calculateNormal(new double[]{ 1, 0, v[3] - v[0]}, new double[]{ 0, 1, v[4] - v[0]}),
		};
		gl.glNormal3d(
			(n[0][0] + n[1][0] + n[2][0] + n[3][0]) / 4.0,
			(n[0][1] + n[1][1] + n[2][1] + n[3][1]) / 4.0,
			(n[0][2] + n[1][2] + n[2][2] + n[3][2]) / 4.0
		);
	}

	//visualize: This is the main visualization function
	void visualize(GL gl) {
		if(textures[TEXTURE_COUNT-1]==-1) return;
		int        i, j, idx; double px, py;
		double z, zscale = Math.max(winWidth, winHeight) * 0.1;
		double/*fftw_real*/  wn = winWidth / (double)(DIM + 1);   // Grid cell width
		double/*fftw_real*/  hn = winHeight / (double)(DIM + 1);  // Grid cell heigh


		if (heightplotSelectPanel.isShadingEnabled()) {
				gl.glEnable(gl.GL_LIGHTING);
		}
		boolean SomethingChangedInTheLightingModel = heightplotSelectPanel.getSomethingChangedInTheLightingModel();

		if(SomethingChangedInTheLightingModel) {
			heightplotSelectPanel.setSomethingChangedInTheLightingModel(false);
			float[] fLightAmbient  = heightplotSelectPanel.getLightAmbientColor();
			float[] fLightDiffuse  = heightplotSelectPanel.getLightDiffuseColor();
			float[] fLightSpecular = { 1.0f, 1.0f, 1.0f, 1.0f };
			float[] fLightPosition = heightplotSelectPanel.getLightPosition();

			FloatBuffer LightAmbient  = FloatBuffer.wrap(fLightAmbient);
			FloatBuffer LightDiffuse  = FloatBuffer.wrap(fLightDiffuse);
			FloatBuffer LightSpecular  = FloatBuffer.wrap(fLightSpecular);
			FloatBuffer LightPosition = FloatBuffer.wrap(fLightPosition);

			if (heightplotSelectPanel.getFlatShading()) {
					gl.glShadeModel(gl.GL_FLAT); // needs param
			}
			else {
					gl.glShadeModel(gl.GL_SMOOTH); // needs param
			}

			gl.glEnable(gl.GL_LIGHT0);
			gl.glLightfv(gl.GL_LIGHT0, gl.GL_AMBIENT,  LightAmbient);
			gl.glLightfv(gl.GL_LIGHT0, gl.GL_DIFFUSE,  LightDiffuse);
			gl.glLightfv(gl.GL_LIGHT0, gl.GL_SPECULAR,  LightSpecular);
			gl.glLightfv(gl.GL_LIGHT0, gl.GL_POSITION, LightPosition);
		}

		gl.glTranslated(winWidth/2.0, 0.0, 0.0);
		gl.glTranslated(0.0, winHeight/2.0, 0.0);
		gl.glRotated(rotationy, 1.0, 0.0, 0.0);
		gl.glRotated(rotationx, 0.0, 0.0, 1.0);
		if(scale > 0)
			gl.glScaled(scale, scale, scale);
		if(scale < 0)
			gl.glScaled(1.0/(-scale), 1.0/(-scale), 1.0/(-scale));
		gl.glTranslated(-winWidth/2.0, 0.0, 0.0);
		gl.glTranslated(0.0, -winHeight/2.0, 0.0);

		gl.glDisable(gl.GL_TEXTURE_2D);
		gl.glEnable(gl.GL_TEXTURE_1D);
		gl.glDisable(gl.GL_BLEND);
		gl.glColor4d(1.0, 1.0, 1.0, 1.0);

		if (draw_smoke) {
			gl.glBindTexture(gl.GL_TEXTURE_1D, textures[TEXTURE_COLORMAP_SMOKE]);
			gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
			for (j = 0; j < DIM - 1; j++) {		//draw smoke
				gl.glBegin(GL.GL_TRIANGLE_STRIP);

				i = 0;
				px = wn + (float)i * wn;
				py = hn + (float)j * hn;
				idx = (j * DIM) + i;
				gl.glTexCoord1d(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill[TEXTURE_COLORMAP_SMOKE]);
				z = getDatasetColor(idx, heightplotSelectPanel);
				if(heightplotSelectPanel.isShadingEnabled()) setNormal(gl, idx, heightplotSelectPanel, zscale);
				gl.glVertex3d(px, py, z * zscale);

				for (i = 0; i < DIM - 1; i++) {
					px = wn + i * wn;
					py = hn + (j + 1) * hn;
					idx = ((j + 1) * DIM) + i;

					gl.glTexCoord1d(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill[TEXTURE_COLORMAP_SMOKE]);
					z = getDatasetColor(idx, heightplotSelectPanel);
					if(heightplotSelectPanel.isShadingEnabled()) setNormal(gl, idx, heightplotSelectPanel,zscale);
					gl.glVertex3d(px, py, z * zscale);
					px = wn + (i + 1) * wn;
					py = hn + j * hn;
					idx = (j * DIM) + (i + 1);
					gl.glTexCoord1d(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill[TEXTURE_COLORMAP_SMOKE]);
					z = getDatasetColor(idx, heightplotSelectPanel);
					if(heightplotSelectPanel.isShadingEnabled()) setNormal(gl, idx, heightplotSelectPanel,zscale);
					gl.glVertex3d(px, py, z * zscale);
				}

				px = wn + (float)(DIM - 1) * wn;
				py = hn + (float)(j + 1) * hn;
				idx = ((j + 1) * DIM) + (DIM - 1);
				gl.glTexCoord1d(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill[TEXTURE_COLORMAP_SMOKE]);
				z = getDatasetColor(idx, heightplotSelectPanel);
				if(heightplotSelectPanel.isShadingEnabled()) setNormal(gl, idx, heightplotSelectPanel,zscale);
				gl.glVertex3d(px, py, z * zscale);
				gl.glEnd();
			}
		}

		gl.glDisable(gl.GL_LIGHTING);

		if(draw_iso_lines) {
			gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE);
			gl.glBindTexture(gl.GL_TEXTURE_1D, textures[TEXTURE_COLORMAP_ISOLINES]);

			gl.glBegin(gl.GL_LINES);
			int    iso_n_lines    = isoLineSelectPanel.getIsoLineCount();
			double iso_low_value  = Math.min(isoLineSelectPanel.getMinIsoValue(), isoLineSelectPanel.getMaxIsoValue());
			double iso_high_value = Math.max(isoLineSelectPanel.getMinIsoValue(), isoLineSelectPanel.getMaxIsoValue());
			double iso_clip_hack  = 128.0 / zscale; // Magix

			for(int y = 0; y < DIM-1; ++y) {
				for(int x = 0; x < DIM-1; ++x) {
					int idx_top_lft = ((x     + DIM) % DIM) + (((y     + DIM) % DIM) * DIM);
					int idx_top_rgt = ((x + 1 + DIM) % DIM) + (((y     + DIM) % DIM) * DIM);
					int idx_btm_lft = ((x     + DIM) % DIM) + (((y + 1 + DIM) % DIM) * DIM);
					int idx_btm_rgt = ((x + 1 + DIM) % DIM) + (((y + 1 + DIM) % DIM) * DIM);

					double top_lft_value = getDatasetColor(idx_top_lft, isoLineSelectPanel);
					double top_rgt_value = getDatasetColor(idx_top_rgt, isoLineSelectPanel);
					double btm_lft_value = getDatasetColor(idx_btm_lft, isoLineSelectPanel);
					double btm_rgt_value = getDatasetColor(idx_btm_rgt, isoLineSelectPanel);

					double htop_lft_value = getDatasetColor(idx_top_lft, heightplotSelectPanel);
					double htop_rgt_value = getDatasetColor(idx_top_rgt, heightplotSelectPanel);
					double hbtm_lft_value = getDatasetColor(idx_btm_lft, heightplotSelectPanel);
					double hbtm_rgt_value = getDatasetColor(idx_btm_rgt, heightplotSelectPanel);

					if(iso_high_value - iso_low_value < 0.001) iso_n_lines = 1; // optimalize ftw
					for(int n = 0; n < iso_n_lines; ++n) {
						double iso_value = iso_low_value + (((iso_high_value - iso_low_value) / (double)(iso_n_lines+1)) * (double)(n+1));
						gl.glTexCoord1d(iso_value);
						// Calculate lookup table index
						int lookup_index = 0;
						if ( top_lft_value > iso_value) lookup_index |= 1;
						if ( top_rgt_value > iso_value) lookup_index |= 2;
						if ( btm_rgt_value > iso_value) lookup_index |= 4;
						if ( btm_lft_value > iso_value) lookup_index |= 8;

						i = 0;
						while( msquare_lookup[lookup_index][i] != -1 ) {
							switch (msquare_lookup[lookup_index][i]) {
								case 0: { // Top edge
									double iso_position = ((iso_value - top_lft_value) / (top_rgt_value-top_lft_value));
									z = iso_position * htop_rgt_value + (1.0-iso_position)*htop_lft_value;
									iso_position *= wn;
									gl.glVertex3d(wn + x * wn + iso_position, hn + y * hn, z * zscale + iso_clip_hack);
								}; break;
								case 1: { // Right edge
									double iso_position = ((iso_value - top_rgt_value) / (btm_rgt_value-top_rgt_value));
									z = iso_position * hbtm_rgt_value + (1.0-iso_position)*htop_rgt_value;
									iso_position *= hn;
									gl.glVertex3d(wn + x * wn + wn, hn + y * hn + iso_position, z * zscale + iso_clip_hack);
								}; break;
								case 2: { // Bottom edge
									double iso_position = ((iso_value - btm_rgt_value) / (btm_lft_value-btm_rgt_value));
									z = iso_position * hbtm_lft_value + (1.0-iso_position)*hbtm_rgt_value;
									iso_position *= wn;
									gl.glVertex3d(wn + x * wn + (wn - iso_position), hn + y * hn + hn, z * zscale + iso_clip_hack);
								}; break;
								case 3: { // Left edge
									double iso_position = ((iso_value - btm_lft_value) / (top_lft_value-btm_lft_value));
									z = iso_position * htop_lft_value + (1.0-iso_position)*hbtm_lft_value;
									iso_position *= hn;
									gl.glVertex3d(wn + x * wn, hn + y * hn + (hn - iso_position), z * zscale + iso_clip_hack);
								}; break;
							}
							++i;
						}
					}
				}
			}
			gl.glEnd();
		}

		if (draw_vecs) {
			double[] vfx, vfy;
			
			if (vectorOptionSelectPanel.getVectorField() == VectorOptionSelectPanel.VECTOR_FIELD_FORCE) {
			    vfx = fx;
			    vfy = fy;
			}
			else {
			    vfx = vx;
			    vfy = vy;
			}
		    
			if (vector_type == VECTOR_TYPE_HEDGEHOG) { // Fixme: These are currently b0rken
				gl.glBegin(GL.GL_LINES);				//draw velocities
				for (i = 0; i < DIM; i++)
					for (j = 0; j < DIM; j++) {
						idx = (j * DIM) + i;
						//direction_to_color(gl, (float)(double)vx[idx],(float)(double)vy[idx],color_dir);
						//set_colormap(gl, getDatasetColor(idx, vectorOptionSelectPanel), vectorOptionSelectPanel);
						z = getDatasetColor(idx, heightplotSelectPanel);
						gl.glTexCoord1d(z * texture_fill[TEXTURE_COLORMAP_SMOKE]);//FIXME: Needs own texture
						gl.glVertex3d(wn + i * wn, hn + j * hn, z * zscale);
						gl.glVertex3d((wn + i * wn) + vec_scale * vfx[idx], (hn + j * hn) + vec_scale * vfy[idx], z * zscale);
					}
				gl.glEnd();
			} else if (vector_type == VECTOR_TYPE_ARROW) {
				gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
				gl.glDisable(gl.GL_TEXTURE_1D);
				gl.glEnable(gl.GL_TEXTURE_2D);
				gl.glBindTexture(gl.GL_TEXTURE_2D, textures[TEXTURE_ARROW_2]);
				gl.glEnable(gl.GL_BLEND);
				float vector_size  = 0.5f * vectorOptionSelectPanel.getVectorSize();
				float vector_scalefactor = vectorOptionSelectPanel.getVectorScaleFactor();
				int gridx = vectorOptionSelectPanel.getVectorGridX();
				int gridy = vectorOptionSelectPanel.getVectorGridY();
				float spacex = (float)((winWidth  - 2 * wn) / (gridx + 0.0f));
				float spacey = (float)((winHeight - 2 * hn) / (gridy + 0.0f));
				gl.glPushMatrix();
				gl.glTranslatef(winWidth / 2.0f, winHeight / 2.0f, 0);
				double maxveclen = vectorOptionSelectPanel.get_longest_vector();
				double maxveclenx = 0.99 * (winWidth / gridx);
				double maxvecleny = 0.99 * (winHeight / gridy);
				double vecscalefact = Math.min(maxveclenx, maxvecleny);
				double vec_clip_hack  = 128.0 / zscale; // Magix
				for (int x = 0 ; x < gridx ; ++x) {
					for (int y = 0 ; y < gridy ; ++y) {
						idx = (int)((x / (float)gridx) * DIM + DIM * (int)(DIM * (y / (float)gridy)));
						double[] result = sampleDataset(x, y, vectorOptionSelectPanel, vectorOptionSelectPanel.getVectorGridX(), vectorOptionSelectPanel.getVectorGridY());
						vectorOptionSelectPanel.update_longest_vector(result[2]);
						double inprod = result[1]/result[2];
						double xdir   = result[0]/Math.abs(result[0]);
						double rotation = (-xdir)*(Math.acos(inprod)/Math.PI*180)+180;
						double size = 0.5 * (vector_size * vector_scalefactor * Math.sqrt(result[2]));
						if((vectorOptionSelectPanel.getScalemode() & vectorOptionSelectPanel.SCALE_SCALE) != 0) {
							size = 0.5 * (result[2] / maxveclen) * (vecscalefact * vector_scalefactor) + vector_size;
						}
// 						size = 0.5 * vecscalefact;
						float[] color = vectorOptionSelectPanel.getGradientColor(result[3]);

						result = sampleDataset(x, y, heightplotSelectPanel, vectorOptionSelectPanel.getVectorGridX(), vectorOptionSelectPanel.getVectorGridY());

						gl.glColor3f(color[0], color[1], color[2]);
							gl.glPushMatrix();
						gl.glTranslatef((float)(spacex*0.5f + spacex * x - winWidth  * 0.5f + wn),
						                (float)(spacey*0.5f + spacey * y - winHeight * 0.5f + hn),
						                0.0f);
						gl.glRotatef((float)rotation, 0, 0, 1);
						gl.glBegin(GL.GL_QUADS); // Can not be moved outside of for-loop because of glTranslatef and glRotatef
						gl.glTexCoord2d(1.0, 0.0);
						gl.glVertex3d( - size, - size, result[4] * zscale + vec_clip_hack);

						gl.glTexCoord2d(0.0, 0.0);
						gl.glVertex3d( + size, - size, result[4] * zscale + vec_clip_hack);

						gl.glTexCoord2d(0.0, 1.0);
						gl.glVertex3d( + size, + size, result[4] * zscale + vec_clip_hack);

						gl.glTexCoord2d(1.0, 1.0);
						gl.glVertex3d( - size, + size, result[4] * zscale + vec_clip_hack);
						gl.glEnd();

						gl.glPopMatrix();
					}
				}
				gl.glPopMatrix();
			}
		}

		gl.glFlush(); // forces all opengl commands to complete. Blocking!!
		++avg_fc;
		long current = System.nanoTime();
		if (current - avg_begin > 500000000) {
			double fps = (((double)(avg_fc_prev * (avg_begin - avg_begin_prev)) + (double)(avg_fc * (current - avg_begin))) / ((double)(current - avg_begin_prev) * 0.5));
			frame.setTitle("Real-time smoke simulation and visualization - " + (int)(fps) + "." + (int)((fps - (int)fps)*10) + " fps");
			avg_begin_prev = avg_begin;
			avg_begin = current;
			avg_fc_prev = avg_fc;
			avg_fc = 0;
// TODO                                        colorOverviewSlider.setMaximum((int)maxvy_lastframe+1);
		}

		smokeColormapSelectPanel.reset_maxdataset_value();
		smokeColormapSelectPanel.reset_mindataset_value();
		vectorOptionSelectPanel.reset_maxdataset_value();
		vectorOptionSelectPanel.reset_mindataset_value();
		vectorOptionSelectPanel.reset_longest_vector();
		isoLineSelectPanel.reset_mindataset_value();
		isoLineSelectPanel.reset_maxdataset_value();
		heightplotSelectPanel.reset_mindataset_value();
		heightplotSelectPanel.reset_maxdataset_value();
	}
	static long avg_begin      = 0;
	static long avg_begin_prev = 0;
	static long avg_fc         = 0;
	static long avg_fc_prev    = 0;


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
		gl.glViewport(0, 0, w, h);
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glLoadIdentity();
// 		glu.gluOrtho2D(0.0, (double)w, 0.0, (double)h);
// 		glu.gluPerspective(92.0, (double)w/(double)h, -Math.max(w,h), Math.max(w,h));
		gl.glOrtho(0.0, (double)w, 0.0, (double)h, -Math.max(w,h), Math.max(w,h) );
		winWidth = w; winHeight = h;
	}

//------ INTERACTION CODE STARTS HERE -----------------------------------------------------------------

//keyboard: Handle key presses
	void keyboard(char key, int x, int y) {
		switch (key) {
			case 'c': color_dir = 1 - color_dir; break;
			case 'S': vec_scale *= 1.2; break;
			case 's': vec_scale *= 0.8; break;
		}
	}



// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
//       cursor movement. Also inject some new matter into the field at the mouse location.
	int lmx = 0, lmy = 0;				//remembers last mouse location
	void drag(int mx, int my) {
		int xi, yi, X, Y; double  dx, dy, len;

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
		System.out.println("S/s:   increase/decrease hedgehog scaling");
		System.out.println("c:     toggle direction coloring on/off");
		new Smoke();
	}

	class SmokeSelectorListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			if (e.getActionCommand().equals("SMOKE_TOGGLE")) {
				draw_smoke = !draw_smoke;
			} else if (e.getActionCommand().equals("VECTOR_TOGGLE")) {
				draw_vecs = !draw_vecs;
			} else if (e.getActionCommand().equals("ISO_LINE_TOGGLE")) {
				draw_iso_lines = !draw_iso_lines;
			}else {
				System.out.println("Smoke: " + e.getActionCommand());
			}
		}
	}

	private JPanel initSimOnOffPanel() {
		JPanel onoffpanel = new JPanel();
		onoffpanel.setBorder(new TitledBorder("Simulation"));
		onoffpanel.setLayout(new BoxLayout(onoffpanel, BoxLayout.X_AXIS));

		JRadioButton simOnButton = new JRadioButton("On");
		simOnButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				frozen = false;
                                panel.repaint(50);
			}
		});

		JRadioButton simOffButton = new JRadioButton("Off");
		simOffButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				frozen = true;
                                panel.repaint(50);
			}
		});
		ButtonGroup simOnOffGroup = new ButtonGroup();
		simOnOffGroup.add(simOnButton);
		simOnOffGroup.add(simOffButton);

		onoffpanel.add(simOnButton);
		onoffpanel.add(simOffButton);

		simOnButton.setSelected(!frozen);
		simOffButton.setSelected(frozen);

		onoffpanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		return onoffpanel;
	}

	private JPanel initSmokeSelectPanel() {
		JCheckBox smokeButton = new JCheckBox("Smoke");
		smokeButton.setMnemonic(KeyEvent.VK_S);
		smokeButton.setActionCommand("SMOKE_TOGGLE");
		smokeButton.addActionListener(new SmokeSelectorListener());

		JCheckBox vectorButton = new JCheckBox("Vectors");
		vectorButton.setMnemonic(KeyEvent.VK_V);
		vectorButton.setActionCommand("VECTOR_TOGGLE");
		vectorButton.setSelected(true);
		vectorButton.addActionListener(new SmokeSelectorListener());

		JCheckBox isoLineButton = new JCheckBox("Iso lines");
		isoLineButton.setMnemonic(KeyEvent.VK_I);
		isoLineButton.setActionCommand("ISO_LINE_TOGGLE");
		isoLineButton.setSelected(true);
		isoLineButton.addActionListener(new SmokeSelectorListener());

		JPanel smokeSelectPanel = new JPanel();
		smokeSelectPanel.setLayout(new GridLayout(3, 1));
		smokeSelectPanel.setBorder(new TitledBorder("Visualizations"));
		smokeSelectPanel.add(smokeButton);
		smokeSelectPanel.add(vectorButton);
		smokeSelectPanel.add(isoLineButton);

		smokeButton.setSelected(draw_smoke);
		vectorButton.setSelected(draw_vecs);
		isoLineButton.setSelected(draw_iso_lines);

		smokeSelectPanel.setLayout(new BoxLayout(smokeSelectPanel, BoxLayout.Y_AXIS));
		smokeSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		return smokeSelectPanel;
	}

	class SimParamListener implements ChangeListener {
		String type = "";
		SimParamListener(String s) {
			type = s;
		}
		public void stateChanged(ChangeEvent e) {
			if (type.equals("DIMENSION")) {
				DIM = (int)((Double)((JSpinner)e.getSource()).getValue()).doubleValue();
				init_simulation(DIM);
			} else if (type.equals("TIMESTEP")) {
				dt = ((Double)((JSpinner)e.getSource()).getValue()).doubleValue();
			} else if (type.equals("VISCOSITY")) {
				visc = ((Double)((JSpinner)e.getSource()).getValue()).doubleValue() / 1000.0d;
			}
		}
	}

	JPanel initSimParamsPanel() {
		JPanel simParamsPanel = new JPanel();
		simParamsPanel.setBorder(new TitledBorder("Simulation Parameters"));
		simParamsPanel.setLayout(new GridLayout(3, 2));
		JSpinner dimension = new JSpinner(new SpinnerNumberModel(76.0, 2.0, 150.0, 2.0));
		JSpinner timestep = new JSpinner(new SpinnerNumberModel(0.4d, 0.04d, 4.0d, 0.01d));
		JSpinner viscosity = new JSpinner(new SpinnerNumberModel(10.0d, 1.0d, 100.0d, 1.0d));
		dimension.addChangeListener(new SimParamListener("DIMENSION"));
		timestep.addChangeListener(new SimParamListener("TIMESTEP"));
		viscosity.addChangeListener(new SimParamListener("VISCOSITY"));

		JPanel test = new JPanel();
		test.add(new JLabel("Dimensions: "));
		test.add(dimension);

		test.add(new JLabel("Time-Step: "));
		test.add(timestep);

		test.add(new JLabel("Viscosity: "));
		test.add(viscosity);

		test.setLayout(new GridLayout(3, 2));
		test.setAlignmentX(Component.LEFT_ALIGNMENT);

		simParamsPanel.setLayout(new BoxLayout(simParamsPanel, BoxLayout.Y_AXIS));
		simParamsPanel.add(test);
		return simParamsPanel;
	}

	private JComponent initOptionPanel(JFrame frame) {
		JTabbedPane tabPane = new JTabbedPane();

		JPanel SimulationOptionPanel = new JPanel();
		SimulationOptionPanel.setLayout(new BoxLayout(SimulationOptionPanel, BoxLayout.Y_AXIS));
		SimulationOptionPanel.add(initSimOnOffPanel());
		SimulationOptionPanel.add(initSmokeSelectPanel());
		SimulationOptionPanel.add(initSimParamsPanel());
		tabPane.addTab("Simulation options", SimulationOptionPanel);

		// Initialize option panel
		smokeColormapSelectPanel = new ColormapSelectPanel(0, 1/*(int)maxvy_lastframe+1*/, 2047, ColormapSelectPanel.COLOR_CUSTOM, frame);
		tabPane.addTab("Smoke options", smokeColormapSelectPanel);

		vectorOptionSelectPanel = new VectorOptionSelectPanel(0, 1/*(int)maxvy_lastframe+1*/, 2047, ColormapSelectPanel.COLOR_CUSTOM, frame);
		tabPane.addTab("Vector options", vectorOptionSelectPanel);

                isoLineSelectPanel = new IsoLineSelectPanel(0, 1/*(int)maxvy_lastframe+1*/, 2047, ColormapSelectPanel.COLOR_CUSTOM, frame);
                tabPane.addTab("ISO line options", isoLineSelectPanel);

                heightplotSelectPanel = new HeightplotSelectPanel(0, 1/*(int)maxvy_lastframe+1*/, 2047, ColormapSelectPanel.COLOR_CUSTOM, frame);
                tabPane.addTab("Heightplot", heightplotSelectPanel);

		tabPane.setSelectedIndex(1);
		return tabPane;
	}

	GLCanvas panel;
	JFrame frame;
	public Smoke() {
		rotationx = 0;
		rotationy = 0;
		init_simulation(DIM);	//initialize the simulation data structures

		// initialize GUI
		frame = new JFrame("Real-time smoke simulation and visualization - 0.0 fps");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		GLCapabilities caps = new GLCapabilities();
		caps.setDoubleBuffered(true);

		// initialize opengl Panel
		for(int i = 0; i < TEXTURE_COUNT; ++i) {
			textures[i] = -1;
		}
		panel = new GLCanvas(caps);
		panel.addGLEventListener(new MyGLEventListener());
		MouseListener ml = new MouseListener();
		panel.addMouseMotionListener(ml);
		panel.addMouseListener(ml);
		panel.addMouseWheelListener(ml);
		panel.addKeyListener(new MyKeyListener());
		panel.setFocusable(true);

		// add panel to window
		frame.setLayout(new BorderLayout());
		JComponent anotherpanel = initOptionPanel(frame);
		anotherpanel.setPreferredSize(new Dimension(256, 0));
		frame.add(anotherpanel, BorderLayout.EAST);
		frame.add(panel, BorderLayout.CENTER);
		frame.setSize(1000, 750);
		frame.doLayout();

		Point     p = new Point();
		Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
		frame.setLocation((int)(d.getWidth() / 2) - (frame.getWidth() / 2), (int)(d.getHeight() / 2) - (frame.getHeight() / 2));



		// show window
		frame.setVisible(true);
	}

	class MyKeyListener extends KeyAdapter {
		public void keyTyped(KeyEvent e) {
			char key = e.getKeyChar();
			Smoke.this.keyboard(key, 0, 0);
		}
	}

	class MouseListener implements MouseInputListener, MouseWheelListener {
		private int mouse_button_state = 0;

		public void mousePressed(MouseEvent e) {
			if(e.getButton() == MouseEvent.BUTTON1) {
				mouse_button_state |= 1;
			}
			if(e.getButton() == MouseEvent.BUTTON2) {
				mouse_button_state |= 2;
			}
			if(e.getButton() == MouseEvent.BUTTON3) {
				mouse_button_state |= 4;
			}
			mousex = e.getX();
			mousex_prev = mousex;
			mousey = e.getY();
			mousey_prev = mousey;
		}

		public void mouseClicked(MouseEvent e) {}
		public void mouseEntered(MouseEvent e) {}
		public void mouseExited(MouseEvent e) {}
		public void mouseMoved(MouseEvent e) {}

		public void mouseWheelMoved(MouseWheelEvent e) {
			int notches = e.getWheelRotation();
			if (notches < 0) { // Mouse wheel moved UP
				scale -= 1.0;
			}
			else { // Mouse wheel moved DOWN
				scale += 1.0;
			}
		}

		public void mouseReleased(MouseEvent e) {
			if(e.getButton() == MouseEvent.BUTTON1) {
				mouse_button_state &= ~1;
			}
			if(e.getButton() == MouseEvent.BUTTON2) {
				mouse_button_state &= ~2;
				rotationx = 0.0;
				rotationy = 0.0;
				scale     = 0.0;
			}
			if(e.getButton() == MouseEvent.BUTTON3) {
				mouse_button_state &= ~4;
			}
		}

		public void mouseDragged(MouseEvent e) {
			if((mouse_button_state&1)!=0) {
				drag(e.getX(), e.getY());
			}
			if((mouse_button_state&4)!=0) {
				mousex = e.getX();
				mousey = e.getY();
				rotationx += ((mousex - mousex_prev) / 10.0) % 360.0;
				rotationy += ((mousey - mousey_prev) / 10.0) % 360.0;
				mousex_prev = mousex;
				mousey_prev = mousey;
			}
		}
	}

	int nextPowerOfTwo(int arg) {
		int x = 1;
		while (x < arg) x <<= 1;
		return x;
	}

	int createTextureFromBuffer(GL gl, Buffer buffer, int internal_format, int pixel_format, int pixel_type, int dimensions, int width, int height) {
		int[] t = new int[1];
		int old_texid = -1;

		gl.glGenTextures(1, t, 0);
		int texid = t[0];
		gl.glBindTexture(dimensions, texid);
		gl.glTexParameteri(dimensions, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST);
		gl.glTexParameteri(dimensions, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST);
		gl.glTexParameteri(dimensions, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP);
		gl.glTexParameteri(dimensions, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP);
		buffer.position(0);
		if (dimensions == gl.GL_TEXTURE_1D) {
			gl.glTexImage1D(dimensions, 0, internal_format, width*height, 0, pixel_format, pixel_type, buffer);
		} else {
			gl.glTexImage2D(dimensions, 0, internal_format, width, height, 0, pixel_format, pixel_type, buffer);
		}
		return texid;
	}

	class MyGLEventListener implements GLEventListener {
		public void init(GLAutoDrawable drawable) {
			GL gl = drawable.getGL();
			gl.setSwapInterval(1); //Meh seems NOP in linux :(
// 			gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE_MINUS_SRC_ALPHA);
			gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA);
			gl.glHint(gl.GL_PERSPECTIVE_CORRECTION_HINT, gl.GL_NICEST);
			gl.glEnable(gl.GL_DEPTH_TEST);
			gl.glDepthFunc(gl.GL_LEQUAL);
			gl.glDepthMask(true);
			gl.glEnable(gl.GL_NORMALIZE);
			textures[TEXTURE_ARROW_1] = loadTexture(gl, "arrow-1.png");
			textures[TEXTURE_ARROW_2] = loadTexture(gl, "arrow-2.png");
			textures[TEXTURE_ARROW_3] = loadTexture(gl, "arrow-3.png");
		}

		public int loadTexture(GL gl, String fileName) {
			/* Load static textures */
			int texture = 0;
			try {
				BufferedImage image = null;
				image = (BufferedImage)ImageIO.read(new File(fileName));
				int[]  iarray = image.getRGB(0, 0, image.getWidth(), image.getHeight(), null, 0, image.getWidth());
				for (int i = 0; i < image.getWidth() * image.getHeight() ; ++i) {
					int c = iarray[i];
					int r = ((c >> 16) & 0xff);
					int g = ((c >>  8) & 0xff);
					int b = ((c >>  0) & 0xff);
					int a = ((c >> 24) & 0xff);//*/
					iarray[i] = (r << 0) | (g << 8) | (b << 16) | (a << 24);
				}
				IntBuffer buffer = IntBuffer.wrap(iarray);
				texture = createTextureFromBuffer(gl, buffer, gl.GL_RGBA, gl.GL_RGBA, gl.GL_UNSIGNED_BYTE, gl.GL_TEXTURE_2D, image.getWidth(), image.getHeight());
				gl.glTexParameteri(gl.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR);
				gl.glTexParameteri(gl.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR);
				return texture;
			} catch (Exception ex) {
				System.out.println("Error while loading textures:");
				System.out.println(ex.getMessage() + " / " + ex.toString());
				//ex.printStackTrace();
				System.exit(1);
			}
			return texture;
		}

		private void regenerate_colormap_texture(GL gl, ColormapSelectPanel panel, int PANEL_ID) {
			int ncolors = panel.getColorCount();
			gl.glEnable(gl.GL_TEXTURE_1D);
			gl.glColor3f(1.0f, 1.0f, 1.0f);
			int n = nextPowerOfTwo(ncolors + 1);
			FloatBuffer texture_data = BufferUtil.newFloatBuffer(n * 3);
			boolean fixup = ncolors + 1 != n;
			for (int i = 0 ; i <= ncolors; ++i) {
				double pos = (double)i / (double)ncolors;
				texture_data.put(panel.getGradientColor(pos));
				if (i == ncolors && fixup) { //Fixup clamping of colors to texture
					fixup = false;
					--i;
				}
			}
			textures[PANEL_ID] = createTextureFromBuffer(gl, texture_data, gl.GL_RGB, gl.GL_RGB, gl.GL_FLOAT, gl.GL_TEXTURE_1D, n, 1);
			texture_fill[PANEL_ID] = (float)((ncolors + 1) / (double)nextPowerOfTwo(ncolors + 1));
			panel.setUpdateGradientTexture(false);
		}

		public void display(GLAutoDrawable drawable) {
			Smoke.this.do_one_simulation_step();
			GL gl = drawable.getGL();

			ColormapSelectPanel[] colormappanels = {smokeColormapSelectPanel, isoLineSelectPanel};
			for(int i = TEXTURE_COLORMAP_SMOKE; i < TEXTURE_COLORMAP_ISOLINES + 1; ++i) {
				if(colormappanels[i].getUpdateGradientTexture()) {
					regenerate_colormap_texture(gl, colormappanels[i - TEXTURE_COLORMAP_SMOKE], i);
				}
			}

			Smoke.this.display(gl);
			gl.glFlush();
		}

		public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height) {
			Smoke.this.reshape(drawable.getGL(), width, height);
		}

		public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged) { }
	}
}
