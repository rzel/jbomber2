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
        boolean draw_iso_lines = false;    //draw the iso lines or not
	static final int VECTOR_TYPE_HEDGEHOG = 0;
	static final int VECTOR_TYPE_ARROW    = VECTOR_TYPE_HEDGEHOG + 1;
	int vector_type = VECTOR_TYPE_ARROW;
	boolean frozen = false;         //toggles on/off the animation
	private static int TEXTURE_COLORMAP = 0;
	private static int TEXTURE_ARROW_1  = 1;
	private static int TEXTURE_ARROW_2  = 2;
	private static int TEXTURE_ARROW_3  = 3;
	private static int TEXTURE_ISOLINES = 4;
	int[] textures = new int[5];

	private ColormapSelectPanel smokeColormapSelectPanel;
	private VectorOptionSelectPanel vectorOptionSelectPanel;
        private IsoLineSelectPanel isoLineSelectPanel;
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
                                double div = 0.5* ((fx[getIdxFromXY(pos[0] - 1, pos[1])] - fx[getIdxFromXY(pos[0] + 1, pos[1])]) +
                                                   (fy[getIdxFromXY(pos[0], pos[1] - 1)] - fy[getIdxFromXY(pos[0], pos[1] + 1)]));
                                dataset_value = (float)div;
                                } break;
			case ColormapSelectPanel.DATASET_V_DIV: {
                                int[] pos = getXYFromIdx(idx);
                                double div = 0.5* ((vx[getIdxFromXY(pos[0] - 1, pos[1])] - vx[getIdxFromXY(pos[0] + 1, pos[1])]) +
                                                   (vy[getIdxFromXY(pos[0], pos[1] - 1)] - vy[getIdxFromXY(pos[0], pos[1] + 1)]));
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

	private double[][][] directional_vectors = new double[DIM][DIM][2];
	private void calculate_directional_vectors(ColormapSelectPanel panel) {
		int[][] indices = {
		                  {-1,-1},
		                  { 0,-1},
		                  { 1,-1},
		                  {-1, 0},
		                  { 1, 0},
		                  {-1, 1},
		                  { 0, 1},
		                  { 1, 1},
		                  };
		double[][] vectors = {
		                     {-Math.sqrt(2), Math.sqrt(2)},
		                     {0.0f, 1.0f},
                                     {Math.sqrt(2), Math.sqrt(2)},
                                     {1.0f, 0.0f},
				     {Math.sqrt(2), -Math.sqrt(2)},
				     {0.0f, -1.0f},
				     {-Math.sqrt(2), -Math.sqrt(2)},
				     {-1.0f, 0.0f},
				     };

		for(int x = 0 ; x < DIM; ++x) {
			for(int y = 0 ; y < DIM; ++y) {
				double dx = 0.0;
				double dy = 0.0;
				double centre = getDatasetColor(x+y*DIM, panel);
				for(int i = 0; i<8; ++i) {
					int idx = (((x + indices[i][0])+DIM) % DIM) + (((y + indices[i][1])+DIM) % DIM) * DIM;
					dx += Math.abs(getDatasetColor(idx, panel)-centre) * vectors[i][0];
					dy += Math.abs(getDatasetColor(idx, panel)-centre) * vectors[i][1];
				}
				directional_vectors[x][y][0] = dx;
				directional_vectors[x][y][1] = dy;
			}
		}
	}

	private double[] sampleDataset(int x, int y, VectorOptionSelectPanel panel) {
		// x and y here represent the (x,y) id of the grid-rectangle we're sampling
		double gridx = vectorOptionSelectPanel.getVectorGridX();
		double gridy = vectorOptionSelectPanel.getVectorGridY();
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
		for (double sample_x = -(cells_per_sample_x / 2.0); sample_x < cells_per_sample_x / 2.0; sample_x += 1.0) {
			for (double sample_y = -(cells_per_sample_y / 2.0); sample_y < cells_per_sample_y / 2.0; sample_y += 1.0) {
				double px, py;
				px = x_sample_centre_pos + sample_x;
				py = y_sample_centre_pos + sample_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell1s  = getDatasetColor((int)(px)+(int)(py)*DIM, vectorOptionSelectPanel);
				double cell1vx = vx[(int)(px)+(int)(py)*DIM];
				double cell1vy = vy[(int)(px)+(int)(py)*DIM];
				px = x_sample_centre_pos + sample_x + nearest_neighbour_x;
				py = y_sample_centre_pos + sample_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell2s  = getDatasetColor((int)(px)+(int)(py)*DIM, vectorOptionSelectPanel);
				double cell2vx = vx[(int)(px)+(int)(py)*DIM];
				double cell2vy = vy[(int)(px)+(int)(py)*DIM];
				px = x_sample_centre_pos + sample_x;
				py = y_sample_centre_pos + sample_y + nearest_neighbour_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell3s  = getDatasetColor((int)(px)+(int)(py)*DIM, vectorOptionSelectPanel);
				double cell3vx = vx[(int)(px)+(int)(py)*DIM];
				double cell3vy = vy[(int)(px)+(int)(py)*DIM];
				px = x_sample_centre_pos + sample_x + nearest_neighbour_x;
				py = y_sample_centre_pos + sample_y + nearest_neighbour_y;
				px = (px + DIM) % DIM;
				py = (py + DIM) % DIM;
				double cell4s  = getDatasetColor((int)(px)+(int)(py)*DIM, vectorOptionSelectPanel);
				double cell4vx = vx[(int)(px)+(int)(py)*DIM];
				double cell4vy = vy[(int)(px)+(int)(py)*DIM];
				avgs += (weight_sx * cell1s + weight_nnx * cell2s) * weight_sy + (weight_sx * cell3s + weight_nnx * cell4s) * weight_nny;
				avgx += (weight_sx * cell1vx + weight_nnx * cell2vx) * weight_sy + (weight_sx * cell3vx + weight_nnx * cell4vx) * weight_nny;
				avgy += (weight_sy * cell1vy + weight_nny * cell2vy) * weight_sy + (weight_sy * cell3vy + weight_nny * cell4vy) * weight_nny;
				++samples;
			}
		}

		double len = Math.sqrt(avgx*avgx + avgy*avgy);
		panel.update_longest_vector(len);
		double[] result = new double[4];
		result[0] = avgx;
		result[1] = avgy;
		result[2] = len;
		result[3] = avgs / ((double)samples);
		return result;
	}

	//visualize: This is the main visualization function
	void visualize(GL gl) {
		int        i, j, idx; double px, py;
		double/*fftw_real*/  wn = winWidth / (double)(DIM + 1);   // Grid cell width
		double/*fftw_real*/  hn = winHeight / (double)(DIM + 1);  // Grid cell heigh

		gl.glDisable(gl.GL_TEXTURE_2D);
		gl.glEnable(gl.GL_TEXTURE_1D);
		gl.glBindTexture(gl.GL_TEXTURE_1D, textures[TEXTURE_COLORMAP]);
		gl.glDisable(gl.GL_BLEND);

		if (draw_smoke) {
			gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
			for (j = 0; j < DIM - 1; j++) {		//draw smoke
				gl.glBegin(GL.GL_TRIANGLE_STRIP);

				i = 0;
				px = wn + (float)i * wn;
				py = hn + (float)j * hn;
				idx = (j * DIM) + i;
				gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
				gl.glVertex2d(px, py);

				for (i = 0; i < DIM - 1; i++) {
					px = wn + i * wn;
					py = hn + (j + 1) * hn;
					idx = ((j + 1) * DIM) + i;

					gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
					gl.glVertex2d(px, py);
					px = wn + (i + 1) * wn;
					py = hn + j * hn;
					idx = (j * DIM) + (i + 1);
					gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
					gl.glVertex2d(px, py);
				}

				px = wn + (float)(DIM - 1) * wn;
				py = hn + (float)(j + 1) * hn;
				idx = ((j + 1) * DIM) + (DIM - 1);
				gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
				gl.glVertex2d(px, py);
				gl.glEnd();
			}
		}

		if (draw_vecs) {
			//calculate_directional_vectors(vectorOptionSelectPanel);
			if (vector_type == VECTOR_TYPE_HEDGEHOG) {
				gl.glBegin(GL.GL_LINES);				//draw velocities
				for (i = 0; i < DIM; i++)
					for (j = 0; j < DIM; j++) {
						idx = (j * DIM) + i;
						//direction_to_color(gl, (float)(double)vx[idx],(float)(double)vy[idx],color_dir);
						//set_colormap(gl, getDatasetColor(idx, vectorOptionSelectPanel), vectorOptionSelectPanel);
						gl.glTexCoord1f(getDatasetColor(idx, vectorOptionSelectPanel) * texture_fill);
						gl.glVertex2d(wn + i * wn, hn + j * hn);
						gl.glVertex2d((wn + i * wn) + vec_scale * vx[idx], (hn + j * hn) + vec_scale * vy[idx]);
					}
				gl.glEnd();
			} else if (vector_type == VECTOR_TYPE_ARROW) {
				gl.glBindTexture(gl.GL_TEXTURE_2D, textures[3]);
				gl.glDisable(gl.GL_TEXTURE_1D);
				gl.glEnable(gl.GL_TEXTURE_2D);
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
				for (int x = 0 ; x < gridx ; ++x) {
					for (int y = 0 ; y < gridy ; ++y) {
						idx = (int)((x / (float)gridx) * DIM + DIM * (int)(DIM * (y / (float)gridy)));
						double[] result = sampleDataset(x, y, vectorOptionSelectPanel);
						double inprod = result[1]/result[2];
						double xdir   = result[0]/Math.abs(result[0]);
						double rotation = (-xdir)*(Math.acos(inprod)/Math.PI*180)+180;
						double size = 0.5 * (vector_size * vector_scalefactor * Math.sqrt(result[2]));
						if((vectorOptionSelectPanel.getScalemode() & vectorOptionSelectPanel.SCALE_SCALE) != 0) {
							size = 0.5 * (result[2] / maxveclen) * vecscalefact;
						}
						size = 0.5 * vecscalefact;
						float[] color = vectorOptionSelectPanel.getGradientColor(result[3]);
						gl.glColor3f(color[0], color[1], color[2]);
						gl.glPushMatrix();
						gl.glTranslatef((float)(spacex*0.5f + spacex * x - winWidth  * 0.5f + wn),
						                (float)(spacey*0.5f + spacey * y - winHeight * 0.5f + hn),
						                0.0f);
						gl.glRotatef((float)rotation, 0, 0, 1);
						gl.glBegin(GL.GL_QUADS); // Can not be moved outside of for-loop because of glTranslatef and glRotatef
						gl.glTexCoord2d(1.0, 0.0);
						gl.glVertex2d( - size, - size);

						gl.glTexCoord2d(0.0, 0.0);
						gl.glVertex2d( + size, - size);

						gl.glTexCoord2d(0.0, 1.0);
						gl.glVertex2d( + size, + size);

						gl.glTexCoord2d(1.0, 1.0);
						gl.glVertex2d( - size, + size);
						gl.glEnd();

						gl.glPopMatrix();
					}
				}
				gl.glPopMatrix();
			}
		}

		//boolean draw_isolines = true;
		if(draw_iso_lines) {
			gl.glEnable(gl.GL_BLEND);
			gl.glDisable(gl.GL_TEXTURE_2D);
			gl.glEnable(gl.GL_TEXTURE_1D);
			gl.glBindTexture(gl.GL_TEXTURE_1D, textures[TEXTURE_ISOLINES]);
			gl.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL);
			for (j = 0; j < DIM - 1; j++) {		//draw smoke
				gl.glBegin(GL.GL_TRIANGLE_STRIP);

				i = 0;
				px = wn + (float)i * wn;
				py = hn + (float)j * hn;
				idx = (j * DIM) + i;
				gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
				gl.glVertex2d(px, py);

				for (i = 0; i < DIM - 1; i++) {
					px = wn + i * wn;
					py = hn + (j + 1) * hn;
					idx = ((j + 1) * DIM) + i;

					gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
					gl.glVertex2d(px, py);
					px = wn + (i + 1) * wn;
					py = hn + j * hn;
					idx = (j * DIM) + (i + 1);
					gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
					gl.glVertex2d(px, py);
				}

				px = wn + (float)(DIM - 1) * wn;
				py = hn + (float)(j + 1) * hn;
				idx = ((j + 1) * DIM) + (DIM - 1);
				gl.glTexCoord1f(getDatasetColor(idx, smokeColormapSelectPanel) * texture_fill);
				gl.glVertex2d(px, py);
				gl.glEnd();
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
		glu.gluOrtho2D(0.0, (double)w, 0.0, (double)h);
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
			} else {
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

		smokeButton.setSelected(true);
		vectorButton.setSelected(false);
		draw_smoke = true;
		draw_vecs = false;
                draw_iso_lines = true;

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

		tabPane.setSelectedIndex(1);
		return tabPane;
	}

	GLCanvas panel;
	JFrame frame;
	public Smoke() {

		init_simulation(DIM);	//initialize the simulation data structures

		// initialize GUI
		frame = new JFrame("Real-time smoke simulation and visualization - 0.0 fps");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		GLCapabilities caps = new GLCapabilities();
		caps.setDoubleBuffered(true);

		// initialize opengl Panel
		panel = new GLCanvas(caps);
		panel.addGLEventListener(new MyGLEventListener());
		panel.addMouseMotionListener(new MouseListener());
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

	class MouseListener extends MouseMotionAdapter {
		public void mouseDragged(MouseEvent e) {
			drag(e.getX(), e.getY());
			//Flow.this.panel.repaint();
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

	float texture_fill = 1.0f;
	class MyGLEventListener implements GLEventListener {
		public void init(GLAutoDrawable drawable) {
			GL gl = drawable.getGL();
			gl.setSwapInterval(1); //Meh seems NOP in linux :(
// 			gl.glBlendFunc(gl.GL_ONE, gl.GL_ONE_MINUS_SRC_ALPHA);
			gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA);
			textures[1] = loadTexture(gl, "arrow-1.png");
			textures[2] = loadTexture(gl, "arrow-2.png");
			textures[3] = loadTexture(gl, "arrow-3.png");
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

		public void display(GLAutoDrawable drawable) {
			Smoke.this.do_one_simulation_step();
			GL gl = drawable.getGL();

			if (smokeColormapSelectPanel.getUpdateGradientTexture()) {
				int ncolors = smokeColormapSelectPanel.getColorCount();
				gl.glEnable(gl.GL_TEXTURE_1D);
				gl.glColor3f(1.0f, 1.0f, 1.0f);
				int n = nextPowerOfTwo(ncolors + 1);
				FloatBuffer texture_data = BufferUtil.newFloatBuffer(n * 3);
				boolean fixup = ncolors + 1 != n;
				for (int i = 0 ; i <= ncolors; ++i) {
					double pos = (double)i / (double)ncolors;
					texture_data.put(smokeColormapSelectPanel.getGradientColor(pos));
					if (i == ncolors && fixup) { //Fixup clamping of colors to texture
						fixup = false;
						--i;
					}
				}
				textures[0] = createTextureFromBuffer(gl, texture_data, gl.GL_RGB, gl.GL_RGB, gl.GL_FLOAT, gl.GL_TEXTURE_1D, n, 1);
				texture_fill = (float)((ncolors + 1) / (double)nextPowerOfTwo(ncolors + 1));
				smokeColormapSelectPanel.setUpdateGradientTexture(true);
			}

// 			if (isolineSelectPanel.getUpdateGradientTexture()) {
			if(isoLineSelectPanel.getUpdateGradientTexture()) {
				double iso_low_value  = Math.min(1.0, Math.max(0.0, isoLineSelectPanel.get_mindataset_value()));
				double iso_high_value = Math.min(1.0, Math.max(0.0, isoLineSelectPanel.get_maxdataset_value()));
				System.out.println("iso_low_value="+iso_low_value);
				System.out.println("iso_high_value="+iso_high_value);
				double wn = winWidth / (double)(DIM + 1);   // Grid cell width
				double hn = winHeight / (double)(DIM + 1);  // Grid cell heigh
				double nn = Math.max(wn,hn);
				int    iso_tex_size   = Math.min(nextPowerOfTwo((int)nn), 2048);
				int    iso_min_texels = (int)(iso_tex_size/nn/(20.0/9.0)+0.5);

				System.out.println("iso_min_texels="+iso_min_texels);
				System.out.println("nn="+nn);

				double iso_n_texels      = ((iso_high_value - iso_low_value) * iso_tex_size);
				double iso_n_low_texels  = (iso_low_value                    * iso_tex_size);
				double iso_n_high_texels = ((1.0-iso_high_value)             * iso_tex_size);

// 				System.out.println("iso_n_low_texels="+iso_n_low_texels);
// 				System.out.println("iso_n_high_texels="+iso_n_high_texels);
// 				System.out.println("iso_n_texels="+iso_n_texels);
// 				System.out.println("iso_tex_size="+iso_tex_size);

				if(iso_n_texels < iso_min_texels ) {
					System.out.println("Adjust");
					double  diff = iso_min_texels - iso_n_texels;
					iso_n_low_texels  -= (diff/2.0);
					iso_n_high_texels -= (diff/2.0);
					iso_n_texels = iso_min_texels;
				}

// 				System.out.println("iso_n_low_texels="+iso_n_low_texels);
// 				System.out.println("iso_n_high_texels="+iso_n_high_texels);
// 				System.out.println("iso_n_texels="+iso_n_texels);
// 				System.out.println("iso_tex_size="+iso_tex_size);

				if(iso_n_low_texels + iso_n_texels + iso_n_high_texels > iso_tex_size) {
					System.out.println("Error: texture size exceeded");
					return;
				}

				float[] intermediate_texbuf = new float[2048*4];
                                //float[] ic = isoLineSelectPanel.getIsoLineColor();
				//float[] fc = {ic[0] , ic[1] , ic[2] , 0.0f};
                                float[] fc = {1.0f, 0.0f, 1.0f, 0.0f};
				for(int i = 0 ; i < 2048*4; i+=4) {
					intermediate_texbuf[i+0] = fc[0];
					intermediate_texbuf[i+1] = fc[1];
					intermediate_texbuf[i+2] = fc[2];
					intermediate_texbuf[i+3] = fc[3];
				}
                                
                                fc = new float[] {1.0f, 0.0f, 1.0f, 0.0f};
				for(int i = 4*(int)(iso_n_low_texels+0.5); i < 4*(int)(iso_n_low_texels+iso_n_texels+0.5); i+=4) {
					intermediate_texbuf[i+0] = fc[0];
					intermediate_texbuf[i+1] = fc[1];
					intermediate_texbuf[i+2] = fc[2];
					intermediate_texbuf[i+3] = fc[3];
				}

				FloatBuffer texture_data = FloatBuffer.wrap(intermediate_texbuf);//BufferUtil.newFloatBuffer(iso_tex_size*4);
				textures[TEXTURE_ISOLINES] = createTextureFromBuffer(gl, texture_data, gl.GL_RGBA, gl.GL_RGBA, gl.GL_FLOAT, gl.GL_TEXTURE_1D, iso_tex_size, 1);
				updateisotexture = false;
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
