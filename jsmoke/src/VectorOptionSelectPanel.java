import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.*;
import javax.swing.event.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.color.*;

import java.util.*;

public class VectorOptionSelectPanel extends ColormapSelectPanel {
// 	private float vector_size     = 24.0f;
// 	private float vector_scale_factor =  1.0f;
	private JSlider vector_size;
	private JLabel  vector_size_label;
	private JSlider vector_scale_factor;
	private JLabel  vector_scale_factor_label;
	private JSlider vector_grid_x;
	private JLabel  vector_grid_x_label;
	private JSlider vector_grid_y;
	private JLabel  vector_grid_y_label;
	public VectorOptionSelectPanel(int minColor, int maxColor, int colorCount, int colormap, JFrame frame) {
		super(minColor, maxColor, colorCount, colormap, frame);

		JPanel optionPanel = new JPanel();
		optionPanel.setBorder(new TitledBorder("Vector Options"));
		optionPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		vector_size         = new JSlider(1, 100, 24);
		vector_scale_factor = new JSlider(1, 1000, 500);
		vector_grid_x       = new JSlider(1, 256, 25);
		vector_grid_y       = new JSlider(1, 256, 25);
		vector_size.addChangeListener(new VectorOptionSelectPanelListener());
		vector_scale_factor.addChangeListener(new VectorOptionSelectPanelListener());
		vector_grid_x.addChangeListener(new VectorOptionSelectPanelListener());
		vector_grid_y.addChangeListener(new VectorOptionSelectPanelListener());
		optionPanel.setLayout(new SpringLayout());

		vector_size_label = new JLabel("Size [" + getVectorSize() + "]:");
		optionPanel.add(vector_size_label);
		optionPanel.add(vector_size);

		vector_scale_factor_label = new JLabel("Scale factor [" + getVectorScaleFactor() + "]:");
		optionPanel.add(vector_scale_factor_label);
		optionPanel.add(vector_scale_factor);

		vector_grid_x_label = new JLabel("Grid X [" + getVectorGridX() + "]:");
		optionPanel.add(vector_grid_x_label);
		optionPanel.add(vector_grid_x);

		vector_grid_y_label = new JLabel("Grid Y [" + getVectorGridY() + "]:");
		optionPanel.add(vector_grid_y_label);
		optionPanel.add(vector_grid_y);
		SpringUtilities.makeCompactGrid(optionPanel, 4, 2,  // rows, cols
		                                             6, 6,  // initX, initY
		                                             6, 6); // xPad, yPad

		add(optionPanel);
	}

	public float getVectorSize() {
		return vector_size.getValue();
	}

	public float getVectorScaleFactor() {
		float f = vector_scale_factor.getValue();
		if(f>=500) {
			f = (f - 400.0f) / 100.0f;
		}
		else {
			f = (float)Math.sqrt(1.0f / (500.0f - f));
		}
		return f;
	}

	public int getVectorGridX() {
			return vector_grid_x.getValue();
	}

	public int getVectorGridY() {
		return vector_grid_y.getValue();
	}

	class VectorOptionSelectPanelListener implements ChangeListener {
		public void stateChanged(ChangeEvent e) {
			if(e.getSource() == vector_size) {
				vector_size_label.setText("Size [" + getVectorSize() + "]:");
			}
			else if(e.getSource() == vector_scale_factor) {
				float f = ((int)(getVectorScaleFactor() * 100.0f))/100.0f;
				vector_scale_factor_label.setText("Scale factor [" + f + "]:");
			}
			else if(e.getSource() == vector_grid_x) {
				vector_grid_x_label.setText("Grid X [" + getVectorGridX() + "]:");
			}
			else if(e.getSource() == vector_grid_y) {
				vector_grid_y_label.setText("Grid Y [" + getVectorGridY() + "]:");
			}
		}
	}
}
