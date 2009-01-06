import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.*;
import javax.swing.event.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.color.*;

import java.util.*;

public class VectorOptionSelectPanel extends ColormapSelectPanel implements ItemListener{
	public static final int VECTOR_FIELD_FORCE = 42;
	public static final int VECTOR_FIELD_VELOCITY = 21;
    
	private static final String HEDGEHOGS = "1d: Hedgehogs";
	private static final String ARROW     = "2d: Arrow";
	private static final String ARROWTAIL = "2d: Arrow tail";
	private static final String ARROWTOP  = "2d: Arrow top";
	private static final String PACMAN    = "2d: Pacman";
	private static final String PYRAMID   = "3d: Pyramid";	
	
	private static final int VT_HEDGEHOGS = 42;
	private static final int VT_ARROW     = 84;
	private static final int VT_ARROWTAIL = 126;
	private static final int VT_ARROWTOP  = 168;
	private static final int VT_PACMAN    = 210;
	private static final int VT_PYRAMID   = 252;
	
	private JSlider vector_size;
	private JLabel  vector_size_label;
	private JSlider vector_scale_factor;
	private JLabel  vector_scale_factor_label;
	private JSlider vector_grid_x;
	private JLabel  vector_grid_x_label;
	private JSlider vector_grid_y;
	private JLabel  vector_grid_y_label;
	private double  longest_vector;
	private double  longest_vector_last;
	
	private JRadioButton forceButton = new JRadioButton("Force field (f)", true);
	private JRadioButton velocityButton = new JRadioButton("Fluid velocity (v)", false);
	
	private JComboBox shapeBox = new JComboBox();
	
	private int vectorType = VT_PACMAN;
	
	public VectorOptionSelectPanel(int minColor, int maxColor, int colorCount, int colormap, JFrame frame) {
		super(minColor, maxColor, colorCount, colormap, frame);

		longest_vector = 0.0;
		longest_vector_last = 1.0;
		
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


		/****************************
		 ** BEGIN Default gradient **
		 ****************************/
		DefaultTableModel tablemodelcolors = (DefaultTableModel)colortable.getModel();
		TableColumnModel cm = colortable.getColumnModel();
		ColorTableRenderer r = new ColorTableRenderer();
		TableColumn c = cm.getColumn(0);
		c.setCellRenderer(r);

		tablemodelcolors.setRowCount(0);
		Object[] o = new Object[1];
		o[0] = new Color(255, 255, 0);
		tablemodelcolors.addRow(o);
		o[0] = new Color(255, 0, 20);
		tablemodelcolors.addRow(o);
		colortable.changeSelection(colortable.getSelectedRow(), 0, false, false);
		colorCountSlider.setMinimum(colortable.getRowCount() - 1);	
		generate_custom_gradient_cache();
		/**************************
		 ** END Default gradient **
		 **************************/

		shapeBox.addItem(HEDGEHOGS);
		shapeBox.addItem(ARROW);
		shapeBox.addItem(ARROWTAIL);
		shapeBox.addItem(ARROWTOP);
		shapeBox.addItem(PACMAN);
		shapeBox.addItem(PYRAMID);
		optionPanel.add(new JLabel("Vector shape:"));
		optionPanel.add(shapeBox);		
		
		ButtonGroup vectorGroup = new ButtonGroup();	
		vectorGroup.add(forceButton);
		vectorGroup.add(velocityButton);
		optionPanel.add(forceButton);
		optionPanel.add(velocityButton);
		
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
		SpringUtilities.makeCompactGrid(optionPanel, 6, 2,  // rows, cols
		                                6, 6,  // initX, initY
		                                6, 6); // xPad, yPad

		add(optionPanel);
	}

	public int getVectorField() {
	    if (forceButton.getModel().isSelected()) {
		return VECTOR_FIELD_FORCE;
	    }
	    else {
		return VECTOR_FIELD_VELOCITY;
	    }
	}
	
	public void update_longest_vector(double v) {
		longest_vector = v > longest_vector ? v : longest_vector;
	}

	public void reset_longest_vector() {
		longest_vector_last = longest_vector;
		longest_vector = 0.0;
	}

	public double get_longest_vector() {
		return longest_vector_last;
	}

	public float getVectorSize() {
		return vector_size.getValue();
	}

	public float getVectorScaleFactor() {
		float f = vector_scale_factor.getValue();
		if (f >= 500) {
			f = (f - 400.0f) / 100.0f;
		} else {
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
	
	public int getVectorType() {
	    return vectorType;
	}

	public void itemStateChanged(ItemEvent e) {
	    if (e.equals(HEDGEHOGS)) {
		vectorType = VT_HEDGEHOGS;
	    }
	    else if (e.equals(ARROW)) {
		vectorType = VT_ARROW;
	    }
	    else if (e.equals(ARROWTAIL)) {
		vectorType = VT_ARROWTAIL;
	    }
	    else if (e.equals(ARROWTOP)) {
		vectorType = VT_ARROWTOP;
	    }
	    else if (e.equals(PACMAN)) {
		vectorType = VT_PACMAN;
	    }
	    else if (e.equals(PYRAMID)) {
		vectorType = VT_PYRAMID;
	    }
	}
	
	class VectorOptionSelectPanelListener implements ChangeListener {
		public void stateChanged(ChangeEvent e) {
			if (e.getSource() == vector_size) {
				vector_size_label.setText("Size [" + getVectorSize() + "]:");
			} else if (e.getSource() == vector_scale_factor) {
				float f = ((int)(getVectorScaleFactor() * 100.0f)) / 100.0f;
				vector_scale_factor_label.setText("Scale factor [" + f + "]:");
			} else if (e.getSource() == vector_grid_x) {
				vector_grid_x_label.setText("Grid X [" + getVectorGridX() + "]:");
			} else if (e.getSource() == vector_grid_y) {
				vector_grid_y_label.setText("Grid Y [" + getVectorGridY() + "]:");
			}
		}
	}
}
