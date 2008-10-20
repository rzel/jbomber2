import javax.swing.*;
import javax.swing.border.*;
import javax.swing.table.*;
import javax.swing.event.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.color.*;

import java.util.*;

public class ColormapSelectPanel extends JPanel implements ActionListener, ChangeListener {
    public static final int COLOR_RAINBOW   = 1;
    public static final int COLOR_GRAYSCALE = 2;
    public static final int COLOR_DEFINED   = 4;
    public static final int COLOR_CUSTOM    = 8;

    public static final int DATASET_RHO = 1;
    public static final int DATASET_F   = 2;
    public static final int DATASET_V   = 4;

    public static final int SCALE_CLAMP = 1;
    public static final int SCALE_SCALE = 2;

		private float maxdataset_value = Float.MIN_VALUE;
		private float mindataset_value = Float.MAX_VALUE;
		private float maxdataset_value_last = 1.0f;
		private float mindataset_value_last = 0.0f;

    private JSliderlessSlider colorPreviewSlider;
    private int colorCount = 2047;
    private ColorTable colortable;
    private ColorSelector colorselector;
    private JLabel colorCountLabel;
    private JSlider colorCountSlider;

    private int dataset   = DATASET_RHO;
    private int scalemode = SCALE_SCALE;
    private int colormap  = COLOR_CUSTOM;

    private boolean update_gradient_texture;

    private UnboundedDoubleSpinnerModel minClampSelectSpinnerModel = new UnboundedDoubleSpinnerModel(0.0);
    private UnboundedDoubleSpinnerModel maxClampSelectSpinnerModel = new UnboundedDoubleSpinnerModel(1.0);    private static float[][] custom_gradient_cache = new float[2048][3];
    int custom_gradient_interpolate_mode = 0;

    public ColormapSelectPanel(int minColor, int maxColor, int colorCount, int colormap, JFrame frame) {
        add(initDatasetSelectPanel());
        add(initScalingSelectPanel());
        add(initColormapPreviewPanel());
        add(initColormapSelectPanel(frame));

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        setAlignmentX(Component.LEFT_ALIGNMENT);

        updateColormap(minColor, maxColor, colorCount, colormap);
    }

		public void update_maxdataset_value(float dataset_value) {
			maxdataset_value = dataset_value > maxdataset_value ? dataset_value : maxdataset_value;
		}

		public void reset_maxdataset_value() {
			maxdataset_value_last = maxdataset_value;
			maxdataset_value = Float.MIN_VALUE;
		}

		public float get_maxdataset_value() {
			return maxdataset_value_last;
		}

		public void update_mindataset_value(float dataset_value) {
			mindataset_value = dataset_value < mindataset_value ? dataset_value : mindataset_value;
		}

		public void reset_mindataset_value() {
			mindataset_value_last = mindataset_value;
			mindataset_value = Float.MAX_VALUE;
		}

		public float get_mindataset_value() {
			return mindataset_value_last;
		}

    public int getDataset() {
         return dataset;
    }

    public int getScalemode() {
         return scalemode;
    }

    public double getMinClamp() {
         return ((Double)minClampSelectSpinnerModel.getValue()).doubleValue();
    }

    public double getMaxClamp() {
         return ((Double)maxClampSelectSpinnerModel.getValue()).doubleValue();
    }

    public boolean getUpdateGradientTexture() {
         return update_gradient_texture;
    }

    public void setUpdateGradientTexture(boolean value) {
         update_gradient_texture = value;
    }

    public int getColorCount() {
         return colorCount;
    }

    public int getColormap() {
         return colormap;
    }

    public float[] getCustomGradient(int index) {
         return custom_gradient_cache[index];
    }

		public void actionPerformed(ActionEvent e) {
			     if(e.getActionCommand() == "SET_DATASET_RHO") {
				dataset = DATASET_RHO;
			}
			else if(e.getActionCommand() == "SET_DATASET_F") {
				dataset = DATASET_F;
			}
			else if(e.getActionCommand() == "SET_DATASET_V") {
				dataset = DATASET_V;
			}
			else if(e.getActionCommand() == "SET_SCALE_CLAMP") {
				scalemode ^= SCALE_CLAMP;
			}
			else if(e.getActionCommand() == "SET_SCALE_SCALE") {
				scalemode ^= SCALE_SCALE;
			}
			else if(e.getActionCommand() == "SET_RAINBOW") {
				float[][] colors = new float[2048][3];

				for(int i = 0; i < colors.length; ++i) {
					rainbow(i / (float)colors.length, colors[i]);
				}
				colormap = COLOR_RAINBOW;
				colorPreviewSlider.setColors(colors, colorCount);
			}
			else if(e.getActionCommand() == "SET_GRAYSCALE") {
				float[][] colors = new float[2048][3];

				for(int i = 0; i < colors.length; ++i) {
					colors[i][0] = colors[i][1] = colors[i][2] = i / (float)colors.length;
				}
				colormap = COLOR_GRAYSCALE;
				colorPreviewSlider.setColors(colors, colorCount);
			}
			else if(e.getActionCommand() == "SET_DEFINED") {
				float[][] colors = new float[2048][3];
				final int NLEVELS = 7;

				for(int i = 0; i < colors.length; ++i) {
					float vy = (float)(i / (double)(colors.length - 1));
					vy *= NLEVELS; vy = (int)(vy); vy/= NLEVELS;
					rainbow(vy, colors[i]);
				}
				colormap = COLOR_DEFINED;
				colorPreviewSlider.setColors(colors, colorCount);
			}
			else if(e.getActionCommand() == "SET_CUSTOM") {
				colormap = COLOR_CUSTOM;
				colorPreviewSlider.setColors(custom_gradient_cache, colorCount);
			}
			else if (e.getActionCommand().equals("COLORTABLE_PICK_COLOR")) {
				int row = colortable.getSelectedRow();
				int column = colortable.getSelectedColumn();
				if(row>=0 && column>=0)
					colorselector.setSelectedColor((Color)colortable.getValueAt(row, column));
				colorselector.setVisible(true);
				generate_custom_gradient_cache();
				update_gradient_texture = true;
			}
			else if (e.getActionCommand().equals("COLORTABLE_ADD_NEW_COLOR")) {
				int row = colortable.getSelectedRow();
				Object[] o=new Object[1];
				o[0]=new Color((int)(Math.random()*256), (int)(Math.random()*256), (int)(Math.random()*256));
				((DefaultTableModel)(colortable.getModel())).insertRow(row, o);
				colortable.changeSelection(row, 0, false, false);
				colorCountSlider.setMinimum(colortable.getRowCount() - 1);
				generate_custom_gradient_cache();
				update_gradient_texture = true;
			}
			else if (e.getActionCommand().equals("COLORTABLE_REMOVE_COLOR")) {
				int row = colortable.getSelectedRow();
				((DefaultTableModel)(colortable.getModel())).removeRow(row);
				colorCountSlider.setMinimum(colortable.getRowCount() -1 );
				generate_custom_gradient_cache();
				update_gradient_texture = true;
			}
			else if (e.getActionCommand().equals("INTERPOLATE_RGB")) {
				custom_gradient_interpolate_mode = 0;
				generate_custom_gradient_cache();
				update_gradient_texture = true;
			}
			else if (e.getActionCommand().equals("INTERPOLATE_HSV")) {
				custom_gradient_interpolate_mode = 1;
				generate_custom_gradient_cache();
				update_gradient_texture = true;
			}
		}

		public void stateChanged(ChangeEvent e) {
			if(e.getSource() == colorCountSlider) {
				int value = ((JSlider)e.getSource()).getValue();
				colorCountLabel.setText("Limit colors to " + (value + 1));
				colorCount = value;
				colorPreviewSlider.setCount(colorCount);
				generate_custom_gradient_cache();
				update_gradient_texture = true;
			}
			else if(e.getSource() ==colorselector ) {
				int row = colortable.getSelectedRow();
				int column = colortable.getSelectedColumn();
				if(row>=0 && column>=0) {
					colortable.setValueAt(e.getSource(), row, column);
				}
				generate_custom_gradient_cache();
				update_gradient_texture = true;
			}
		}

    private JPanel initDatasetSelectPanel() {
        JRadioButton rhoButton = new JRadioButton("rho");
        rhoButton.setSelected(true);
        rhoButton.addActionListener(this);
				rhoButton.setActionCommand("SET_DATASET_RHO");

        JRadioButton fButton = new JRadioButton("|f|");
        fButton.addActionListener(this);
				fButton.setActionCommand("SET_DATASET_F");

				JRadioButton vButton = new JRadioButton("|v|");
				vButton.addActionListener(this);
				vButton.setActionCommand("SET_DATASET_V");

        ButtonGroup datasetSelectGroup = new ButtonGroup();
        datasetSelectGroup.add(rhoButton);
        datasetSelectGroup.add(fButton);
        datasetSelectGroup.add(vButton);

        JPanel datasetSelectPanel = new JPanel();
        datasetSelectPanel.setLayout(new BoxLayout(datasetSelectPanel, BoxLayout.X_AXIS));
        datasetSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        datasetSelectPanel.setBorder(new TitledBorder("Dataset"));
        datasetSelectPanel.add(rhoButton);
        datasetSelectPanel.add(fButton);
        datasetSelectPanel.add(vButton);
        return datasetSelectPanel;
    }

    private JPanel initScalingSelectPanel() {
        JCheckBox clampButton = new JCheckBox("clamping");
        clampButton.setMnemonic(KeyEvent.VK_C);
        clampButton.setActionCommand("SET_SCALE_CLAMP");
        clampButton.addActionListener(this);

        JCheckBox scaleButton = new JCheckBox("scaling");
        scaleButton.setMnemonic(KeyEvent.VK_S);
        scaleButton.setActionCommand("SET_SCALE_SCALE");
        scaleButton.setSelected(true);
        scaleButton.addActionListener(this);

        colorCountLabel = new JLabel("Limit colors to 2047");
        colorCountSlider = new JSlider(JSlider.HORIZONTAL, 1, 2047, 2047);
        colorCountSlider.addChangeListener(this);

        JPanel clampSelectPanel = new JPanel();
        clampSelectPanel.setLayout(new GridLayout(3,2));
        clampSelectPanel.add(clampButton);
        clampSelectPanel.add(new JLabel());
        clampSelectPanel.add(new JLabel("Low: ", SwingConstants.RIGHT));
        clampSelectPanel.add(new JSpinner(minClampSelectSpinnerModel));
        clampSelectPanel.add(new JLabel("High: ", SwingConstants.RIGHT));
        clampSelectPanel.add(new JSpinner(maxClampSelectSpinnerModel));

        JPanel scaleSelectPanel = new JPanel();
        scaleSelectPanel.setLayout(new BoxLayout(scaleSelectPanel, BoxLayout.Y_AXIS));
        scaleButton.setAlignmentX(Component.LEFT_ALIGNMENT);
        clampSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        colorCountLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
        colorCountSlider.setAlignmentX(Component.LEFT_ALIGNMENT);
        scaleSelectPanel.setBorder(new TitledBorder("Scaling:"));
        scaleSelectPanel.add(scaleButton);
        scaleSelectPanel.add(clampSelectPanel);
        scaleSelectPanel.add(colorCountLabel);
        scaleSelectPanel.add(colorCountSlider);
        return scaleSelectPanel;
    }

    private JPanel initColormapPreviewPanel() {
        colorPreviewSlider = new JSliderlessSlider(new DefaultBoundedRangeModel());
        colorPreviewSlider.setMajorTickSpacing(1);
        colorPreviewSlider.setPaintTicks(true);
        colorPreviewSlider.setPaintLabels(true);

        JPanel colormapPreviewPanel = new JPanel();
        colormapPreviewPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        colormapPreviewPanel.setLayout(new BoxLayout(colormapPreviewPanel, BoxLayout.Y_AXIS));
        colormapPreviewPanel.setBorder(new TitledBorder("Colormap preview:"));
        colormapPreviewPanel.add(colorPreviewSlider);
        return colormapPreviewPanel;
    }

    private JPanel initColormapSelectPanel(JFrame frame) {
       // Initialize colormap selection
        JRadioButton rainbowButton = new JRadioButton("Rainbow");
        rainbowButton.addActionListener(this);
				rainbowButton.setActionCommand("SET_RAINBOW");

        JRadioButton grayscaleButton = new JRadioButton("Grayscale");
        grayscaleButton.addActionListener(this);
				grayscaleButton.setActionCommand("SET_GRAYSCALE");

        JRadioButton definedButton   = new JRadioButton("Defined");
        definedButton.addActionListener(this);
				definedButton.setActionCommand("SET_DEFINED");

        JRadioButton customButton   = new JRadioButton("Custom");
        customButton.setSelected(true);
        customButton.addActionListener(this);
				customButton.setActionCommand("SET_CUSTOM");

        ButtonGroup colorMapSelectGroup = new ButtonGroup();
        colorMapSelectGroup.add(rainbowButton);
        colorMapSelectGroup.add(grayscaleButton);
        colorMapSelectGroup.add(definedButton);
        colorMapSelectGroup.add(customButton);

        rainbowButton.setAlignmentX(Component.LEFT_ALIGNMENT);
        grayscaleButton.setAlignmentX(Component.LEFT_ALIGNMENT);
        definedButton.setAlignmentX(Component.LEFT_ALIGNMENT);
        customButton.setAlignmentX(Component.LEFT_ALIGNMENT);

        JPanel colormapSelectPanel = new JPanel();
        colormapSelectPanel.setLayout(new BoxLayout(colormapSelectPanel, BoxLayout.Y_AXIS));
        colormapSelectPanel.setBorder(new TitledBorder("Colormaps"));

        colormapSelectPanel.add(rainbowButton);
        colormapSelectPanel.add(grayscaleButton);
        colormapSelectPanel.add(definedButton);
        colormapSelectPanel.add(customButton);

        JPanel y = initCustomColorPanel(frame);
        y.setBorder(new EmptyBorder(0,24,0,0));
        y.setAlignmentX(Component.LEFT_ALIGNMENT);
        colormapSelectPanel.add(y);

        colormapSelectPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        return colormapSelectPanel;
    }

    public void updateColormap(int minColor,
                               int maxColor,
                               int colorCount,
                               int colormap) {
        colorPreviewSlider.setMinimum(minColor);
        colorPreviewSlider.setMaximum(maxColor);
        this.colorCount = colorCount;
        this.colormap = colormap;
        setColors(colormap, colorCount);
    }

    private void setColors(int colormap, int colorCount) {
        float[][] colors = new float[2048][3];

        switch (colormap) {
            case COLOR_RAINBOW:
                 for(int i = 0; i < colors.length; ++i) {
                    rainbow(i / (float)colors.length, colors[i]);
                }
                break;
            case COLOR_GRAYSCALE:
                for(int i = 0; i < colors.length; ++i) {
                    colors[i][0] = colors[i][1] = colors[i][2] = i / (float)colors.length;
                }
                break;
            case COLOR_DEFINED:
                final int NLEVELS = 7;

                for(int i = 0; i < colors.length; ++i) {
                    float vy = (float)(i / (double)(colors.length - 1));
                    vy *= NLEVELS; vy = (int)(vy); vy/= NLEVELS;
                    rainbow(vy, colors[i]);
                }

                break;
            case COLOR_CUSTOM:
                colors = custom_gradient_cache;
                break;
        }
        colorPreviewSlider.setColors(colors, colorCount);
    }

    public static void rainbow(float value,float[] color) {
        final float dx=0.8f;
        if (value<0) value=0; if (value>1) value=1;
        value = (6-2*dx)*value+dx;
        color[0] = Math.max(0.0f,(3-Math.abs(value-4)-Math.abs(value-5))/2);
        color[1] = Math.max(0.0f,(4-Math.abs(value-2)-Math.abs(value-4))/2);
        color[2] = Math.max(0.0f,(3-Math.abs(value-1)-Math.abs(value-2))/2);
    }

    private void generate_custom_gradient_cache() {
        int n = colortable.getRowCount();
        ColorSpace cs = ((Color)colortable.getValueAt(0, 0)).getColorSpace(); // Assume all colors use the same color space
        Color a = null, b = null;
        float[] a_cc, b_cc;
        int cm = -1;
        for(int k = 0; k < 2048; ++k) {
            float f = (float)(k * 1.0/2048.0);
            float r = (n-1) * f;
            int m = (int)r;

            if(m != cm) {
                a = (Color)colortable.getValueAt(m, 0);
                b = (Color)colortable.getValueAt(Math.min(m+1, n-1), 0);
                cm = m;
            }

            float h = r - m;
            float l = 1 - h;

            switch(custom_gradient_interpolate_mode) {
                case 0:
                    custom_gradient_cache[k][0] = (a.getRed()   * l + b.getRed()   * h) / 255.0f;
                    custom_gradient_cache[k][1] = (a.getGreen() * l + b.getGreen() * h) / 255.0f;
                    custom_gradient_cache[k][2] = (a.getBlue()  * l + b.getBlue()  * h) / 255.0f;
                    break;
                case 1:
                    a_cc = a.getColorComponents(null);
                    b_cc = b.getColorComponents(null);

                    a_cc = cs.toCIEXYZ(a_cc);
                    b_cc = cs.toCIEXYZ(b_cc);
                    a_cc[0] = (a_cc[0] * l + b_cc[0] * h);
                    a_cc[1] = (a_cc[1] * l + b_cc[1] * h);
                    a_cc[2] = (a_cc[2] * l + b_cc[2] * h);
                    a_cc = cs.fromCIEXYZ(a_cc);

                    float[] rgb = cs.toRGB(a_cc);
                    custom_gradient_cache[k][0] = a_cc[0];
                    custom_gradient_cache[k][1] = a_cc[1];
                    custom_gradient_cache[k][2] = a_cc[2];
                    break;
            }
        }
    }

    private JPanel initCustomColorPanel(JFrame frame) {
        colorselector = new ColorSelector(frame);
        colorselector.addChangeListener(this);

        DefaultTableModel tablemodelcolors= new DefaultTableModel(1,1);

        tablemodelcolors.addTableModelListener(new ColorTableTableModelListener());
        colortable = new ColorTable(tablemodelcolors);
        colortable.addMouseListener(new ColorTableRightClick(this));

        colortable.getColumnModel().getColumn(colortable.convertColumnIndexToView(0)).setHeaderValue("Colors");
        colortable.doLayout();

        TableColumnModel cm = colortable.getColumnModel();
        ColorTableRenderer r = new ColorTableRenderer();
        TableColumn c = cm.getColumn(0);
        c.setCellRenderer(r);

        tablemodelcolors.setRowCount(0);
        Object[] o = new Object[1];
        o[0] = new Color(0,0,0);
        tablemodelcolors.addRow(o);
        o[0] = new Color(0,255,255);
        tablemodelcolors.addRow(o);
        o[0] = new Color(255,0,255);
        tablemodelcolors.addRow(o);
        colorCountSlider.setMinimum(colortable.getRowCount() - 1 );
        generate_custom_gradient_cache();

        JScrollPane scroll  = new JScrollPane(colortable);
        JPanel panel = new JPanel();
        panel.setBorder(new TitledBorder("Current gradient"));
        BoxLayout panellayout = new BoxLayout(panel, BoxLayout.X_AXIS);
        panel.setLayout(panellayout);

        JRadioButton colorInterpolateModeRGB = new JRadioButton("Interpolate in RGB space");
        JRadioButton colorInterpolateModeHSV = new JRadioButton("Interpolate in HSV space");
        colorInterpolateModeRGB.setActionCommand("INTERPOLATE_RGB");
        colorInterpolateModeHSV.setActionCommand("INTERPOLATE_HSV");
        colorInterpolateModeRGB.addActionListener(this);
        colorInterpolateModeHSV.addActionListener(this);
        colorInterpolateModeRGB.setSelected(custom_gradient_interpolate_mode == 0);
        colorInterpolateModeHSV.setSelected(custom_gradient_interpolate_mode == 0);
        ButtonGroup colorInterpolateModeGroup = new ButtonGroup();
        colorInterpolateModeGroup.add(colorInterpolateModeRGB);
        colorInterpolateModeGroup.add(colorInterpolateModeHSV);
        JPanel colorInterpolateModePanel = new JPanel();
        colorInterpolateModePanel.setLayout(new BoxLayout(colorInterpolateModePanel, BoxLayout.Y_AXIS));
        colorInterpolateModeRGB.setAlignmentX(Component.LEFT_ALIGNMENT);
        colorInterpolateModeHSV.setAlignmentX(Component.LEFT_ALIGNMENT);
        colorInterpolateModePanel.add(colorInterpolateModeRGB);
        colorInterpolateModePanel.add(colorInterpolateModeHSV);
        colorInterpolateModePanel.setAlignmentX(Component.LEFT_ALIGNMENT);
        scroll.setAlignmentX(Component.LEFT_ALIGNMENT);
        panel.add(colorInterpolateModePanel);
        panel.add(scroll);
        panel.setLayout(new BoxLayout(panel, BoxLayout.Y_AXIS));
        return panel;
    }
}

class UnboundedDoubleSpinnerModel implements SpinnerModel {
    private Vector listeners = new Vector();

    private double value;

    public UnboundedDoubleSpinnerModel(double value) {
        this.value = value;
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
    }

    public Object getValue() {
        return new Double(value);
    }
}